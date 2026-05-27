#!/usr/bin/env python3
"""
QC Assessment Script
Evaluates fastp results and mapping rates to determine which samples pass QC.
"""

import json
import re
import os
from pathlib import Path


def parse_fastp_data(fastp_files):
    """Parse fastp JSON output to extract quality metrics."""
    metrics = {}

    for fastp_json in fastp_files:
        sample = Path(fastp_json).stem.replace("_fastp", "")
        with open(fastp_json) as f:
            data = json.load(f)

        summary = data.get("summary", {}).get("after_filtering", {})
        total_reads = summary.get("total_reads", 0)
        total_bases = summary.get("total_bases", 0)
        average_length = summary.get("mean_length") or (
            total_bases / total_reads if total_reads else 0
        )

        q30_rate = summary.get("q30_rate")
        if q30_rate is None:
            q30_bases = summary.get("q30_bases", 0)
            q30_rate = q30_bases / total_bases if total_bases else 0

        gc_content = summary.get("gc_content", 0)
        if isinstance(gc_content, (int, float)) and gc_content <= 1:
            gc_content = gc_content * 100

        metrics[sample] = {
            "Total reads": total_reads,
            "Average length": average_length,
            "Q30%": q30_rate * 100 if q30_rate is not None else 0,
            "GC content %": gc_content,
            "Passed reads": data.get("filtering_result", {}).get("passed_filter_reads", 0),
        }

    return metrics


def parse_flagstat_data(flagstat_files):
    """Parse samtools flagstat output to extract mapping rates."""
    mapping = {}
    pattern = re.compile(r"(\d+) \+ \d+ mapped \((\d+\.\d+)%")

    for flagstat_file in flagstat_files:
        sample = Path(flagstat_file).stem.replace(".flagstat", "")
        rate = 0.0
        try:
            with open(flagstat_file) as f:
                for line in f:
                    m = pattern.search(line)
                    if m:
                        rate = float(m.group(2))
                        break
        except FileNotFoundError:
            pass
        mapping[sample] = rate

    return mapping


def assess_qc(metrics, mapping_rates, thresholds):
    """Assess samples against QC thresholds."""
    results = {}

    for sample, data in metrics.items():
        status = "PASS"
        issues = []

        length = float(data.get("Average length", 0))
        if length < thresholds.get("min_read_length", 50):
            status = "FAIL"
            issues.append(f"Low read length: {length:.1f} bp")

        q30 = float(data.get("Q30%", 0))
        if q30 < thresholds.get("min_q30_pct", 75):
            status = "FAIL"
            issues.append(f"Low Q30 rate: {q30:.1f}%")

        gc = float(data.get("GC content %", 0))
        if gc < 30 or gc > 50:
            issues.append(f"Unusual GC content: {gc:.1f}% (expected ~38%)")

        if data.get("Passed reads", 0) == 0:
            status = "FAIL"
            issues.append("No reads passed fastp filtering")

        mapped_pct = mapping_rates.get(sample, 0.0)
        min_mapped = thresholds.get("min_mapped_pct", 20.0)
        if mapped_pct < min_mapped:
            status = "FAIL"
            issues.append(f"Low mapping rate: {mapped_pct:.2f}% (threshold: {min_mapped}%)")

        results[sample] = {"status": status, "issues": issues}

    return results


def main():
    thresholds = {
        "min_read_length": 50,   # bp — discard reads shorter than this
        "min_q30_pct": 75,       # % bases with Phred ≥30
        "min_mapped_pct": 50.0,  # % reads mapped; <50% flags poor sample quality
    }

    try:
        fastp_files = snakemake.input.fastp
        flagstat_files = snakemake.input.flagstats
        output_file = snakemake.output[0]
    except NameError:
        fastp_files = sorted(Path("results/qc").glob("*_fastp.json"))
        flagstat_files = sorted(Path("results/mapping").glob("*.flagstat"))
        output_file = "results/qc/qc_summary.txt"

    metrics = parse_fastp_data(fastp_files)
    mapping_rates = parse_flagstat_data(flagstat_files)
    results = assess_qc(metrics, mapping_rates, thresholds)

    with open(output_file, "w") as f:
        f.write("=" * 70 + "\n")
        f.write("BMS503 QUALITY CONTROL ASSESSMENT REPORT\n")
        f.write("=" * 70 + "\n\n")

        for sample, result in sorted(results.items()):
            f.write(f"Sample: {sample}\n")
            f.write(f"  Status: {result['status']}\n")
            if result["issues"]:
                for issue in result["issues"]:
                    f.write(f"    - {issue}\n")
            else:
                f.write(f"    - No issues detected\n")
            f.write("\n")

        passed = sum(1 for r in results.values() if r["status"] == "PASS")
        total = len(results)
        f.write(f"\nSummary: {passed}/{total} samples passed QC\n")
        f.write("=" * 70 + "\n")

    print(f"QC assessment complete. {passed}/{total} samples passed. Report saved to {output_file}")


if __name__ == "__main__":
    main()
