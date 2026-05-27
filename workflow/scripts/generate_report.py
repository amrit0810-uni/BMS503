#!/usr/bin/env python3
"""
Diagnostic Report Generation Script
Creates final diagnostic report in HTML and text formats.

Report structure:
  1. Header
  2. AUTOMATED INTERPRETATION SUMMARY  ← top
  3. Run statistics
  4. Section 1 — Quality Control table
  5. Section 2 — Variant Analysis
  6. Section 3 — Phylogenetic Analysis
  7. Section 4 — Software & Tool Versions  ← bottom
  8. Footer

New-sample detection: samples absent from data/reference_db/master_alignment.fasta
(populated by the previous run's update_master_files rule) are automatically flagged
as [NEW] in the QC table and given focused attention in the interpretation.
No config field needed — just drop new FASTQs into data/raw/ and re-run.
"""

from datetime import datetime
import os
import re
import json
import gzip
import subprocess
import textwrap
from pathlib import Path


# ── DATA PARSING ──────────────────────────────────────────────────────────────

def count_variants_in_vcf(vcf_file):
    count = 0
    try:
        with gzip.open(vcf_file, 'rt') as f:
            for line in f:
                if not line.startswith('#'):
                    count += 1
    except Exception:
        pass
    return count


def parse_fastp_metrics(qc_dir="results/qc"):
    metrics = {}
    for fastp_json in Path(qc_dir).glob("*_fastp.json"):
        sample = fastp_json.stem.replace("_fastp", "")
        try:
            with open(fastp_json) as f:
                data = json.load(f)
            bf           = data.get("summary", {}).get("before_filtering", {})
            af           = data.get("summary", {}).get("after_filtering", {})
            fr           = data.get("filtering_result", {})
            before_reads = bf.get("total_reads", 0)
            after_reads  = af.get("total_reads", 0)
            passed_reads = fr.get("passed_filter_reads", 0)
            avg_length   = af.get("mean_length") or (
                af.get("total_bases", 0) / after_reads if after_reads else 0)
            q30_rate     = af.get("q30_rate", 0)
            gc_content   = af.get("gc_content", 0)
            if gc_content <= 1:
                gc_content *= 100
            metrics[sample] = {
                "before_reads": before_reads,
                "after_reads":  after_reads,
                "passed_reads": passed_reads,
                "avg_length":   avg_length,
                "q30_rate":     q30_rate * 100 if q30_rate else 0,
                "gc_content":   gc_content,
            }
        except Exception as e:
            print(f"Error parsing {fastp_json}: {e}")
    return metrics


def parse_qc_summary(qc_summary_file):
    status = {}
    current_sample = None
    try:
        with open(qc_summary_file) as f:
            for line in f:
                line = line.rstrip()
                if line.startswith("Sample: "):
                    current_sample = line[len("Sample: "):]
                    status[current_sample] = {"status": "PASS", "reasons": []}
                elif line.strip().startswith("Status: ") and current_sample:
                    status[current_sample]["status"] = line.strip()[len("Status: "):]
                elif line.strip().startswith("- ") and current_sample:
                    reason = line.strip()[2:]
                    if reason != "No issues detected":
                        status[current_sample]["reasons"].append(reason)
    except Exception as e:
        print(f"Warning: could not parse QC summary: {e}")
    return status


def parse_flagstat_files(mapping_dir="results/mapping"):
    rates = {}
    for f in Path(mapping_dir).glob("*.flagstat"):
        sample = f.stem
        try:
            for line in open(f):
                if "mapped (" in line:
                    m = re.search(r'\((\d+\.?\d*)%', line)
                    if m:
                        rates[sample] = float(m.group(1))
                    break
        except Exception as e:
            print(f"Warning: could not parse {f}: {e}")
    return rates


def parse_coverage_files(coverage_files):
    coverage = {}
    for cov_file in coverage_files:
        sample = Path(cov_file).stem.replace(".coverage", "")
        try:
            with open(cov_file) as f:
                headers = f.readline().strip().split('\t')
                idx     = {h: i for i, h in enumerate(headers)}
                mean_depth = 0.0
                breadth    = 0.0
                for line in f:
                    fields = line.strip().split('\t')
                    if not fields:
                        continue
                    try:
                        mean_depth += float(fields[idx['meandepth']])
                        breadth     = max(breadth, float(fields[idx['coverage']]))
                    except (KeyError, ValueError, IndexError):
                        pass
            coverage[sample] = {"mean_depth": mean_depth, "breadth": breadth}
        except Exception as e:
            print(f"Warning: could not parse {cov_file}: {e}")
            coverage[sample] = {"mean_depth": None, "breadth": None}
    return coverage


def parse_nextclade_tsv(tsv_file):
    lineages = {}
    try:
        with open(tsv_file) as f:
            headers = f.readline().strip().split('\t')
            idx     = {h: i for i, h in enumerate(headers)}
            for line in f:
                fields = line.strip().split('\t')
                if not fields or len(fields) < 2:
                    continue
                name = fields[idx['seqName']]
                lineages[name] = {
                    'pango': (fields[idx['Nextclade_pango']]
                              if 'Nextclade_pango' in idx and idx['Nextclade_pango'] < len(fields)
                              else 'N/A'),
                    'clade': (fields[idx['clade']]
                              if 'clade' in idx and idx['clade'] < len(fields)
                              else 'N/A'),
                    'nc_qc': (fields[idx['qc.overallStatus']]
                              if 'qc.overallStatus' in idx and idx['qc.overallStatus'] < len(fields)
                              else 'N/A'),
                }
    except Exception as e:
        print(f"Warning: could not parse Nextclade TSV: {e}")
    return lineages


# ── NEW-SAMPLE DETECTION ──────────────────────────────────────────────────────

def get_new_samples(current_samples,
                    tracking_file="data/reference_db/processed_samples.txt",
                    ref_db_dir="data/reference_db"):
    """
    Return (new_sample_set, is_first_run).

    Primary source: data/reference_db/processed_samples.txt — written by
    update_master_files at the END of each run, so it reflects PREVIOUS runs only.

    Fallback: if processed_samples.txt is absent (e.g. update_master_files never
    completed or the file was deleted), scan data/reference_db/*_consensus.fasta.
    Any sample that already has a consensus file there was processed in a prior run.

    is_first_run is True only when there is genuinely no prior data at all —
    no tracking file AND no consensus files in the reference database.
    """
    tracking = Path(tracking_file)
    if tracking.exists():
        try:
            prev = set(tracking.read_text().splitlines()) - {""}
            return set(current_samples) - prev, False
        except Exception as e:
            print(f"Warning: could not read tracking file: {e}")

    # Fallback: derive previously-processed set from existing consensus files
    db_dir = Path(ref_db_dir)
    if db_dir.is_dir():
        prev = {p.stem.replace("_consensus", "")
                for p in db_dir.glob("*_consensus.fasta")}
        if prev:
            return set(current_samples) - prev, False

    # Truly no prior data — genuine first run
    return set(current_samples), True


# ── TOOL VERSION DETECTION ────────────────────────────────────────────────────

def get_tool_versions():
    """Return list of (tool, version, purpose) by querying each installed tool."""
    def _run(cmd, pattern):
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
            out = r.stdout + r.stderr
            m   = re.search(pattern, out)
            return m.group(1) if m else (out.strip().split('\n')[0] or "detected")
        except Exception:
            return "not found"

    return [
        ("fastp",     _run(["fastp",     "--version"],  r"fastp\s+v?(\S+)"),
         "Read QC, adapter trimming, JSON metrics"),
        ("BWA-MEM2",  _run(["bwa-mem2", "version"],      r"(\d+\.\d+\S*)"),
         "Reference-guided read alignment (BWA-MEM2)"),
        ("SAMtools",  _run(["samtools",  "--version"],  r"samtools\s+(\S+)"),
         "BAM sorting, indexing, flagstat, coverage"),
        ("bcftools",  _run(["bcftools",  "--version"],  r"bcftools\s+(\S+)"),
         "Variant calling, VCF generation"),
        ("Nextclade", _run(["nextclade", "--version"],  r"(\d+\.\d+\.\d+)"),
         "Reference-guided alignment, Pango lineage assignment"),
        ("VeryFastTree", _run(["VeryFastTree", "-help"], r"FastTree\s+version\s+(\S+)"),
         "Maximum-likelihood phylogenetic tree (GTR+CAT)"),
    ]


# ── MODULAR INTERPRETATION BUILDERS ──────────────────────────────────────────
#
# Each function returns an HTML string for one thematic paragraph.
# build_interpretation() assembles them in order.

def _interp_qc(total, pass_count, failing_s, mean_map, mean_depth, qc_status):
    """QC pass/fail summary sentence."""
    parts = []
    parts.append(
        f"Analysis of {total} sample{'s' if total != 1 else ''} "
        f"yielded <strong>{pass_count}/{total}</strong> passing quality control."
    )
    if mean_map is not None:
        depth_clause = (f" and mean genome coverage depth of "
                        f"<strong>{mean_depth:.0f}&times;</strong>"
                        if mean_depth is not None else "")
        parts.append(
            f"Passing samples achieved a mean reference mapping rate of "
            f"<strong>{mean_map:.1f}%</strong>{depth_clause}."
        )
    if failing_s:
        fail_items = []
        for s in failing_s:
            r = (qc_status or {}).get(s, {}).get("reasons", [])
            fail_items.append(f"<em>{s}</em> ({r[0] if r else 'low quality'})")
        parts.append(
            f"Sample{'s' if len(failing_s) > 1 else ''} "
            f"{', '.join(fail_items)} failed QC and "
            f"{'were' if len(failing_s) > 1 else 'was'} excluded from "
            f"phylogenetic analysis."
        )
    else:
        parts.append("All samples passed QC and were included in phylogenetic analysis.")
    return " ".join(parts)


def _interp_variants(total_variants, total, avg_vars):
    """Variant diversity sentence."""
    if total_variants == 0:
        return ""
    return (
        f"A total of <strong>{total_variants}</strong> variant sites were identified "
        f"across {total} samples (mean {avg_vars:.1f} per sample), reflecting "
        f"nucleotide diversity relative to the reference genome."
    )


def _interp_new_samples(new_samples, is_first_run, qc_status, variant_counts,
                         mapping_rates, coverage_data, nextclade_data, organism):
    """
    Focused paragraph for samples added in this run.
    On first run: brief note that no baseline exists.
    On subsequent runs: per-sample QC + lineage summary for each new sample.
    """
    if is_first_run:
        return (
            f"This is the first pipeline run &mdash; no prior baseline exists. "
            f"All samples are newly processed. On subsequent runs, any sample "
            f"absent from the previous master alignment will be highlighted as "
            f"<strong>[NEW]</strong> in the QC table and summarised here."
        )
    if not new_samples:
        return "No new samples were detected compared to the previous run."

    sorted_new = sorted(new_samples)
    names_str  = ", ".join(f"<em>{s}</em>" for s in sorted_new)
    lines = [
        f"<strong>{len(new_samples)} new sample"
        f"{'s' if len(new_samples) > 1 else ''} added in this run</strong>: {names_str}."
    ]

    for sample in sorted_new:
        qc_st     = (qc_status or {}).get(sample, {}).get("status", "UNKNOWN")
        m_map     = (mapping_rates or {}).get(sample)
        m_depth   = ((coverage_data or {}).get(sample, {}).get("mean_depth"))
        m_vars    = (variant_counts or {}).get(sample, 0)
        nc        = (nextclade_data or {}).get(sample, {})
        m_lineage = nc.get("pango", "N/A")
        m_clade   = nc.get("clade",  "N/A")

        map_str   = f"{m_map:.1f}%"          if m_map   is not None else "N/A"
        depth_str = f"{m_depth:.0f}&times;"  if m_depth is not None else "N/A"

        if qc_st == "PASS":
            if m_lineage not in ("N/A", "unknown", ""):
                lineage_note = (
                    f"classified as {organism} Pango lineage "
                    f"<strong>{m_lineage}</strong> (clade {m_clade})"
                )
            else:
                lineage_note = (
                    "lineage unclassified &mdash; organism may differ from reference "
                    "(low-confidence Nextclade hit; consider Kraken2 for identification)"
                )
            lines.append(
                f"&bull; '<strong>{sample}</strong>': QC PASS, "
                f"mapping {map_str}, depth {depth_str}, "
                f"{lineage_note}, {m_vars} variant site{'s' if m_vars != 1 else ''}. "
                f"Phylogenetic placement shown in tree below."
            )
        else:
            fail_r = (qc_status or {}).get(sample, {}).get("reasons", [])
            if m_map is not None and m_map < 30:
                lines.append(
                    f"&bull; '<strong>{sample}</strong>': very low mapping rate ({map_str}) "
                    f"&mdash; likely <em>not</em> {organism}. Excluded from phylogenetic "
                    f"analysis. Taxonomic classification (e.g. Kraken2) is recommended to "
                    f"identify this organism."
                )
            else:
                lines.append(
                    f"&bull; '<strong>{sample}</strong>': QC FAIL "
                    f"({'; '.join(fail_r) if fail_r else 'see QC table'}). "
                    f"Excluded from phylogenetic analysis."
                )

    return " ".join(lines)


def _interp_phylogeny(pass_count):
    """Phylogenetic conclusion sentence."""
    if pass_count < 2:
        return ""
    return (
        f"Phylogenetic relationships among <strong>{pass_count}</strong> passing "
        f"samples were inferred using VeryFastTree (GTR+CAT model). The resulting tree "
        f"provides a genomic epidemiology snapshot; samples sharing a recent common "
        f"ancestor may indicate linked transmission events or shared geographic origin."
    )


def build_interpretation(qc_metrics, qc_status, variant_counts, total_variants,
                          pass_count, mapping_rates, coverage_data, nextclade_data,
                          organism, new_samples, is_first_run):
    """
    Assemble per-theme paragraphs into the full auto-interpretation block.
    Returns (html_text, plain_text).
    """
    total     = len(qc_metrics)
    failing_s = [s for s, d in (qc_status or {}).items() if d.get("status") == "FAIL"]
    passing_s = [s for s, d in (qc_status or {}).items() if d.get("status") == "PASS"]

    pass_rates  = [v for k, v in (mapping_rates or {}).items()
                   if k in passing_s and v is not None]
    mean_map    = sum(pass_rates) / len(pass_rates) if pass_rates else None

    pass_depths = [v.get("mean_depth") for k, v in (coverage_data or {}).items()
                   if k in passing_s and v.get("mean_depth") is not None]
    mean_depth  = sum(pass_depths) / len(pass_depths) if pass_depths else None

    avg_vars = total_variants / total if total else 0

    paragraphs = [
        _interp_qc(total, pass_count, failing_s, mean_map, mean_depth, qc_status),
        _interp_variants(total_variants, total, avg_vars),
        _interp_new_samples(new_samples, is_first_run, qc_status, variant_counts,
                            mapping_rates, coverage_data, nextclade_data, organism),
        _interp_phylogeny(pass_count),
    ]

    html_para  = " ".join(p for p in paragraphs if p)
    plain_para = re.sub(r'<[^>]+>', '', html_para).replace("&times;", "x") \
                                                   .replace("&mdash;", " - ") \
                                                   .replace("&bull;", "*") \
                                                   .replace("&amp;",  "&") \
                                                   .replace("&minus;", "-")
    return html_para, plain_para


# ── HTML REPORT ───────────────────────────────────────────────────────────────

def generate_html_report(output_file, qc_summary, qc_metrics, variants_dir,
                          tree_svg, organism="Pathogen",
                          qc_status=None, mapping_rates=None,
                          nextclade_data=None, coverage_data=None,
                          new_samples=None, is_first_run=True):

    new_samples = new_samples or set()

    # Variant counts
    total_variants = 0
    variant_counts = {}
    if Path(variants_dir).exists():
        for vcf in Path(variants_dir).glob("*.vcf.gz"):
            sample = vcf.stem.replace(".vcf", "")
            c = count_variants_in_vcf(str(vcf))
            variant_counts[sample] = c
            total_variants += c

    pass_count   = sum(1 for d in (qc_status or {}).values() if d.get("status") == "PASS")
    avg_variants = total_variants / len(qc_metrics) if qc_metrics else 0

    # Tree (inline SVG → self-contained HTML)
    tree_embed = (Path(tree_svg).read_text()
                  if tree_svg and Path(tree_svg).exists()
                  else "<p><em>Phylogenetic tree image not available.</em></p>")

    # QC table rows
    qc_table_rows = ""
    for sample in sorted(qc_metrics.keys()):
        metrics      = qc_metrics[sample]
        sample_qc    = (qc_status or {}).get(sample, {})
        qc_st        = sample_qc.get("status", "PASS" if metrics["passed_reads"] > 0 else "FAIL")
        status_class = "status-pass" if qc_st == "PASS" else "status-fail"
        map_pct      = (mapping_rates or {}).get(sample)
        map_str      = f"{map_pct:.1f}%" if map_pct is not None else "N/A"
        nc           = (nextclade_data or {}).get(sample, {})
        lineage      = nc.get("pango", "N/A")
        clade        = nc.get("clade",  "N/A")
        cov          = (coverage_data or {}).get(sample, {})
        mean_depth   = cov.get("mean_depth")
        breadth      = cov.get("breadth")
        depth_str    = f"{mean_depth:.1f}&times;" if mean_depth is not None else "N/A"
        breadth_str  = f"{breadth:.1f}%"          if breadth    is not None else "N/A"
        is_new       = (sample in new_samples) and (not is_first_run)
        row_style    = ' style="background-color:#fff9c4;font-weight:bold;"' if is_new else ""

        qc_table_rows += f"""        <tr{row_style}>
            <td>{sample}{"&nbsp;[NEW]" if is_new else ""}</td>
            <td class="{status_class}">{qc_st}</td>
            <td>{metrics['before_reads']:,}</td>
            <td>{metrics['after_reads']:,}</td>
            <td>{metrics['passed_reads']:,}</td>
            <td>{metrics['avg_length']:.1f}</td>
            <td>{metrics['q30_rate']:.1f}%</td>
            <td>{metrics['gc_content']:.1f}%</td>
            <td>{map_str}</td>
            <td>{depth_str}</td>
            <td>{breadth_str}</td>
            <td>{lineage}</td>
            <td>{clade}</td>
            <td>{variant_counts.get(sample, 0)}</td>
        </tr>\n"""

    # Auto-interpretation (top of report)
    interp_html, _ = build_interpretation(
        qc_metrics, qc_status, variant_counts, total_variants,
        pass_count, mapping_rates, coverage_data, nextclade_data,
        organism, new_samples, is_first_run,
    )

    # Tool versions table
    tool_rows = ""
    for tool, version, purpose in get_tool_versions():
        tool_rows += (f"        <tr>"
                      f"<td><strong>{tool}</strong></td>"
                      f"<td><code>{version}</code></td>"
                      f"<td>{purpose}</td>"
                      f"</tr>\n")

    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>BMS503 Genomic Surveillance &mdash; Diagnostic Report</title>
    <style>
        body           {{ font-family: Arial, sans-serif; margin: 20px; background-color: #f5f5f5; }}
        .header        {{ background-color: #003366; color: white; padding: 20px; border-radius: 5px; }}
        .header p      {{ margin: 4px 0; }}
        .section       {{ background-color: white; padding: 20px; margin-top: 20px;
                          border-left: 4px solid #003366; border-radius: 3px; overflow-x: auto; }}
        .section h2    {{ color: #003366; margin-top: 0; }}
        .interp-box    {{ background-color: #e8f5e9; border-left: 6px solid #2e7d32;
                          padding: 18px 22px; margin-top: 20px; border-radius: 3px; }}
        .interp-box h2 {{ color: #1b5e20; margin-top: 0; font-size: 1.05em; }}
        .interp-box p  {{ margin: 0; line-height: 1.75; }}
        table          {{ width: 100%; border-collapse: collapse; margin-top: 10px; font-size: 12px; }}
        table th       {{ background-color: #e0e0e0; padding: 8px; text-align: left;
                          font-weight: bold; border-bottom: 2px solid #999; }}
        table td       {{ padding: 8px; border-bottom: 1px solid #ddd; vertical-align: top; }}
        table tr:nth-child(even) {{ background-color: #f9f9f9; }}
        .status-pass   {{ color: #2e7d32; font-weight: bold; }}
        .status-fail   {{ color: #c62828; font-weight: bold; }}
        .summary       {{ background-color: #f0f0f0; padding: 12px; border-radius: 3px; margin: 10px 0; }}
        .summary p     {{ margin: 4px 0; }}
        pre            {{ background: #1e1e1e; color: #d4d4d4; padding: 16px; border-radius: 4px;
                          font-size: 12px; overflow-x: auto; white-space: pre; line-height: 1.5; }}
        img, svg       {{ max-width: 100%; height: auto; margin: 20px 0; }}
        .footer        {{ font-size: 12px; color: #666; margin-top: 40px;
                          padding-top: 20px; border-top: 1px solid #ddd; }}
        code           {{ background: #f4f4f4; padding: 1px 4px; border-radius: 3px; font-size: 0.9em; }}
    </style>
</head>
<body>

<div class="header">
    <h1>BMS503 Genomic Surveillance Pipeline</h1>
    <p>Diagnostic Report &mdash; {organism} Sequence Analysis</p>
    <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
</div>

<!-- ══ AUTOMATED INTERPRETATION (top of report) ══ -->
<div class="interp-box">
    <h2>&#x1F4CB; Automated Interpretation Summary</h2>
    <p>{interp_html}</p>
</div>

<!-- ══ RUN STATISTICS ══ -->
<div class="section">
    <h2>Run Statistics</h2>
    <div class="summary">
        <p><strong>Total Samples Analyzed:</strong> {len(qc_metrics)}</p>
        <p><strong>Samples Passing QC:</strong> {pass_count}</p>
        <p><strong>New Samples This Run:</strong> {len(new_samples) if not is_first_run else "First run &mdash; all {len(qc_metrics)} samples newly processed"}</p>
        <p><strong>Total Variants Identified:</strong> {total_variants}</p>
        <p><strong>Average Variants per Sample:</strong> {avg_variants:.1f}</p>
    </div>
</div>

<!-- ══ SECTION 1: QC TABLE ══ -->
<div class="section">
    <h2>1. Quality Control Summary</h2>
    <table>
        <tr>
            <th>Sample</th>
            <th>Status</th>
            <th>Reads Before</th>
            <th>Reads After</th>
            <th>Passed QC</th>
            <th>Avg Length</th>
            <th>Q30 Rate</th>
            <th>GC Content</th>
            <th>Mapped %</th>
            <th>Mean Depth</th>
            <th>Breadth</th>
            <th>Lineage (Pango)</th>
            <th>Clade</th>
            <th>Variants</th>
        </tr>
{qc_table_rows}    </table>
    <p style="font-size:11px;color:#888;margin-top:8px;">
        [NEW] = sample absent from previous run&apos;s master alignment (added this run). On first run, no [NEW] tags are shown as all samples are newly processed.
        Highlighted rows are automatically detected &mdash; no config required.
        Mean Depth = mean read depth across the reference.
        Breadth = fraction of reference positions covered at &ge;1&times;.
    </p>
</div>

<!-- ══ SECTION 2: VARIANTS ══ -->
<div class="section">
    <h2>2. Variant Analysis</h2>
    <p>Total variant sites identified: <strong>{total_variants}</strong>
       (mean {avg_variants:.1f} per sample).
       Variants filtered at minimum allele frequency 5% and minimum quality 30.</p>
    <p>Per-sample VCF files: <code>results/variants/</code></p>
</div>

<!-- ══ SECTION 3: PHYLOGENETICS ══ -->
<div class="section">
    <h2>3. Phylogenetic Analysis</h2>
    <p>Tree inferred from Nextclade reference-guided alignment using VeryFastTree (GTR+CAT).
       Branch lengths = substitutions per site. Only QC-passing samples included.</p>
    {tree_embed}
    <p><strong>Methods &amp; Limitations:</strong> VeryFastTree (GTR+CAT) is suitable for rapid
    genomic surveillance but does not model recombination or gene conversion. Bootstrap
    values are local SH-like supports. For formal inference in a clinical setting,
    IQ-TREE with a time-calibrated model or the Nextstrain augur pipeline is preferred.</p>
</div>

<!-- ══ SECTION 4: TOOL VERSIONS (bottom) ══ -->
<div class="section">
    <h2>4. Software &amp; Tool Versions</h2>
    <p>All tools installed in a single conda environment (<code>bms503-all</code>),
       defined in <code>config/env_all.yml</code>. Versions recorded at report-generation time.</p>
    <table>
        <tr><th>Tool</th><th>Version</th><th>Role in Pipeline</th></tr>
{tool_rows}    </table>
</div>

<div class="footer">
    <p>Generated by the BMS503 Genomic Surveillance Pipeline.
       For questions contact the BMS503 team.</p>
</div>

</body>
</html>"""

    with open(output_file, 'w') as f:
        f.write(html_content)
    print(f"HTML report generated: {output_file}")


# ── TEXT REPORT ───────────────────────────────────────────────────────────────

def generate_text_report(output_file, qc_summary, organism="Pathogen",
                          qc_metrics=None, qc_status=None, variant_counts=None,
                          total_variants=0, pass_count=0, mapping_rates=None,
                          coverage_data=None, nextclade_data=None,
                          new_samples=None, is_first_run=True):

    new_samples = new_samples or set()
    _, interp_plain = build_interpretation(
        qc_metrics or {}, qc_status, variant_counts or {}, total_variants,
        pass_count, mapping_rates, coverage_data, nextclade_data,
        organism, new_samples, is_first_run,
    )

    with open(output_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("BMS503 GENOMIC SURVEILLANCE PIPELINE — DIAGNOSTIC REPORT\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Organism:  {organism}\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        f.write("AUTOMATED INTERPRETATION SUMMARY\n")
        f.write("-" * 80 + "\n")
        for para in textwrap.wrap(interp_plain, width=78):
            f.write(para + "\n")
        f.write("\n")

        f.write("QUALITY CONTROL RESULTS\n")
        f.write("-" * 80 + "\n")
        if os.path.exists(qc_summary):
            with open(qc_summary) as qc_file:
                f.write(qc_file.read())
        else:
            f.write("QC summary not available\n")
        f.write("\n")

        f.write("PIPELINE STAGES\n")
        f.write("-" * 80 + "\n")
        f.write("  1. Quality Control      fastp\n")
        f.write("  2. Read Mapping         BWA-MEM2 -> SAMtools sort / index / coverage\n")
        f.write("  3. Variant Calling      bcftools mpileup | call | view\n")
        f.write("  4. Alignment            Nextclade (reference-guided, Pango lineage)\n")
        f.write("  5. Phylogenetic Tree    VeryFastTree (GTR+CAT model)\n\n")

        f.write("SOFTWARE VERSIONS\n")
        f.write("-" * 80 + "\n")
        for tool, version, purpose in get_tool_versions():
            f.write(f"  {tool:<14} {version:<22} {purpose}\n")
        f.write("\n")

        f.write("RESULTS LOCATION\n")
        f.write("-" * 80 + "\n")
        f.write("  results/qc/         Quality control outputs\n")
        f.write("  results/mapping/    BAM files, flagstat, coverage\n")
        f.write("  results/variants/   VCF files\n")
        f.write("  results/phylogeny/  Alignment, tree, Nextclade TSV\n")
        f.write("  results/reports/    This report\n\n")

        f.write("=" * 80 + "\n")

    print(f"Text report generated: {output_file}")


# ── ENTRY POINT ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    try:
        qc_summary   = snakemake.input.qc
        output_html  = snakemake.output.html
        output_txt   = snakemake.output.txt
        tree_svg     = snakemake.output.svg
        variants_dir = "results/variants"

        # Generate SVG before building the report that embeds it
        import sys as _sys
        import os as _os
        _sys.path.insert(0, _os.path.dirname(__file__))
        from visualize_tree import tree_to_svg
        tree_to_svg(snakemake.input.tree, tree_svg)

        try:
            import yaml
            with open("config/config.yaml") as f:
                config_data = yaml.safe_load(f)
            organism = config_data.get("organism", "Pathogen")
        except Exception:
            organism = "SARS-CoV-2"

        qc_metrics     = parse_fastp_metrics("results/qc")
        qc_status      = parse_qc_summary(qc_summary)
        mapping_rates  = parse_flagstat_files("results/mapping")
        coverage_data  = parse_coverage_files(snakemake.input.coverage)
        nextclade_data = parse_nextclade_tsv(snakemake.input.tsv)

        variant_counts = {}
        if Path(variants_dir).exists():
            for vcf in Path(variants_dir).glob("*.vcf.gz"):
                sample = vcf.stem.replace(".vcf", "")
                variant_counts[sample] = count_variants_in_vcf(str(vcf))
        total_variants = sum(variant_counts.values())
        pass_count     = sum(1 for d in (qc_status or {}).values() if d.get("status") == "PASS")

        new_samples, is_first_run = get_new_samples(list(qc_metrics.keys()))

        generate_html_report(
            output_html, qc_summary, qc_metrics, variants_dir,
            tree_svg, organism,
            qc_status=qc_status,
            mapping_rates=mapping_rates,
            nextclade_data=nextclade_data,
            coverage_data=coverage_data,
            new_samples=new_samples,
            is_first_run=is_first_run,
        )
        generate_text_report(
            output_txt, qc_summary, organism,
            qc_metrics=qc_metrics,
            qc_status=qc_status,
            variant_counts=variant_counts,
            total_variants=total_variants,
            pass_count=pass_count,
            mapping_rates=mapping_rates,
            coverage_data=coverage_data,
            nextclade_data=nextclade_data,
            new_samples=new_samples,
            is_first_run=is_first_run,
        )

    except (NameError, TypeError):
        qc_metrics = parse_fastp_metrics("results/qc")
        organism   = "SARS-CoV-2"
        new_samples, is_first_run = get_new_samples(list(qc_metrics.keys()))
        generate_html_report(
            "results/reports/diagnostic_report.html",
            "results/qc/qc_summary.txt",
            qc_metrics, "results/variants",
            "results/phylogeny/phylogenetic_tree.nwk.svg",
            organism,
            new_samples=new_samples,
            is_first_run=is_first_run,
        )
        generate_text_report(
            "results/reports/diagnostic_report.txt",
            "results/qc/qc_summary.txt",
            organism,
            qc_metrics=qc_metrics,
            new_samples=new_samples,
            is_first_run=is_first_run,
        )
