###############################################
# BMS503 Genomic Surveillance Pipeline
###############################################

import os
import re
import subprocess
from pathlib import Path

# Pre-activate bms503-all for every rule shell regardless of invocation method
_conda_base = subprocess.run(
    ["conda", "info", "--base"], capture_output=True, text=True
).stdout.strip()
_env_lib = f"{_conda_base}/envs/bms503-all/lib"
# bcftools (bioconda) links against libgsl but expects cblas symbols that live in
# libgslcblas — preloading it makes them visible without patching the binary.
_gslcblas = f"{_env_lib}/libgslcblas.so.0"
_preload = f"export LD_PRELOAD={_gslcblas}; " if os.path.exists(_gslcblas) else ""
shell.prefix(f"export PATH={_conda_base}/envs/bms503-all/bin:$PATH; {_preload}")

# bcftools (bioconda) is compiled against libgsl.so.25 but conda installs GSL 2.7
# (libgsl.so.27). Create a symlink so the dynamic linker finds the library.
_gsl25 = f"{_env_lib}/libgsl.so.25"
_gsl27 = f"{_env_lib}/libgsl.so.27"
if not os.path.exists(_gsl25) and os.path.exists(_gsl27):
    os.symlink(_gsl27, _gsl25)

# ============================================================
# HELPER FUNCTIONS
# ============================================================

def get_reference_fasta():
    """Auto-detect reference FASTA file in data/reference/ directory"""
    ref_dir = config.get("reference_dir", "data/reference")
    if "reference_fasta" in config and os.path.exists(config["reference_fasta"]):
        return config["reference_fasta"]
    if os.path.isdir(ref_dir):
        for f in os.listdir(ref_dir):
            if f.endswith(('.fasta', '.fa', '.fna')):
                return os.path.join(ref_dir, f)
    raise FileNotFoundError(f"No reference genome found in {ref_dir}/")

REFERENCE = get_reference_fasta()

def get_raw_fastq_files(sample):
    """Return raw FASTQ input paths for a sample."""
    data_dir = config.get("data_dir", "data/raw")
    pair_patterns = [
        (f"{sample}_R1.fastq.gz", f"{sample}_R2.fastq.gz"),
        (f"{sample}_1.fastq.gz", f"{sample}_2.fastq.gz"),
        (f"{sample}.R1.fastq.gz", f"{sample}.R2.fastq.gz"),
    ]
    for r1, r2 in pair_patterns:
        r1_path = os.path.join(data_dir, r1)
        r2_path = os.path.join(data_dir, r2)
        if os.path.exists(r1_path) and os.path.exists(r2_path):
            return [r1_path, r2_path]
    for single in [f"{sample}.fastq.gz", f"{sample}.fq.gz", f"{sample}.fastq", f"{sample}.fq"]:
        path = os.path.join(data_dir, single)
        if os.path.exists(path):
            return [path]
    return []


def get_processed_fastq_files(sample):
    """Prefer fastp-trimmed reads for mapping."""
    raw = get_raw_fastq_files(sample)
    if len(raw) == 2:
        return [
            f"results/qc/{sample}_R1.trimmed.fastq.gz",
            f"results/qc/{sample}_R2.trimmed.fastq.gz",
        ]
    elif len(raw) == 1:
        return [f"results/qc/{sample}.trimmed.fastq.gz"]
    return []


def get_sample_list(data_dir):
    """Extract unique sample names from FASTQ files."""
    samples = set()
    for f in os.listdir(data_dir):
        if not f.endswith(('.fastq.gz', '.fastq', '.fq.gz', '.fq')):
            continue

        match = re.match(r"^(.+)(_R?[12])\.(fastq|fq)(\.gz)?$", f)
        if match:
            samples.add(match.group(1))
            continue

        match = re.match(r"^(.+)\.(R[12])\.(fastq|fq)(\.gz)?$", f)
        if match:
            samples.add(match.group(1))
            continue

        clean_name = f
        for suffix in [".fastq.gz", ".fq.gz", ".fastq", ".fq"]:
            if clean_name.endswith(suffix):
                clean_name = clean_name[: -len(suffix)]
        samples.add(clean_name)

    return sorted(samples)

SAMPLES = get_sample_list(config.get("data_dir", "data/raw"))

# ============================================================
# RULE: ALL
# ============================================================

rule all:
    input:
        expand("results/qc/{sample}_fastp.html", sample=SAMPLES),
        expand("results/qc/{sample}_fastp.json", sample=SAMPLES),
        "results/qc/qc_summary.txt",
        expand("results/mapping/{sample}.bam", sample=SAMPLES),
        expand("results/mapping/{sample}.bam.bai", sample=SAMPLES),
        expand("results/mapping/{sample}.flagstat", sample=SAMPLES),
        expand("results/mapping/{sample}.coverage.txt", sample=SAMPLES),
        expand("results/variants/{sample}.vcf.gz", sample=SAMPLES),
        "results/phylogeny/aligned_sequences.fasta",
        "results/phylogeny/nextclade.tsv",
        "results/phylogeny/phylogenetic_tree.nwk",
        "results/reports/diagnostic_report.html",
        "results/reports/diagnostic_report.txt",
        expand("data/reference_db/{sample}_consensus.fasta", sample=SAMPLES),
        "data/reference_db/master_alignment.fasta",
        "data/reference_db/master_tree.nwk"

# ============================================================
# STAGE 1: QUALITY CONTROL
# ============================================================

rule fastp:
    input:
        raw=lambda wc: get_raw_fastq_files(wc.sample)
    output:
        trimmed1="results/qc/{sample}_R1.trimmed.fastq.gz",
        trimmed2="results/qc/{sample}_R2.trimmed.fastq.gz",
        html="results/qc/{sample}_fastp.html",
        json="results/qc/{sample}_fastp.json"
    log:
        "logs/fastp/{sample}.log"
    threads: 2
    shell:
        """
        if [ -n "{input.raw[1]}" ]; then
            fastp -i {input.raw[0]} -I {input.raw[1]} -o {output.trimmed1} -O {output.trimmed2} \
                  -h {output.html} -j {output.json} -w {threads} > {log} 2>&1
        else
            fastp -i {input.raw[0]} -o {output.trimmed1} \
                  -h {output.html} -j {output.json} -w {threads} > {log} 2>&1
            touch {output.trimmed2}
        fi
        """

rule qc_threshold:
    input:
        fastp=expand("results/qc/{sample}_fastp.json", sample=SAMPLES),
        flagstats=expand("results/mapping/{sample}.flagstat", sample=SAMPLES)
    output:
        "results/qc/qc_summary.txt"
    log:
        "logs/qc_threshold.log"
    script:
        "workflow/scripts/qc_assessment.py"

# ============================================================
# STAGE 2: READ MAPPING AND ALIGNMENT
# ============================================================

rule index_reference:
    input:
        REFERENCE
    output:
        f"{REFERENCE}.fai"
    log:
        "logs/index_reference.log"
    shell:
        "samtools faidx {input} > {log} 2>&1"

rule bwa_index:
    input:
        REFERENCE
    output:
        expand(f"{REFERENCE}.{{ext}}", ext=["0123", "amb", "ann", "bwt.2bit.64", "pac"])
    log:
        "logs/bwa_index.log"
    shell:
        "bwa-mem2 index {input} > {log} 2>&1"

rule bwa_map:
    input:
        ref=REFERENCE,
        idx=expand(f"{REFERENCE}.{{ext}}", ext=["0123", "amb", "ann", "bwt.2bit.64", "pac"]),
        reads=lambda wc: get_processed_fastq_files(wc.sample)
    output:
        "results/mapping/{sample}.bam"
    log:
        "logs/bwa_map/{sample}.log"
    threads: 4
    shell:
        """
        bwa-mem2 mem -t {threads} {input.ref} {input.reads} 2> {log} |
        samtools sort -@ {threads} -o {output} - >> {log} 2>&1
        """

rule index_bam:
    input:
        "results/mapping/{sample}.bam"
    output:
        "results/mapping/{sample}.bam.bai"
    log:
        "logs/index_bam/{sample}.log"
    threads: 2
    shell:
        "samtools index -@ {threads} {input} > {log} 2>&1"

rule post_mapping:
    input:
        "results/mapping/{sample}.bam"
    output:
        flagstat="results/mapping/{sample}.flagstat",
        coverage="results/mapping/{sample}.coverage.txt"
    log:
        "logs/post_mapping/{sample}.log"
    threads: 2
    shell:
        """
        samtools flagstat -@ {threads} {input} > {output.flagstat} 2> {log}
        samtools coverage {input} > {output.coverage} 2>> {log}
        """

# ============================================================
# STAGE 3: VARIANT CALLING
# ============================================================

rule call_variants:
    input:
        bam="results/mapping/{sample}.bam",
        bai="results/mapping/{sample}.bam.bai",
        ref=REFERENCE
    output:
        vcf="results/variants/{sample}.vcf.gz",
        idx="results/variants/{sample}.vcf.gz.csi"
    log:
        "logs/call_variants/{sample}.log"
    threads: 2
    params:
        min_af=config.get("variants", {}).get("min_allele_frequency", 0.05)
    shell:
        """
        bcftools mpileup --threads {threads} -Ou -d 10000 -f {input.ref} {input.bam} |
        bcftools call -mv -Oz |
        bcftools view --threads {threads} --min-af {params.min_af}:nref -Oz -o {output.vcf} > {log} 2>&1
        bcftools index --threads {threads} {output.vcf} >> {log} 2>&1
        """

# ============================================================
# STAGE 4: PHYLOGENETIC ANALYSIS
# ============================================================

rule consensus_sequence:
    input:
        bam="results/mapping/{sample}.bam",
        ref=REFERENCE,
        vcf="results/variants/{sample}.vcf.gz",
        vcf_index="results/variants/{sample}.vcf.gz.csi"
    output:
        consensus="results/phylogeny/{sample}_consensus.fasta",
        db_copy="data/reference_db/{sample}_consensus.fasta"
    log:
        "logs/consensus/{sample}.log"
    shell:
        """
        # Skip bcftools consensus when VCF is empty — it segfaults on headerless input
        VARIANT_COUNT=$(bcftools view -H {input.vcf} 2>/dev/null | wc -l)
        if [ "$VARIANT_COUNT" -eq 0 ]; then
            echo "No variants called for {wildcards.sample} — using reference sequence" > {log}
            echo ">{wildcards.sample}" > {output.consensus}
            tail -n +2 {input.ref} >> {output.consensus}
        else
            set +e
            bcftools consensus -f {input.ref} {input.vcf} > {output.consensus}.tmp 2> {log}
            BCFTOOLS_EXIT=$?
            set -e
            if [ $BCFTOOLS_EXIT -eq 0 ]; then
                echo ">{wildcards.sample}" > {output.consensus}
                tail -n +2 {output.consensus}.tmp >> {output.consensus}
                rm {output.consensus}.tmp
            else
                echo "Warning: bcftools consensus failed for {wildcards.sample}, using reference sequence (Exit code: $BCFTOOLS_EXIT)" >> {log}
                echo ">{wildcards.sample}" > {output.consensus}
                tail -n +2 {input.ref} >> {output.consensus}
            fi
        fi
        cp {output.consensus} {output.db_copy}
        """

rule combine_sequences:
    """
    Combine the reference genome with consensus sequences from QC-passing samples
    only.  Failing samples remain in the QC table but are excluded from alignment
    and the tree so that low-quality sequences (e.g. reference-fallback consensus
    from bcftools crash) do not introduce spurious phylogenetic signal.
    """
    input:
        ref=REFERENCE,
        saved=expand("data/reference_db/{sample}_consensus.fasta", sample=SAMPLES),
        qc_summary="results/qc/qc_summary.txt"
    output:
        "results/phylogeny/combined_sequences.fasta"
    log:
        "logs/combine_sequences.log"
    shell:
        """
        python3 - <<'EOF'
import sys, shutil

passing = []
current = None
with open("{input.qc_summary}") as f:
    for line in f:
        line = line.strip()
        if line.startswith("Sample: "):
            current = line[8:]
        elif "Status: PASS" in line and current:
            passing.append(current)
            current = None

snakemake_samples = "{input.saved}".split()
snakemake_samples = [s.split("/")[-1].replace("_consensus.fasta", "") for s in snakemake_samples]

with open("{output}", "wb") as out:
    with open("{input.ref}", "rb") as r:
        shutil.copyfileobj(r, out)
    for sample in sorted(passing):
        path = f"data/reference_db/{{sample}}_consensus.fasta"
        try:
            with open(path, "rb") as f:
                shutil.copyfileobj(f, out)
        except FileNotFoundError:
            print(f"Warning: consensus not found for {{sample}}", file=sys.stderr)

total = sum(1 for line in open("{output}") if line.startswith(">"))
print(f"Combined {{total}} sequences ({{len(passing)}} passing samples + reference)")
print(f"Excluded: {{set(snakemake_samples) - set(passing)}}")
EOF
        """

rule nextclade_align:
    """
    Align consensus sequences against the SARS-CoV-2 reference and assign
    Pango lineages using Nextclade. Uses the local dataset downloaded by
    setup.sh; falls back to a live download if the directory is absent.
    """
    input:
        combined="results/phylogeny/combined_sequences.fasta"
    output:
        aligned="results/phylogeny/aligned_sequences.fasta",
        tsv="results/phylogeny/nextclade.tsv"
    log:
        "logs/nextclade_align.log"
    threads: workflow.cores
    params:
        dataset_dir=config.get("nextclade_dataset_dir", "")
    shell:
        """
        if [ -n "{params.dataset_dir}" ] && [ -d "{params.dataset_dir}" ]; then
            nextclade run \
                --input-dataset {params.dataset_dir} \
                --output-fasta {output.aligned} \
                --output-tsv {output.tsv} \
                --jobs {threads} \
                {input.combined} > {log} 2>&1
        else
            nextclade run \
                --dataset-name sars-cov-2 \
                --output-fasta {output.aligned} \
                --output-tsv {output.tsv} \
                --jobs {threads} \
                {input.combined} > {log} 2>&1
        fi
        """

rule build_tree:
    """
    Build phylogenetic tree.  On first run: standard VeryFastTree.
    On subsequent runs with no new taxa: use master tree as --intree starting
    topology (fast).  When new taxa are present, build from scratch — VeryFastTree
    --intree requires a strictly binary tree and cannot accept polytomies that
    would result from naively grafting new leaves onto the root.
    """
    input:
        aligned="results/phylogeny/aligned_sequences.fasta",
        ref=REFERENCE
    output:
        "results/phylogeny/phylogenetic_tree.nwk"
    log:
        "logs/build_tree.log"
    threads: workflow.cores
    shell:
        """
        set +e
        if [ -f "data/reference_db/master_tree.nwk" ]; then
            echo "Master tree found — checking for new/stale taxa" > {log}
            python3 workflow/scripts/graft_new_taxa.py \
                data/reference_db/master_tree.nwk \
                {input.aligned} \
                /tmp/bms503_starting_tree.nwk >> {log} 2>&1
            GRAFT_EXIT=$?
            set -e
            if [ $GRAFT_EXIT -eq 0 ]; then
                echo "No new taxa — using master tree as starting topology (fast mode)" >> {log}
                VeryFastTree -nt -gtr -threads {threads} -intree /tmp/bms503_starting_tree.nwk \
                    < {input.aligned} > {output} 2>> {log}
            else
                echo "New taxa detected — building tree from scratch" >> {log}
                VeryFastTree -nt -gtr -threads {threads} < {input.aligned} > {output} 2>> {log}
            fi
        else
            set -e
            echo "No master tree found — building tree from scratch" > {log}
            VeryFastTree -nt -gtr -threads {threads} < {input.aligned} > {output} 2>> {log}
        fi

        # Root with reference sequence as outgroup so the tree is interpretable
        REF_NAME=$(grep '^>' {input.ref} | head -1 | sed 's/^>//' | awk '{{print $1}}')
        echo "Rooting tree with outgroup: $REF_NAME" >> {log}
        python3 - <<PYEOF >> {log} 2>&1
from Bio import Phylo
ref_name = "$REF_NAME"
tree_file = "{output}"
try:
    tree = Phylo.read(tree_file, "newick")
    tree.root_with_outgroup({{"name": ref_name}})
    Phylo.write(tree, tree_file, "newick")
    print(f"Tree rooted with outgroup: {{ref_name}}")
except Exception as e:
    print(f"Warning: outgroup rooting failed ({{e}}); tree left unrooted")
PYEOF
        """

# ============================================================
# STAGE 5: DIAGNOSTIC REPORTING
# ============================================================

rule generate_report:
    input:
        qc="results/qc/qc_summary.txt",
        tree="results/phylogeny/phylogenetic_tree.nwk",
        aligned="results/phylogeny/aligned_sequences.fasta",
        variants=expand("results/variants/{sample}.vcf.gz", sample=SAMPLES),
        flagstats=expand("results/mapping/{sample}.flagstat", sample=SAMPLES),
        coverage=expand("results/mapping/{sample}.coverage.txt", sample=SAMPLES),
        tsv="results/phylogeny/nextclade.tsv"
    output:
        html="results/reports/diagnostic_report.html",
        txt="results/reports/diagnostic_report.txt",
        svg="results/phylogeny/phylogenetic_tree.nwk.svg"
    log:
        "logs/generate_report.log"
    script:
        "workflow/scripts/generate_report.py"

# ============================================================
# STAGE 6: UPDATE REFERENCE DATABASE MASTER FILES
# ============================================================

rule update_master_files:
    """
    Persist the final alignment and tree from this run as the master files
    used by subsequent runs for incremental alignment and tree building.
    Runs after generate_report to ensure it is the very last pipeline step.
    """
    input:
        alignment="results/phylogeny/aligned_sequences.fasta",
        tree="results/phylogeny/phylogenetic_tree.nwk",
        report="results/reports/diagnostic_report.html"
    output:
        master_alignment="data/reference_db/master_alignment.fasta",
        master_tree="data/reference_db/master_tree.nwk"
    log:
        "logs/update_master_files.log"
    shell:
        """
        cp {input.alignment} {output.master_alignment}
        cp {input.tree} {output.master_tree}
        SEQ_COUNT=$(grep -c '^>' {output.master_alignment})
        echo "Master alignment updated: $SEQ_COUNT sequences" > {log}
        echo "Master tree updated from {input.tree}" >> {log}

        # Record every sample name processed in this or any previous run.
        # generate_report reads this file (before it is written here) to detect
        # new samples — so it always reflects the state BEFORE the current run.
        python3 - <<'PYEOF'
from pathlib import Path
db = Path("data/reference_db")
# Every processed sample has a consensus file here regardless of QC status
samples = sorted(f.stem.replace("_consensus", "") for f in db.glob("*_consensus.fasta"))
tracking = db / "processed_samples.txt"
with open(tracking, "w") as f:
    for s in samples:
        f.write(s + "\\n")
print(f"Tracked {{len(samples)}} samples in {{tracking}}")
PYEOF
        """

# ============================================================
# END OF PIPELINE
# ============================================================

print("Pipeline configuration loaded successfully.")
print(f"Detected {len(SAMPLES)} sample(s). Reference genome: {REFERENCE}")
