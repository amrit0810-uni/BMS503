# BMS503 Genomic Surveillance Pipeline

**Version v1.8 — Production Ready**  
**See [CHANGELOG.md](CHANGELOG.md) for the complete history of all bug fixes and improvements.**

---

## Table of Contents

1. [Prerequisites & Installation](#prerequisites--installation)
2. [Quick Start](#quick-start)
3. [Pipeline Architecture](#pipeline-architecture)
4. [Prepare Your Data](#prepare-your-data)
5. [Run the Pipeline](#run-the-pipeline)
6. [View Results](#view-results)
7. [Configuration](#configuration)
8. [Organism Guide](#organism-guide)
9. [Performance Guide](#performance-guide)
10. [Running Multiple Datasets](#running-multiple-datasets)
11. [Sharing the Pipeline](#sharing-the-pipeline)
12. [Troubleshooting](#troubleshooting)
13. [Known Limitations](#known-limitations)
14. [GenAI Contribution Statement](#genai-contribution-statement)
15. [Statement of Contributions](#statement-of-contributions)

---

## Prerequisites & Installation

**Requirements**: Linux/Unix or WSL, [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Anaconda, 5–10 GB free disk space

If you don't have conda, install Miniconda first:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

### Setup (one time only)

Run `setup.sh` — it creates the conda environment, applies any runtime fixes, and sets up the required directories:

```bash
bash setup.sh
```

`setup.sh` uses mamba if available (faster) and falls back to conda. It is safe to re-run.

---

## Quick Start

```bash
# 1. First time on a new machine
bash setup.sh

# 2. Copy your data
cp /path/to/your/fastq/*.fastq.gz data/raw/
cp /path/to/reference.fasta       data/reference/

# 3. Run
./run_pipeline.sh
```

**View results**: Open `results/reports/diagnostic_report.html` in your browser.

See `QUICKSTART.md` for a concise reference card.

---

## Pipeline Architecture

5 automated stages, organism-agnostic:

```
Input FASTQ Files
       ↓
┌──────────────────────────────────────────┐
│ STAGE 1: QUALITY CONTROL                 │
│ - fastp (read QC, filtering, trimming)   │
│ - QC threshold assessment                │
│ - Post-mapping QC gate (mapping rate)    │
└──────────────────────────────────────────┘
       ↓
┌──────────────────────────────────────────┐
│ STAGE 2: READ MAPPING                    │
│ - Reference genome indexing (BWA-MEM2)   │
│ - Read alignment (BWA-MEM2)              │
│ - BAM sorting, indexing, coverage        │
└──────────────────────────────────────────┘
       ↓
┌──────────────────────────────────────────┐
│ STAGE 3: VARIANT CALLING                 │
│ - Pileup generation (bcftools mpileup)   │
│ - Variant identification (bcftools call) │
│ - Compressed VCF output                  │
└──────────────────────────────────────────┘
       ↓
┌──────────────────────────────────────────┐
│ STAGE 4: PHYLOGENETIC ANALYSIS           │
│ - Consensus sequence generation          │
│ - Reference-guided alignment (Nextclade) │
│ - Pango lineage assignment               │
│ - Incremental tree building (VeryFastTree) │
│ - Tree visualization (SVG)               │
└──────────────────────────────────────────┘
       ↓
┌──────────────────────────────────────────┐
│ STAGE 5: DIAGNOSTIC REPORTING            │
│ - Per-sample QC status + failure reasons │
│ - Coverage depth & breadth               │
│ - Variant counts, phylogenetic tree      │
│ - HTML + text report                     │
└──────────────────────────────────────────┘
```

### Tools

| Stage | Tool | Version | Purpose |
|-------|------|---------|---------|
| QC | fastp | 1.3.3 | Read QC, filtering, JSON metrics |
| Mapping | BWA-MEM2 | 2.2.1 | Align reads to reference |
| Mapping | SAMtools | 1.21 | BAM manipulation, coverage |
| Variants | bcftools | 1.21 | Variant calling |
| Phylogeny | Nextclade | 3.21.2 | Reference-guided alignment + lineage |
| Phylogeny | VeryFastTree | 4.0.5 | Phylogenetic tree inference (multi-threaded) |
| Workflow | Snakemake | 8.30.0 | Pipeline orchestration |
| Visualization | Biopython | 1.79 | Tree parsing and SVG output |

### Project Structure

```
BMS503/
├── Snakefile                      # Workflow definition
├── setup.sh                       # First-time setup (run once)
├── run_pipeline.sh                # Main entry script
├── config/
│   ├── config.yaml                # Pipeline settings
│   └── env_all.yml                # Conda environment (bms503-all) — all tools in one
├── workflow/scripts/
│   ├── qc_assessment.py           # QC metrics parsing + mapping-rate gate
│   ├── generate_report.py         # HTML/text report generation
│   ├── visualize_tree.py          # Tree SVG visualisation
│   └── graft_new_taxa.py          # Incremental tree helper
├── data/
│   ├── raw/                       # FASTQ files go here
│   ├── reference/                 # Reference genome goes here
│   └── reference_db/              # Accumulated consensus sequences (auto-populated)
├── results/                       # Pipeline outputs (auto-created)
│   ├── qc/
│   ├── mapping/
│   ├── variants/
│   ├── phylogeny/
│   └── reports/
├── logs/                          # Execution logs (auto-created)
└── CHANGELOG.md                   # Complete history of all changes
```

---

## Prepare Your Data

### FASTQ files → `data/raw/`

```bash
cp /path/to/your/fastq/*.fastq.gz data/raw/
ls -lah data/raw/   # Verify
```

Supported: `.fastq.gz`, `.fastq`, `.fq.gz`, `.fq`  
Paired-end (`_R1`/`_R2`) and single-end files are both auto-detected.

### Reference genome → `data/reference/`

```bash
cp /path/to/reference.fasta data/reference/
ls -lah data/reference/   # Verify
```

Supported: `.fasta`, `.fa`, `.fna`, `.ffn`

### Verify setup

Run `bash setup.sh` — it is safe to re-run and will report whether the environment already exists.

Expected:
```
✅ FASTQ files: X file(s) found
✅ Reference genome: 1 file(s) found
✅ Configuration: config/config.yaml exists
✅ Pipeline: Snakefile ready
```

---

## Run the Pipeline

```bash
# Run using all available cores (default, fastest)
./run_pipeline.sh

# Limit cores on a shared system
./run_pipeline.sh --cores 8

# Dry-run (preview, no execution)
./run_pipeline.sh --dry-run
```

The `bms503-all` environment is activated automatically by `run_pipeline.sh`. If you invoke `snakemake` directly, prepend the environment's bin to your PATH first:

```bash
export PATH="$(conda info --base)/envs/bms503-all/bin:$PATH"
snakemake -s Snakefile --cores $(nproc) --configfile config/config.yaml
```

### Monitor progress

```bash
watch -n 5 'ls -lah results/*'
```

---

## View Results

| Result | Location |
|--------|----------|
| **Main Report** | `results/reports/diagnostic_report.html` ⭐ |
| Per-sample QC | `results/qc/*_fastp.html` |
| QC summary | `results/qc/qc_summary.txt` |
| Phylogenetic Tree | `results/phylogeny/phylogenetic_tree.nwk.svg` |
| Nextclade lineages | `results/phylogeny/nextclade.tsv` |
| Variant Calls | `results/variants/*.vcf.gz` |
| Alignments | `results/mapping/*.bam` |
| Coverage | `results/mapping/*.coverage.txt` |

Open `results/reports/diagnostic_report.html` in your browser for the full analysis summary.

---

## Configuration

Edit `config/config.yaml`:

```yaml
# Organism name (shown in report)
organism: "SARS-CoV-2"

# QC thresholds (adjust for organism and platform)
qc_thresholds:
  min_read_length: 50       # Minimum read length (bp)
  min_mean_quality: 30      # Phred quality score
  max_n_content: 5          # Max % N bases
  min_mapped_reads: 50      # Min % reads mapped (post-mapping gate)
```

---

## Organism Guide

Only three things change between organisms: the organism name in config, the QC thresholds, and the reference genome.

### Switching organisms

1. Update `config/config.yaml`
2. Place the correct reference in `data/reference/`
3. Place your FASTQ files in `data/raw/`
4. Run: `./run_pipeline.sh`

### Virus configurations

#### SARS-CoV-2 (genome ~29.9 kb, Illumina/ONT)
```yaml
organism: "SARS-CoV-2"
qc_thresholds:
  min_read_length: 50
  min_mean_quality: 30
  max_n_content: 5
  min_mapped_reads: 50
```
Reference: NCBI NC_045512.2 (Wuhan-Hu-1)

#### Influenza A/B (genome ~13.6 kb, 8 segments)
```yaml
organism: "Influenza A/H1N1"
qc_thresholds:
  min_read_length: 40
  min_mean_quality: 28
  max_n_content: 5
  min_mapped_reads: 50
```
Reference: NCBI FluSurver / WHO FluNet  
Note: 8 RNA segments — consider analysing each separately or concatenated.

#### Monkeypox (genome ~197 kb)
```yaml
organism: "Monkeypox virus (MPXV)"
qc_thresholds:
  min_read_length: 75
  min_mean_quality: 30
  max_n_content: 3
  min_mapped_reads: 60
```
Reference: NCBI NC_063383.1

#### Ebola (genome ~19 kb)
```yaml
organism: "Ebola virus (EBOV)"
qc_thresholds:
  min_read_length: 50
  min_mean_quality: 30
  max_n_content: 5
  min_mapped_reads: 50
```

### Bacteria configurations

#### Mycobacterium tuberculosis (genome ~4.4 Mb, high GC ~65%)
```yaml
organism: "Mycobacterium tuberculosis"
qc_thresholds:
  min_read_length: 100
  min_mean_quality: 32
  max_n_content: 2
  min_mapped_reads: 60
```
Reference: NCBI NC_000962.3 (H37Rv)

#### Salmonella enterica (genome ~4.6–4.9 Mb)
```yaml
organism: "Salmonella enterica"
qc_thresholds:
  min_read_length: 80
  min_mean_quality: 30
  max_n_content: 3
  min_mapped_reads: 55
```

#### Escherichia coli (genome ~4.6 Mb)
```yaml
organism: "Escherichia coli"
qc_thresholds:
  min_read_length: 75
  min_mean_quality: 30
  max_n_content: 3
  min_mapped_reads: 55
```

### Fungi configurations

#### Candida auris (genome ~12 Mb)
```yaml
organism: "Candida auris"
qc_thresholds:
  min_read_length: 100
  min_mean_quality: 30
  max_n_content: 3
  min_mapped_reads: 60
```

### QC thresholds by platform

| Platform | min_read_length | min_mean_quality |
|----------|-----------------|------------------|
| Illumina (short-read) | 50–100 bp | 30–35 |
| Oxford Nanopore | 5000 bp | 10–15 |
| PacBio | 5000 bp | 20–25 |

### QC thresholds by genome size

| Genome size | min_read_length | min_mean_quality | max_n_content |
|------------|-----------------|------------------|---------------|
| < 20 kb (viruses) | 40–50 bp | 28–30 | 5% |
| 20–100 kb (large viruses) | 50–75 bp | 30–32 | 3–5% |
| 1–5 Mb (bacteria) | 75–100 bp | 30–32 | 2–3% |
| 5–20 Mb (fungi) | 100–150 bp | 30–32 | 2% |

### Performance by organism (16 vCPUs, observed)

| Organism | Genome size | Samples | Total time |
|----------|------------|---------|------------|
| SARS-CoV-2 | 29.9 kb | 16 | ~2 min 30s |
| SARS-CoV-2 | 29.9 kb | 27 | ~2 min 50s |
| Influenza | 13.6 kb | — | Expected faster (smaller genome) |
| M. tuberculosis | 4.4 Mb | — | Significantly longer (150× larger genome) |
| Candida auris | 12 Mb | — | Significantly longer (400× larger genome) |

Timings for non-SARS-CoV-2 organisms are estimates; actual runtime scales with genome size and read depth.

### Where to find reference genomes

- NCBI GenBank / RefSeq: https://www.ncbi.nlm.nih.gov/
- GISAID (respiratory viruses): https://www.gisaid.org/
- European Nucleotide Archive: https://www.ebi.ac.uk/ena/

---

## Performance Guide

### Threading architecture

| Stage | Threads/job | Notes |
|-------|-------------|-------|
| fastp (QC) | 2 | Runs one job per sample in parallel |
| BWA-MEM2 (mapping) | 4 | Runs one job per sample in parallel |
| bcftools (variants) | 2 | Runs one job per sample in parallel |
| Nextclade (alignment) | all cores | Single job, scales with core count |
| VeryFastTree (tree) | all cores | Multi-threaded; incremental mode gives additional speed-up |

With 16 cores, each parallel stage fills all cores cleanly (e.g. 4 concurrent BWA-MEM2 jobs × 4 threads = 16).

### Memory usage

| Stage | Approximate |
|-------|-------------|
| fastp | 1–2 GB |
| BWA-MEM2 mapping | 2–4 GB |
| bcftools | 1–2 GB |
| Nextclade | 2–4 GB |
| VeryFastTree | 0.5–1 GB |
| **Total** | **~8–16 GB** |

### Time vs sample count (SARS-CoV-2, 16 vCPUs)

| Samples | Time | Notes |
|---------|------|-------|
| 16 | ~2 min 30s | Observed |
| 27 | ~2 min 50s | Observed |
| ~50 | ~4–5 min | Estimated |
| ~100 | ~8–10 min | Estimated |

Time scales slowly with sample count because per-sample stages (QC, mapping, variant calling) run fully in parallel — the bottleneck is Nextclade alignment and VeryFastTree, both of which are single multi-threaded jobs that grow with alignment size, not sample count linearly.

### Optimisation tips

```bash
# Default — uses all available cores automatically
./run_pipeline.sh

# Leave headroom on shared systems
./run_pipeline.sh --cores 14

# Monitor during a run
watch -n 1 'top -b -n1 | head -20'
free -h
```

Placing `data/` on SSD or NVMe gives the biggest single-hardware improvement for the mapping stage.

To inspect or adjust per-stage thread counts:
```bash
grep -n "threads:" Snakefile
```

---

## Running Multiple Datasets

### Sequential

```bash
# Dataset 1
./run_pipeline.sh

# Clear and reload for Dataset 2
snakemake --delete-all-output
rm data/raw/* data/reference/*
cp dataset2/*.fastq.gz data/raw/
cp dataset2/reference.fasta data/reference/
sed -i 's/organism: .*/organism: "Influenza A\/H1N1"/' config/config.yaml
./run_pipeline.sh
```

### Parallel (separate directories)

```bash
cp -r . ../BMS503_covid
cp -r . ../BMS503_flu
# Analyse each organism in its own directory
```

---

## Sharing the Pipeline

### GitHub repository

The pipeline is hosted at:

**https://github.com/amrit0810-uni/BMS503**

Previous iteration archives (zip snapshots) are stored in `previous_iterations/`.

### Downloading the pipeline (new team members)

**Option A — Clone (recommended, keeps git history):**

```bash
git clone https://github.com/amrit0810-uni/BMS503.git
cd BMS503
bash setup.sh
```

**Option B — Download as zip (no git required):**

1. Go to https://github.com/amrit0810-uni/BMS503
2. Click **Code → Download ZIP**
3. Extract and run:

```bash
unzip BMS503-main.zip
cd BMS503-main
bash setup.sh   # use bash, not ./ — ensures permissions are set
```

**Option C — Previous iteration zip:**

Visit `previous_iterations/` on GitHub to download a specific archived version.

### What is tracked in git

```
✅ Snakefile
✅ config/config.yaml, config/env_all.yml
✅ workflow/scripts/*.py
✅ setup.sh, run_pipeline.sh, clean.sh
✅ README.md, QUICKSTART.md, CHANGELOG.md
✅ previous_iterations/   — archived zip snapshots
```

Data files are excluded automatically by `.gitignore`:

```
❌ data/            — all sequencing data, references, databases
❌ results/         — auto-generated outputs
❌ .snakemake/      — Snakemake cache
❌ *.pdf            — assessment documents
```

### Pushing updates

```bash
git add <changed files>
git commit -m "describe your changes"
git push
```

### Distribution checklist

- [ ] Run `bash setup.sh` on target machine (creates env + fixes permissions)
- [ ] `config/config.yaml` has correct organism and thresholds
- [ ] Tested on a fresh system
- [ ] `.gitignore` prevents large data commits

---

## Troubleshooting

### "No FASTQ files found"
- Check: `ls -la data/raw/`
- Extension must be `.fastq.gz`, `.fastq`, `.fq.gz`, or `.fq`
- Test: `gunzip -t data/raw/*.fastq.gz`

### "No reference genome found"
- Check: `ls -la data/reference/`
- Extension must be `.fasta`, `.fa`, `.fna`, or `.ffn`
- Test: `head -n 5 data/reference/*`

### "Snakemake not found" or "Tool not found" (bwa-mem2, samtools, fastp, etc.)
The `bms503-all` environment may not have been created yet:
```bash
mamba env create -f config/env_all.yml
# or force-overwrite an existing one:
mamba env create -f config/env_all.yml --force
```

### "Conda environment creation fails"
```bash
conda update -n base conda
mamba env create -f config/env_all.yml   # Re-create the environment
```

### Out of memory / CPU maxed
```bash
./run_pipeline.sh --cores 2   # Reduce parallelism
free -h                        # Check available RAM
```

### Pipeline hangs or seems stuck
```bash
htop        # Check CPU
df -h       # Check disk
# Stop with Ctrl+C and restart
```

### Different results on different machines
Force a fresh environment:
```bash
conda env remove -n bms503-all
mamba env create -f config/env_all.yml
./run_pipeline.sh
```

### High N content in consensus sequences
Low coverage or repetitive regions. Relax threshold if acceptable:
```yaml
qc_thresholds:
  max_n_content: 10
```

### Very few variants called
Choose a closer reference genome, or relax variant calling thresholds in config.yaml.

### Low mapping rate
Check reference is appropriate:
```bash
head -n 2 data/reference/reference_genome.fasta
samtools flagstat results/mapping/sample.bam
```

### Tree construction fails
```bash
# Check alignment
grep -c "^>" results/phylogeny/aligned_sequences.fasta
# Check Nextclade
nextclade --version
```

---

## Known Limitations

- **VeryFastTree**: Infers trees using the GTR+CAT model, which approximates rate variation across sites. Bootstrap values are localSH-like supports, not traditional bootstrap replicates. Results are appropriate for surveillance and epidemiological clustering but should be interpreted with caution for formal phylogenetic inference.
- **Samples with 0 variants**: May indicate poor mapping or failed QC — correctly identified by the post-mapping QC gate.
- **Nextclade lineage assignment**: Only works for organisms with an available Nextclade dataset (SARS-CoV-2, Influenza, RSV, etc.). For other organisms, the alignment stage still runs but lineage fields will be empty.

### Success indicators

- [ ] fastp QC reports: `ls results/qc/*_fastp.html`
- [ ] BAM files: `ls results/mapping/*.bam`
- [ ] Variants: `ls results/variants/*.vcf.gz`
- [ ] Tree: `ls results/phylogeny/phylogenetic_tree.nwk`
- [ ] Report: `ls results/reports/diagnostic_report.html`

---

## GenAI Contribution Statement

**GenAI Tools Used**: Microsoft Copilot, Claude Code (Anthropic)

### Microsoft Copilot

**How Used**:
- Initial ideation and conceptual design of the pipeline architecture
- Early brainstorming of the 5-stage workflow (QC → mapping → variant calling → phylogeny → reporting)
- Suggesting appropriate bioinformatics tools for each stage

### Claude Code (Anthropic)

**How Used**:
- Drafting and debugging Snakemake workflow rules and pipeline structure
- Creating and fixing Python helper scripts (QC assessment, report generation, tree visualisation, incremental tree grafting)
- Writing and updating shell scripts (`setup.sh`, `run_pipeline.sh`, `clean.sh`)
- Debugging all critical pipeline bugs across all phases (see [CHANGELOG.md](CHANGELOG.md))
- Writing and updating all documentation

**Validation & Testing**:
- All code was reviewed and tested by group members
- Each Snakemake rule was validated against bioinformatics best practices
- Pipeline tested end-to-end with 26 real SARS-CoV-2 samples: 23/26 passing QC
- Report accuracy and variant counts verified against raw output files

**What GenAI Did NOT Do**:
- Generate biological interpretations or conclusions
- Create fake results or data
- Replace group understanding of the analysis workflow

---

## Statement of Contributions

### Group Members and Roles

[**To be completed by group**: Each member should be listed with their primary contributions]

```
Member Name: [Name]
  Lead:       [Components where they had primary responsibility]
  Supporting: [Components they helped develop]
  Validation: [Testing and verification work]
```

---

**Last Updated**: May 19, 2026  
**Pipeline version**: v1.8
