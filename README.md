# BMS503 Genomic Surveillance Pipeline

**Version v1.8 — Production Ready**  
**See [CHANGELOG.md](CHANGELOG.md) for the complete history of all bug fixes and improvements.**

---

## Table of Contents

1. [Prerequisites & Installation](#prerequisites--installation)
2. [Quick Start](#quick-start)
3. [Pipeline Architecture](#pipeline-architecture)
4. [Tool Selection Rationale](#tool-selection-rationale)
5. [Prepare Your Data](#prepare-your-data)
6. [Run the Pipeline](#run-the-pipeline)
7. [View Results](#view-results)
8. [Configuration](#configuration)
9. [Performance Guide](#performance-guide)
10. [Sharing the Pipeline](#sharing-the-pipeline)
11. [Troubleshooting](#troubleshooting)
12. [Known Limitations](#known-limitations)
13. [GenAI Contribution Statement](#genai-contribution-statement)
14. [Statement of Contributions](#statement-of-contributions)
15. [References](#references)

---

## Prerequisites & Installation

**Requirements**: Linux/Unix or WSL, [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Anaconda, mamba, 5–10 GB free disk space

### Get the code

Clone the repository and enter the project directory:

```bash
git clone https://github.com/amrit0810-uni/BMS503.git
cd BMS503
```

If you don't have git, you can also download a zip from the repository page — see [Sharing the Pipeline](#sharing-the-pipeline) for details.

### Install conda (if not already installed)

If you don't have conda, install Miniconda first:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

Then install mamba (required — conda alone is too slow for this environment):

```bash
conda install -n base -c conda-forge mamba
```

### Setup (one time only)

Run `setup.sh` — it creates the conda environment, applies any runtime fixes, and sets up the required directories:

```bash
bash setup.sh
```

`setup.sh` requires mamba and will exit with instructions if it is not found. It is safe to re-run.

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

5 automated stages for SARS-CoV-2 genomic surveillance:

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

## Tool Selection Rationale

Each tool was selected against available alternatives for the specific demands of viral genomic surveillance: speed, correctness on near-haploid genomes, and suitability for Illumina short reads. File types are shown to illustrate how data flows through the pipeline.

### Quality Control — fastp

| | Files |
|---|---|
| **Input** | `.fastq.gz` / `.fastq` — raw paired-end or single-end reads |
| **Output** | `.fastq.gz` (trimmed reads), `.json` (QC metrics), `.html` (QC report) |

**Chosen over**: Trimmomatic, FastQC + Trim Galore

fastp performs adapter trimming, quality filtering, and QC reporting in a single pass, producing structured JSON output that the pipeline parses directly for per-sample QC metrics ([Chen, 2025](#ref-chen-2025)). Trimmomatic requires a separate tool (FastQC) for reporting and outputs plain text; Trim Galore is a wrapper around Cutadapt and FastQC that adds dependencies without adding capability relevant to this workflow. fastp is substantially faster than both alternatives due to multi-threading and a single-pass design.

### Read Mapping — BWA-MEM2

| | Files |
|---|---|
| **Input** | `.fastq.gz` (trimmed reads), `.fasta` (reference genome + index files: `.0123`, `.bwt.2bit.64`, `.amb`, `.ann`, `.pac`) |
| **Output** | `.bam` (sorted aligned reads) |

**Chosen over**: BWA (0.7.17), Bowtie2, Minimap2

BWA-MEM2 is the direct successor to BWA-MEM, producing identical alignments while running ~3× faster through SIMD vectorisation ([Vasimuddin et al., 2019](#ref-vasimuddin-2019)). It is the established standard for Illumina short-read mapping in variant-calling workflows. Bowtie2 is optimised for RNA-seq and has lower sensitivity for structural variants. Minimap2 is designed for long reads (ONT/PacBio) and is not recommended for Illumina short reads.

### BAM Processing — SAMtools

| | Files |
|---|---|
| **Input** | `.bam` (sorted alignments) |
| **Output** | `.bam.bai` (index), `.flagstat` (mapping statistics), `.coverage.txt` (depth and breadth per position) |

**Chosen over**: Picard, sambamba

SAMtools is the de facto standard for BAM manipulation in bioinformatics pipelines, maintained by the same team as bcftools and sharing a common C library (htslib) ([Danecek et al., 2021](#ref-danecek-2021)). This tight integration eliminates format-conversion overhead between the mapping and variant-calling stages. Picard is Java-based (slower startup, higher memory overhead) and primarily targets duplicate marking, which is not required for viral surveillance. sambamba offers faster sorting but lacks the breadth of utilities SAMtools provides.

### Variant Calling — bcftools

| | Files |
|---|---|
| **Input** | `.bam` + `.bam.bai` (indexed alignments), `.fasta` (reference) |
| **Output** | `.vcf.gz` (compressed variant calls), `.vcf.gz.csi` (index) |

**Chosen over**: GATK HaplotypeCaller, FreeBayes, DeepVariant

bcftools mpileup/call is optimised for haploid and near-haploid genomes — the correct model for viral surveillance, where a single consensus sequence per sample is the target ([Danecek et al., 2021](#ref-danecek-2021)). GATK HaplotypeCaller is designed for diploid human genomes and introduces unnecessary complexity (gVCF mode, joint genotyping) for this use case. FreeBayes supports haploid mode but is significantly slower. DeepVariant requires GPU resources and is trained on human germline data, making it poorly suited for viral sequencing.

### Consensus Generation — bcftools consensus

| | Files |
|---|---|
| **Input** | `.vcf.gz` (variant calls), `.fasta` (reference) |
| **Output** | `.fasta` (per-sample consensus sequence) |

Variant calls are applied to the reference sequence to produce a per-sample consensus FASTA. Samples with zero variants (failed mapping or poor coverage) fall back to a reference copy rather than crashing, which allows the pipeline to continue and flag those samples in the QC report.

### Alignment & Lineage Assignment — Nextclade

| | Files |
|---|---|
| **Input** | `.fasta` (all consensus sequences, multi-sample) |
| **Output** | `.fasta` (reference-guided multiple sequence alignment), `.tsv` (lineage assignments, QC flags, mutation calls) |

**Chosen over**: MAFFT + Pangolin (separate tools), Augur

Nextclade performs reference-guided multiple sequence alignment and Pango lineage assignment in a single step, using curated datasets maintained by the Nextstrain team ([Aksamentov et al., 2021](#ref-aksamentov-2021)). Running MAFFT and Pangolin separately would require two tool invocations, an intermediate file, and separate dataset management. Augur (the Nextstrain CLI) offers similar functionality but is designed for the full Nextstrain pipeline and carries more dependencies than are needed here.

### Phylogenetic Tree Inference — VeryFastTree

| | Files |
|---|---|
| **Input** | `.fasta` (multiple sequence alignment from Nextclade) |
| **Output** | `.nwk` (Newick format phylogenetic tree), `.svg` (tree visualisation) |

**Chosen over**: IQ-TREE2, RAxML-NG, FastTree2

VeryFastTree is a multi-threaded reimplementation of FastTree2 that uses the GTR+CAT model — appropriate for surveillance-scale datasets where speed and scalability matter more than exhaustive model selection ([Piñeiro et al., 2020](#ref-pineiro-2020); [Piñeiro & Pichel, 2024](#ref-pineiro-2024)). FastTree2 itself is single-threaded and does not scale to all available cores. IQ-TREE2 provides better statistical support (proper ultrafast bootstrap) and automatic model selection via ModelFinder, but is 3–10× slower for equivalent sample counts; for 17 SARS-CoV-2 samples this adds roughly 1–5 minutes to the tree step. RAxML-NG is similarly rigorous and similarly slower. For formal publication-quality phylogenetics, IQ-TREE2 would be the preferred choice.

### Workflow Orchestration — Snakemake

| | Files |
|---|---|
| **Input** | `Snakefile` (workflow rules), `config/config.yaml` (parameters), `config/env_all.yml` (environment) |
| **Output** | Coordinates all stages; produces no files directly |

**Chosen over**: Nextflow, Make

Snakemake is Python-based, making it directly accessible to bioinformaticians already using Python, and integrates natively with conda environments — allowing the entire tool stack to be declared in a single `env_all.yml` and activated automatically ([Mölder et al., 2021](#ref-molder-2021)). Its rule-based, file-driven execution model means only rules whose inputs have changed are re-run, which is critical for iterative surveillance workflows where new samples are added to an existing dataset. Nextflow uses a Groovy DSL that is less familiar and requires more boilerplate for conda integration. Make lacks conda awareness and is not designed for bioinformatics scheduling.

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

Quick checks before running:

```bash
ls data/raw/*.fastq.gz        # Should list your FASTQ files
ls data/reference/*.fasta     # Should show NC_045512.fasta (or similar)
```

`run_pipeline.sh` will report what it finds and exit with a clear message if anything is missing.

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

The pipeline is pre-configured for SARS-CoV-2 Illumina short-read data. The defaults in `config/config.yaml` are set for SARS-CoV-2 and do not need to be changed for standard runs.

### QC and variant calling thresholds

| Parameter | Default | Basis |
|-----------|---------|-------|
| `min_read_length` | 50 bp | See below |
| `min_q30_pct` | 75% | See below |
| `min_mapped_reads` | 50% | See below |
| `min_allele_frequency` | 0.1 (10%) | See below |

### Threshold rationale and sources

**`min_read_length: 50 bp`**
Illumina paired-end reads for SARS-CoV-2 are typically 150–250 bp. A 50 bp minimum discards adapter dimers and severely degraded fragments while retaining all reads of biological value. This is consistent with the ARTIC Network SARS-CoV-2 sequencing protocol ([Tyson et al., 2020](#ref-tyson-2020)) and the PHA4GE bioinformatics pipeline recommendations for SARS-CoV-2 ([Griffiths et al., 2022](#ref-griffiths-2022)), and is widely adopted in published surveillance pipelines including the COVID-19 Genomics UK (COG-UK) consortium workflow ([COVID-19 Genomics UK Consortium, 2020](#ref-coguk-2020)).

**`min_q30_pct: 75%`**
Illumina's platform specification defines a passing run as one in which ≥75% of bases achieve Phred quality Q30 (1 error per 1,000 bases). This threshold is cited in the WHO *Genomic Sequencing of SARS-CoV-2* guidance ([World Health Organization, 2021](#ref-who-2021)) as the recommended base-quality standard for surveillance sequencing, and is consistent with PHA4GE quality control recommendations for pathogen genomic data ([Griffiths et al., 2022](#ref-griffiths-2022)).

**`min_mapped_reads: 50%`**
The WHO ([World Health Organization, 2021](#ref-who-2021)) and COG-UK sequencing standards recommend a minimum mapping rate of 70–80% for a high-quality sample. The pipeline uses 50% as the minimum fail gate — samples between 50–70% are flagged for review but not discarded, acknowledging that samples with lower viral load or partial degradation may still yield usable consensus sequences. Samples below 50% are failed at the post-mapping QC gate, consistent with PHA4GE guidance on minimum data quality for inclusion in genomic epidemiology analyses ([Griffiths et al., 2022](#ref-griffiths-2022)).

**`min_allele_frequency: 0.1`**
A 10% allele frequency threshold is standard for SARS-CoV-2 variant calling. The ARTIC amplicon sequencing pipeline ([Tyson et al., 2020](#ref-tyson-2020)) and Nextstrain analysis workflows apply a 5–10% minor allele frequency filter to distinguish true low-frequency variants from Illumina sequencing error, which has a typical substitution error rate of ~0.1–1% per base ([Schirmer et al., 2016](#ref-schirmer-2016)). A 10% threshold provides a conservative margin above the instrument noise floor while capturing genuine variants.

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
- [ ] `config/config.yaml` thresholds are appropriate for your samples
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
The `bms503-all` environment may not have been created yet. Re-run setup:
```bash
bash setup.sh
```
To force-recreate an existing but broken environment:
```bash
conda env remove -n bms503-all
bash setup.sh
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
bash setup.sh
./run_pipeline.sh
```

### High N content in consensus sequences
Low coverage or repetitive regions. Relax threshold if acceptable:
```yaml
qc_thresholds:
  max_n_content: 10
```

### Very few variants called
Ensure the reference genome is NC_045512.2 (Wuhan-Hu-1). If correct, the sample may be highly similar to the reference — check coverage depth in `results/mapping/*.coverage.txt`.

### Low mapping rate
Confirm the reference is NC_045512.2 (SARS-CoV-2 Wuhan-Hu-1):
```bash
head -n 1 data/reference/*.fasta          # Should show >NC_045512.2
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
- **Nextclade lineage assignment**: Uses the SARS-CoV-2 dataset downloaded by `setup.sh`. Pango lineage fields will be empty for samples that fail Nextclade's internal QC (e.g. very low coverage or highly divergent sequences).

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

---

## References

<a id="ref-aksamentov-2021"></a>Aksamentov, I., Roemer, C., Hodcroft, E. B., & Neher, R. A. (2021). Nextclade: Clade assignment, mutation calling and quality control for viral genomes. *Journal of Open Source Software*, *6*(67), Article 3773. https://doi.org/10.21105/joss.03773

<a id="ref-chen-2025"></a>Chen, S. (2025). fastp 1.0: An ultra-fast all-round tool for FASTQ data quality control and preprocessing. *iMeta*, *4*(5), Article e70078. https://doi.org/10.1002/imt2.70078

<a id="ref-cock-2009"></a>Cock, P. J. A., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., Friedberg, I., Hamelryck, T., Kauff, F., Wilczynski, B., & de Hoon, M. J. L. (2009). Biopython: Freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*, *25*(11), 1422–1423. https://doi.org/10.1093/bioinformatics/btp163

<a id="ref-coguk-2020"></a>COVID-19 Genomics UK (COG-UK) Consortium. (2020). An integrated national scale SARS-CoV-2 genomic surveillance network. *The Lancet Microbe*, *1*(3), e99–e100. https://doi.org/10.1016/S2666-5247(20)30054-9

<a id="ref-danecek-2021"></a>Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*, *10*(2), Article giab008. https://doi.org/10.1093/gigascience/giab008

<a id="ref-griffiths-2022"></a>Griffiths, E. J., Timme, R. E., Page, A. J., Alikhan, N.-F., Fornika, D., Maguire, F., Mendes, C. I., Tausch, S. H., Black, A., Connor, T. R., Tyson, G. H., Aanensen, D. M., Alcock, B., Campos, J., Christoffels, A., da Silva, A. G., Grunt, S., Haas, W., Hodcroft, E. B., & Hsiao, W. W. L. (2022). The PHA4GE SARS-CoV-2 contextual data specification for open genomic epidemiology. *GigaScience*, *11*, giac003. https://doi.org/10.1093/gigascience/giac003

<a id="ref-molder-2021"></a>Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., & Köster, J. (2021). Sustainable data analysis with Snakemake. *F1000Research*, *10*, Article 33. https://doi.org/10.12688/f1000research.29032.2

<a id="ref-pineiro-2020"></a>Piñeiro, C., Abuín, J. M., & Pichel, J. C. (2020). VeryFastTree: Speeding up the estimation of phylogenies for large alignments through parallelization and vectorization strategies. *Bioinformatics*, *36*(17), 4658–4659. https://doi.org/10.1093/bioinformatics/btaa582

<a id="ref-pineiro-2024"></a>Piñeiro, C., & Pichel, J. C. (2024). Efficient phylogenetic tree inference for massive taxonomic datasets: Harnessing the power of a server to analyze 1 million taxa. *GigaScience*, *13*, 1–12. https://doi.org/10.1093/gigascience/giae004

<a id="ref-schirmer-2016"></a>Schirmer, M., D'Amore, R., Ijaz, U. Z., Hall, N., & Quince, C. (2016). Illumina error profiles: Resolving fine-scale variation in metagenomic sequencing data. *BMC Bioinformatics*, *17*, Article 125. https://doi.org/10.1186/s12859-016-0984-y

<a id="ref-tyson-2020"></a>Tyson, J. R., James, P., Stoddart, D., Sparks, N., Sherrat, A., Noakes, C. J., & Slater, H. (2020). *Improvements to the ARTIC multiplex PCR method for SARS-CoV-2 genome sequencing using nanopore* [Preprint]. bioRxiv. https://doi.org/10.1101/2020.09.04.283077

<a id="ref-vasimuddin-2019"></a>Vasimuddin, M., Misra, S., Li, H., & Aluru, S. (2019). Efficient architecture-aware acceleration of BWA-MEM for multicore systems. *2019 IEEE International Parallel and Distributed Processing Symposium (IPDPS)*, 314–324. https://doi.org/10.1109/IPDPS.2019.00041

<a id="ref-who-2021"></a>World Health Organization. (2021). *Genomic sequencing of SARS-CoV-2: A guide to implementation for maximum impact on public health*. World Health Organization. https://www.who.int/publications/i/item/9789240018440
