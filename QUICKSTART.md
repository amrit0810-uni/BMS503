# BMS503 Genomic Surveillance Pipeline — Quick Start

## First time on a new machine

```bash
# Step 1 — one-time setup (creates conda environment + directories)
bash setup.sh

# Step 2 — add your data
cp /path/to/your/samples/*.fastq.gz  data/raw/
cp /path/to/reference.fasta          data/reference/

# Step 3 — run
./run_pipeline.sh
```

That's it. Results appear in `results/reports/diagnostic_report.html`.

---

## Every run after that

```bash
./run_pipeline.sh
```

Snakemake skips everything already done — only new or changed files are processed.

---

## Adding new samples mid-project

Drop new FASTQ files into `data/raw/` and re-run:

```bash
cp /path/to/new_sample_1.fastq.gz data/raw/
cp /path/to/new_sample_2.fastq.gz data/raw/
./run_pipeline.sh
```

New samples are flagged `[NEW]` in the report and receive focused interpretation.

---

## Useful options

| Command | What it does |
|---------|-------------|
| `./run_pipeline.sh` | Run using all available CPU cores |
| `./run_pipeline.sh --cores 8` | Limit to 8 cores (shared systems) |
| `./run_pipeline.sh --dry-run` | Preview steps without executing |
| `./clean.sh` | Remove all outputs for a clean transfer |
| `bash setup.sh` | Re-run setup (safe to repeat) |

---

## Where results are saved

| Directory | Contents |
|-----------|----------|
| `results/qc/` | fastp HTML reports, QC summary |
| `results/mapping/` | BAM files, flagstat, coverage |
| `results/variants/` | VCF files per sample |
| `results/phylogeny/` | Alignment, Newick tree, Nextclade TSV |
| `results/reports/` | **Diagnostic report (HTML + TXT)** |

---

## Requirements

- Linux / WSL (Ubuntu 20.04+)
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Anaconda
- 5–10 GB free disk space for the conda environment
- 8–16 GB RAM recommended for large sample sets
