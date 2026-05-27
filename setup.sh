#!/bin/bash
#
# BMS503 Genomic Surveillance Pipeline — First-Time Setup
# Run once after extracting the zip:  bash setup.sh
#
# This script:
#   1. Makes all pipeline scripts executable
#   2. Creates the bms503-all conda environment (skips if already exists)
#   3. Applies a runtime fix for a bcftools/GSL library mismatch
#   4. Creates required data directories
#   5. Downloads the Nextclade SARS-CoV-2 dataset (skips if already present)
#
# After this completes, run the pipeline with:  ./run_pipeline.sh

set -euo pipefail

# Always run from the directory this script lives in
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# ── Colours ────────────────────────────────────────────────────────────────────
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BOLD='\033[1m'
NC='\033[0m'

print_ok()   { echo -e "${GREEN}  ✔  $1${NC}"; }
print_info() { echo -e "${GREEN}  →  $1${NC}"; }
print_warn() { echo -e "${YELLOW}  !  $1${NC}"; }
print_err()  { echo -e "${RED}  ✘  $1${NC}" >&2; }

echo ""
echo -e "${BOLD}BMS503 Genomic Surveillance Pipeline — Setup${NC}"
echo "────────────────────────────────────────────"
echo ""

# ── Step 1: Fix script permissions ─────────────────────────────────────────────
print_info "Setting script permissions..."
chmod +x "$SCRIPT_DIR"/*.sh
print_ok "All .sh files are now executable"
echo ""

# ── Step 2: Check conda ─────────────────────────────────────────────────────────
if ! command -v conda &> /dev/null; then
    print_err "conda not found. Install Miniconda first:"
    echo ""
    echo "  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    echo "  bash Miniconda3-latest-Linux-x86_64.sh"
    echo "  source ~/.bashrc"
    echo ""
    exit 1
fi

CONDA_BASE=$(conda info --base)
print_ok "conda found: $CONDA_BASE"

# Require mamba for reliable environment creation
if ! command -v mamba &> /dev/null; then
    print_err "mamba not found. Install it first:"
    echo ""
    echo "  conda install -n base -c conda-forge mamba"
    echo ""
    echo "Then re-run:  bash setup.sh"
    exit 1
fi
PKG_CMD="mamba"
print_ok "mamba found"
echo ""

# ── Step 3: Create conda environment ───────────────────────────────────────────
ENV_DIR="${CONDA_BASE}/envs/bms503-all"

if [ -d "$ENV_DIR" ]; then
    print_ok "bms503-all environment already exists — skipping creation"
else
    print_info "Creating bms503-all environment (this takes 5–10 minutes)..."
    echo ""
    $PKG_CMD env create -f config/env_all.yml
    echo ""
    print_ok "bms503-all environment created"
fi
echo ""

# ── Step 4: Apply bcftools/GSL library fix ──────────────────────────────────────
# bcftools (bioconda) is compiled against libgsl.so.25 but conda installs GSL 2.7
# (libgsl.so.27). A symlink bridges the gap; LD_PRELOAD supplies the cblas symbols.
ENV_LIB="${ENV_DIR}/lib"
if [ -f "${ENV_LIB}/libgsl.so.27" ] && [ ! -e "${ENV_LIB}/libgsl.so.25" ]; then
    ln -sf "${ENV_LIB}/libgsl.so.27" "${ENV_LIB}/libgsl.so.25"
    print_ok "Applied bcftools/GSL compatibility fix"
else
    print_ok "bcftools/GSL compatibility — no action needed"
fi
echo ""

# ── Step 5: Create required directories ────────────────────────────────────────
print_info "Creating data directories..."
mkdir -p data/raw data/reference data/reference_db data/nextclade_db logs
print_ok "data/raw/         ← place FASTQ files here"
print_ok "data/reference/   ← place reference genome here"
print_ok "data/reference_db/ (pipeline database — do not edit)"
print_ok "logs/             (pipeline logs)"
echo ""

# ── Step 6: Download Nextclade SARS-CoV-2 dataset ──────────────────────────────
NEXTCLADE_DB="data/nextclade_db/sars-cov-2"
if [ -f "${NEXTCLADE_DB}/pathogen.json" ]; then
    print_ok "Nextclade SARS-CoV-2 dataset already present — skipping download"
else
    print_info "Downloading Nextclade SARS-CoV-2 dataset (requires internet)..."
    "${ENV_DIR}/bin/nextclade" dataset get \
        --name sars-cov-2 \
        --output-dir "$NEXTCLADE_DB" 2>&1 | tail -3
    if [ -f "${NEXTCLADE_DB}/pathogen.json" ]; then
        print_ok "Nextclade dataset downloaded to ${NEXTCLADE_DB}"
    else
        print_warn "Nextclade dataset download failed — pipeline will download it at runtime instead"
    fi
fi
echo ""

# ── Done ───────────────────────────────────────────────────────────────────────
echo -e "${BOLD}${GREEN}Setup complete.${NC}"
echo ""
echo "Next steps:"
echo "  1. Copy FASTQ files to       data/raw/"
echo "  2. Copy reference genome to  data/reference/"
echo "  3. Run the pipeline:         ./run_pipeline.sh"
echo ""
echo "See QUICKSTART.md for more details."
echo ""
