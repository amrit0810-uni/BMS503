#!/bin/bash
#
# BMS503 Genomic Surveillance Pipeline - Main Entry Point
# Usage: ./run_pipeline.sh [--dry-run] [--cores N]
# First time on a new machine: bash setup.sh

set -euo pipefail

# Always run from the directory this script lives in
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PIPELINE_DIR"

# Clean invalid conda cache
if [ -d ".snakemake/conda" ]; then
    find .snakemake/conda -maxdepth 1 -type d -empty -delete
fi

# Configuration
DEFAULT_CORES=$(nproc)
DRY_RUN=false

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Functions
print_header() {
    echo -e "${GREEN}===============================================${NC}"
    echo -e "${GREEN}$1${NC}"
    echo -e "${GREEN}===============================================${NC}"
}

print_error() {
    echo -e "${RED}ERROR: $1${NC}" >&2
}

print_warning() {
    echo -e "${YELLOW}WARNING: $1${NC}"
}

print_info() {
    echo -e "${GREEN}INFO: $1${NC}"
}

show_usage() {
    cat << EOF
Usage: ./run_pipeline.sh [OPTIONS]

Options:
    --dry-run           Show what would be executed without running it
    --cores N           Limit CPU cores (default: all available via nproc)
    --help              Display this help message

Examples:
    ./run_pipeline.sh                  # Use all available cores (fastest)
    ./run_pipeline.sh --cores 8        # Limit to 8 cores (shared system)
    ./run_pipeline.sh --dry-run

First time on a new machine: bash setup.sh

EOF
}

# Parse command-line arguments
CORES=$DEFAULT_CORES
while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --cores)
            CORES="$2"
            shift 2
            ;;
        --help)
            show_usage
            exit 0
            ;;
        *)
            print_error "Unknown option: $1"
            show_usage
            exit 1
            ;;
    esac
done

# Print header
print_header "BMS503 Genomic Surveillance Pipeline ($CORES cores)"

# Record start time
START_TIME=$(date +%s)

# Check prerequisites
print_info "Checking prerequisites..."

CONDA_BASE=$(conda info --base 2>/dev/null || true)
if [ -z "$CONDA_BASE" ]; then
    print_error "conda not found. Run bash setup.sh first."
    exit 1
fi

ENV_DIR="${CONDA_BASE}/envs/bms503-all"
if [ ! -d "$ENV_DIR" ]; then
    print_error "bms503-all environment not found. Run bash setup.sh first."
    exit 1
fi

if [ ! -f "Snakefile" ]; then
    print_error "Snakefile not found. Make sure you're in the pipeline root directory."
    exit 1
fi

# Ensure required directories exist
mkdir -p data/raw data/reference data/reference_db logs

# bcftools (bioconda) is compiled against libgsl.so.25 but conda installs GSL 2.7.
# Symlink bridges the version gap; LD_PRELOAD supplies the cblas symbols.
ENV_LIB="${ENV_DIR}/lib"
if [ -f "${ENV_LIB}/libgsl.so.27" ] && [ ! -e "${ENV_LIB}/libgsl.so.25" ]; then
    ln -sf "${ENV_LIB}/libgsl.so.27" "${ENV_LIB}/libgsl.so.25"
fi
if [ -f "${ENV_LIB}/libgslcblas.so.0" ]; then
    export LD_PRELOAD="${ENV_LIB}/libgslcblas.so.0:${LD_PRELOAD:-}"
fi

# Check data availability
DATA_COUNT=$(find data/raw -type f \( -name "*.fastq*" -o -name "*.fq*" \) 2>/dev/null | wc -l)

if [ "$DATA_COUNT" -eq 0 ]; then
    print_warning "No FASTQ files found in data/raw/"
    print_info "To run the pipeline, copy your FASTQ files to data/raw/"
    print_info "Supported formats: .fastq.gz, .fastq, .fq.gz, .fq"
    echo ""
else
    print_info "Found $DATA_COUNT FASTQ file(s) in data/raw/"
fi

# Check reference availability
REF_COUNT=$(find data/reference -type f \( -name "*.fasta" -o -name "*.fa" -o -name "*.fna" -o -name "*.ffn" \) 2>/dev/null | wc -l)

if [ "$REF_COUNT" -eq 0 ]; then
    print_warning "No reference genome found in data/reference/"
    print_info "To run the pipeline, copy your reference genome to data/reference/"
    print_info "Supported formats: .fasta, .fa, .fna, .ffn"
    echo ""
else
    REF_FILE=$(find data/reference -type f \( -name "*.fasta" -o -name "*.fa" -o -name "*.fna" -o -name "*.ffn" \) | head -1)
    print_info "Auto-detected reference genome: $(basename $REF_FILE)"
fi

# Verify organism is configured
ORGANISM=$(grep "^organism:" config/config.yaml | head -1 | cut -d'"' -f2)
if [ -z "$ORGANISM" ]; then
    print_warning "Organism not specified in config/config.yaml. Using default name."
else
    print_info "Organism: $ORGANISM"
fi

# Exit if critical files missing
if [ "$DATA_COUNT" -eq 0 ] || [ "$REF_COUNT" -eq 0 ]; then
    print_error "Cannot proceed: Missing data or reference files"
    echo ""
    print_info "Setup instructions:"
    echo "  1. Copy FASTQ files to data/raw/"
    echo "  2. Copy reference genome to data/reference/"
    echo "  3. Run ./run_pipeline.sh again"
    exit 1
fi

SNAKEFILE="Snakefile"

# Build Snakemake command
SNAKEMAKE_CMD="snakemake -s $SNAKEFILE --cores $CORES --configfile config/config.yaml"

if [ "$DRY_RUN" = true ]; then
    SNAKEMAKE_CMD="$SNAKEMAKE_CMD --dry-run"
    print_info "Running in DRY-RUN mode (no commands will be executed)"
else
    SNAKEMAKE_CMD="$SNAKEMAKE_CMD -p"
fi

# Execute pipeline
print_info "Starting pipeline with $CORES cores..."
print_info "Command: $SNAKEMAKE_CMD"
echo ""

export PATH="$(conda info --base)/envs/bms503-all/bin:$PATH"

if eval "$SNAKEMAKE_CMD"; then
    # Calculate elapsed time
    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))
    HOURS=$((ELAPSED / 3600))
    MINUTES=$(((ELAPSED % 3600) / 60))
    SECONDS=$((ELAPSED % 60))

    print_header "Pipeline Completed Successfully"
    echo ""
    print_info "Results are available in the following directories:"
    echo "  - Quality Control:    results/qc/"
    echo "  - Mapping Results:    results/mapping/"
    echo "  - Variant Calls:      results/variants/"
    echo "  - Phylogenetic Trees: results/phylogeny/"
    echo "  - Reports:            results/reports/"
    echo ""
    print_info "View the diagnostic report: results/reports/diagnostic_report.html"
    echo ""
    printf "${GREEN}Total Execution Time: ${HOURS}h ${MINUTES}m ${SECONDS}s${NC}\n"
else
    print_error "Pipeline execution failed. Check error messages above."
    exit 1
fi
