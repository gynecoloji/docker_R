#!/bin/bash
# =============================================================================
# Run Seurat Pipeline with Singularity
# =============================================================================
#
# This script runs the Seurat pipeline using your existing R container
# pulled from Docker Hub.
#
# Author: Ji (gynecoloji)
# Date: 2025-01-03
#
# Usage:
#   bash scripts/run_singularity_local.sh [INPUT_FILE] [OUTPUT_DIR]
#
# Example:
#   bash scripts/run_singularity_local.sh data/pbmc.rds results/pbmc_analysis
#
# =============================================================================

set -e  # Exit on error

# =============================================================================
# CONFIGURATION
# =============================================================================

# Singularity container
SIF_FILE="containers/bioinformatic_r_4.5.1.sif"

# Default values
INPUT_FILE="${1:-data/test_data.rds}"
OUTPUT_DIR="${2:-results/analysis}"
CONFIG_FILE="${3:-seurat_pipeline/config/default_config.yaml}"
CORES="${4:-4}"

# Project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# =============================================================================
# VALIDATE INPUTS
# =============================================================================

echo "========================================================================"
echo "           Seurat Pipeline - Singularity Execution"
echo "========================================================================"
echo ""

# Check if container exists
if [ ! -f "$SIF_FILE" ]; then
    echo "ERROR: Singularity container not found: $SIF_FILE"
    echo ""
    echo "Please build the container first:"
    echo "  bash scripts/build_singularity.sh"
    echo ""
    exit 1
fi

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: Input file not found: $INPUT_FILE"
    echo ""
    echo "Usage: bash scripts/run_singularity_local.sh INPUT_FILE [OUTPUT_DIR]"
    echo ""
    exit 1
fi

# Check if config exists
if [ ! -f "$CONFIG_FILE" ]; then
    echo "ERROR: Config file not found: $CONFIG_FILE"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"
mkdir -p logs

echo "Configuration:"
echo "  Container:  $SIF_FILE"
echo "  Input:      $INPUT_FILE"
echo "  Output:     $OUTPUT_DIR"
echo "  Config:     $CONFIG_FILE"
echo "  Cores:      $CORES"
echo ""

# =============================================================================
# DETERMINE SINGULARITY COMMAND
# =============================================================================

if command -v singularity &> /dev/null; then
    SING_CMD="singularity"
elif command -v apptainer &> /dev/null; then
    SING_CMD="apptainer"
else
    echo "ERROR: Neither Singularity nor Apptainer found"
    exit 1
fi

echo "Using: $SING_CMD"
echo ""

# =============================================================================
# RUN PIPELINE
# =============================================================================

echo "Starting pipeline..."
echo "========================================================================"
echo ""

# Get absolute paths
ABS_INPUT=$(realpath "$INPUT_FILE")
ABS_OUTPUT=$(realpath "$OUTPUT_DIR")
ABS_CONFIG=$(realpath "$CONFIG_FILE")
ABS_PROJECT=$(realpath "$PROJECT_ROOT")

# Run with Singularity
$SING_CMD exec \
    --bind "$ABS_PROJECT:/workspace" \
    --pwd /workspace \
    "$SIF_FILE" \
    Rscript /workspace/seurat_pipeline/run_seurat_pipeline.R \
        --config "/workspace/$CONFIG_FILE" \
        --input "/workspace/$INPUT_FILE" \
        --output "/workspace/$OUTPUT_DIR" \
        --cores "$CORES"

EXIT_CODE=$?

# =============================================================================
# REPORT RESULTS
# =============================================================================

echo ""
echo "========================================================================"
if [ $EXIT_CODE -eq 0 ]; then
    echo "                    Pipeline Completed Successfully!"
    echo "========================================================================"
    echo ""
    echo "Results saved to: $OUTPUT_DIR"
    echo ""
    echo "Generated files:"
    echo "  - QC plots:      $OUTPUT_DIR/qc/"
    echo "  - Visualizations: $OUTPUT_DIR/plots/"
    echo "  - Objects:       $OUTPUT_DIR/objects/"
    echo "  - Logs:          logs/pipeline.log"
    echo ""
else
    echo "                    Pipeline Failed (Exit code: $EXIT_CODE)"
    echo "========================================================================"
    echo ""
    echo "Check logs for details: logs/pipeline.log"
    echo ""
fi

exit $EXIT_CODE
