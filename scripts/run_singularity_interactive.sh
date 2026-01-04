#!/bin/bash
# =============================================================================
# Interactive R Session with Singularity
# =============================================================================
#
# This script starts an interactive R session inside the Singularity container
# with your project directories mounted.
#
# Author: Ji (gynecoloji)
# Date: 2025-01-03
#
# Usage:
#   bash scripts/run_singularity_interactive.sh
#
# =============================================================================

set -e

# Configuration
SIF_FILE="containers/bioinformatic_r_4.5.1.sif"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

echo "========================================================================"
echo "           Interactive R Session - Singularity Container"
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

# Determine Singularity command
if command -v singularity &> /dev/null; then
    SING_CMD="singularity"
elif command -v apptainer &> /dev/null; then
    SING_CMD="apptainer"
else
    echo "ERROR: Neither Singularity nor Apptainer found"
    exit 1
fi

echo "Container: $SIF_FILE"
echo "Command:   $SING_CMD"
echo ""
echo "Working directory will be mounted at: /workspace"
echo ""
echo "You can access your files as:"
echo "  data/           -> /workspace/data/"
echo "  results/        -> /workspace/results/"
echo "  seurat_pipeline/ -> /workspace/seurat_pipeline/"
echo ""
echo "Useful commands inside R:"
echo "  # Source pipeline functions"
echo "  source('/workspace/seurat_pipeline/R/utils.R')"
echo "  source('/workspace/seurat_pipeline/R/qc_functions.R')"
echo "  source('/workspace/seurat_pipeline/R/processing_functions.R')"
echo "  "
echo "  # Load your data"
echo "  seurat_obj <- readRDS('/workspace/data/your_data.rds')"
echo "  "
echo "  # Use pipeline functions interactively"
echo "  seurat_obj <- add_qc_metrics(seurat_obj, species='human')"
echo ""
echo "Press Ctrl+D or type 'quit()' to exit R"
echo ""
echo "========================================================================"
echo ""

# Get absolute path
ABS_PROJECT=$(realpath "$PROJECT_ROOT")

# Start interactive R session
$SING_CMD exec \
    --bind "$ABS_PROJECT:/workspace" \
    --pwd /workspace \
    "$SIF_FILE" \
    R --no-save --quiet
