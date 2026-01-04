#!/bin/bash
# =============================================================================
# Build Singularity Container from Docker Hub
# =============================================================================
#
# This script pulls your existing Docker image and converts it to Singularity
# format for use on HPC systems.
#
# Author: Ji (gynecoloji)
# Date: 2025-01-03
#
# Usage:
#   bash scripts/build_singularity.sh
#
# =============================================================================

set -e  # Exit on error

# Configuration
DOCKER_IMAGE="docker://gynecoloji/bioinformatic_r_4_5_1:v2"
SIF_FILE="containers/bioinformatic_r_4.5.1.sif"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo -e "${BLUE}======================================================================${NC}"
echo -e "${BLUE}       Building Singularity Container for Seurat Pipeline${NC}"
echo -e "${BLUE}======================================================================${NC}"
echo ""

# Create containers directory
mkdir -p "$(dirname "$SIF_FILE")"

echo -e "${GREEN}[1/3] Checking Singularity installation...${NC}"
if ! command -v singularity &> /dev/null; then
    echo -e "${YELLOW}Singularity not found. Checking for Apptainer...${NC}"
    if command -v apptainer &> /dev/null; then
        SING_CMD="apptainer"
        echo -e "${GREEN}Using Apptainer (Singularity)${NC}"
    else
        echo "ERROR: Neither Singularity nor Apptainer found."
        echo "Please install Singularity/Apptainer first."
        exit 1
    fi
else
    SING_CMD="singularity"
    echo -e "${GREEN}Using Singularity${NC}"
fi

echo ""
echo -e "${GREEN}[2/3] Pulling Docker image from Docker Hub...${NC}"
echo "Source: $DOCKER_IMAGE"
echo "Output: $SIF_FILE"
echo ""
echo "This may take several minutes (downloading ~2.5 GB)..."
echo ""

# Pull and convert
$SING_CMD pull "$SIF_FILE" "$DOCKER_IMAGE"

echo ""
echo -e "${GREEN}[3/3] Verifying container...${NC}"

# Check file size
FILE_SIZE=$(du -h "$SIF_FILE" | cut -f1)
echo "Container size: $FILE_SIZE"

# Test R version
echo ""
echo "Testing R installation..."
$SING_CMD exec "$SIF_FILE" R --version | head -n 1

# Test Seurat availability
echo ""
echo "Testing Seurat package..."
$SING_CMD exec "$SIF_FILE" Rscript -e "library(Seurat); cat('Seurat version:', as.character(packageVersion('Seurat')), '\n')"

echo ""
echo -e "${GREEN}======================================================================${NC}"
echo -e "${GREEN}                   Build Complete!${NC}"
echo -e "${GREEN}======================================================================${NC}"
echo ""
echo "Container file: $SIF_FILE"
echo ""
echo "Quick start:"
echo "  # Interactive R session"
echo "  singularity shell $SIF_FILE"
echo ""
echo "  # Run pipeline"
echo "  singularity exec $SIF_FILE \\"
echo "    Rscript seurat_pipeline/run_seurat_pipeline.R \\"
echo "    --config seurat_pipeline/config/default_config.yaml"
echo ""
echo "See scripts/run_singularity_*.sh for more examples"
echo ""
