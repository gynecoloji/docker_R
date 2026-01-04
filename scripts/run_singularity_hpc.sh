#!/bin/bash
#SBATCH --job-name=seurat_pipeline
#SBATCH --output=logs/slurm_%j.out
#SBATCH --error=logs/slurm_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --partition=general
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@example.com

# =============================================================================
# SLURM Job Script for Seurat Pipeline
# =============================================================================
#
# This script runs the Seurat pipeline on an HPC cluster using SLURM and
# Singularity with your existing Docker Hub image.
#
# Author: Ji (gynecoloji)
# Date: 2025-01-03
#
# Usage:
#   sbatch scripts/run_singularity_hpc.sh
#
# To customize:
#   1. Edit SBATCH parameters above (time, CPUs, memory, partition)
#   2. Edit CONFIGURATION section below (input file, output dir, etc.)
#   3. Submit: sbatch scripts/run_singularity_hpc.sh
#
# =============================================================================

# Exit on error
set -e

# =============================================================================
# CONFIGURATION
# =============================================================================

# Container location (adjust if needed)
SIF_FILE="containers/bioinformatic_r_4.5.1.sif"

# Input/Output (EDIT THESE FOR YOUR ANALYSIS)
INPUT_FILE="data/my_dataset.rds"
OUTPUT_DIR="results/my_analysis_${SLURM_JOB_ID}"
CONFIG_FILE="seurat_pipeline/config/default_config.yaml"

# Processing parameters
CORES=${SLURM_CPUS_PER_TASK:-20}

# Project directory (automatically detected)
PROJECT_DIR="$SLURM_SUBMIT_DIR"

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================

echo "========================================================================"
echo "              Seurat Pipeline - HPC SLURM Execution"
echo "========================================================================"
echo ""
echo "SLURM Job Information:"
echo "  Job ID:        $SLURM_JOB_ID"
echo "  Job Name:      $SLURM_JOB_NAME"
echo "  Node:          $SLURM_NODELIST"
echo "  CPUs:          $SLURM_CPUS_PER_TASK"
echo "  Memory:        ${SLURM_MEM_PER_NODE}MB"
echo "  Partition:     $SLURM_JOB_PARTITION"
echo "  Submit Dir:    $SLURM_SUBMIT_DIR"
echo ""

# Load Singularity module (adjust for your HPC)
if command -v module &> /dev/null; then
    echo "Loading Singularity module..."
    # Uncomment and adjust for your HPC system:
    # module load singularity
    # module load apptainer
fi

# Determine Singularity command
if command -v singularity &> /dev/null; then
    SING_CMD="singularity"
    echo "Using: singularity ($(singularity --version))"
elif command -v apptainer &> /dev/null; then
    SING_CMD="apptainer"
    echo "Using: apptainer ($(apptainer --version))"
else
    echo "ERROR: Neither Singularity nor Apptainer found"
    echo "Please load the appropriate module first"
    exit 1
fi

echo ""

# =============================================================================
# VALIDATE INPUTS
# =============================================================================

echo "Configuration:"
echo "  Project Dir:   $PROJECT_DIR"
echo "  Container:     $SIF_FILE"
echo "  Input:         $INPUT_FILE"
echo "  Output:        $OUTPUT_DIR"
echo "  Config:        $CONFIG_FILE"
echo "  Cores:         $CORES"
echo ""

# Change to project directory
cd "$PROJECT_DIR"

# Check if container exists
if [ ! -f "$SIF_FILE" ]; then
    echo "ERROR: Singularity container not found: $SIF_FILE"
    echo ""
    echo "Build the container first:"
    echo "  bash scripts/build_singularity.sh"
    echo ""
    exit 1
fi

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: Input file not found: $INPUT_FILE"
    exit 1
fi

# Check if config exists
if [ ! -f "$CONFIG_FILE" ]; then
    echo "ERROR: Config file not found: $CONFIG_FILE"
    exit 1
fi

# Create output directories
mkdir -p "$OUTPUT_DIR"
mkdir -p logs

echo "All files validated successfully."
echo ""

# =============================================================================
# SET ENVIRONMENT VARIABLES
# =============================================================================

# R environment variables
export OMP_NUM_THREADS=$CORES
export R_MAX_VSIZE=100Gb

# Singularity environment variables
export SINGULARITY_BIND="$PROJECT_DIR:/workspace"
export SINGULARITYENV_OMP_NUM_THREADS=$CORES

echo "Environment variables set."
echo ""

# =============================================================================
# RUN PIPELINE
# =============================================================================

echo "========================================================================"
echo "Starting Seurat pipeline..."
echo "Start time: $(date)"
echo "========================================================================"
echo ""

# Record start time
START_TIME=$(date +%s)

# Run the pipeline
$SING_CMD exec \
    --bind "$PROJECT_DIR:/workspace" \
    --pwd /workspace \
    --cleanenv \
    "$SIF_FILE" \
    Rscript /workspace/seurat_pipeline/run_seurat_pipeline.R \
        --config "/workspace/$CONFIG_FILE" \
        --input "/workspace/$INPUT_FILE" \
        --output "/workspace/$OUTPUT_DIR" \
        --cores "$CORES"

# Capture exit code
EXIT_CODE=$?

# Record end time
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
HOURS=$((RUNTIME / 3600))
MINUTES=$(((RUNTIME % 3600) / 60))
SECONDS=$((RUNTIME % 60))

# =============================================================================
# REPORT RESULTS
# =============================================================================

echo ""
echo "========================================================================"
echo "End time: $(date)"
echo "Runtime: ${HOURS}h ${MINUTES}m ${SECONDS}s"
echo "========================================================================"
echo ""

if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ Pipeline completed successfully!"
    echo ""
    echo "Results saved to: $OUTPUT_DIR"
    echo ""
    echo "Generated files:"
    echo "  - QC plots:       $OUTPUT_DIR/qc/"
    echo "  - Visualizations: $OUTPUT_DIR/plots/"
    echo "  - Objects:        $OUTPUT_DIR/objects/"
    echo "  - Logs:           logs/pipeline.log"
    echo "  - SLURM output:   logs/slurm_${SLURM_JOB_ID}.out"
    echo ""
    
    # Copy important files to a summary location
    SUMMARY_DIR="$OUTPUT_DIR/summary"
    mkdir -p "$SUMMARY_DIR"
    
    echo "Creating job summary..."
    cat > "$SUMMARY_DIR/job_info.txt" << EOF
SLURM Job Summary
=================
Job ID:       $SLURM_JOB_ID
Job Name:     $SLURM_JOB_NAME
Node:         $SLURM_NODELIST
CPUs:         $SLURM_CPUS_PER_TASK
Memory:       ${SLURM_MEM_PER_NODE}MB
Partition:    $SLURM_JOB_PARTITION
Start Time:   $(date -d @$START_TIME)
End Time:     $(date -d @$END_TIME)
Runtime:      ${HOURS}h ${MINUTES}m ${SECONDS}s
Exit Code:    $EXIT_CODE

Input File:   $INPUT_FILE
Output Dir:   $OUTPUT_DIR
Config File:  $CONFIG_FILE
Container:    $SIF_FILE

Status:       SUCCESS
EOF
    
    echo "Job summary saved to: $SUMMARY_DIR/job_info.txt"
    
else
    echo "✗ Pipeline failed with exit code: $EXIT_CODE"
    echo ""
    echo "Check logs for details:"
    echo "  - Pipeline log: logs/pipeline.log"
    echo "  - SLURM output: logs/slurm_${SLURM_JOB_ID}.out"
    echo "  - SLURM error:  logs/slurm_${SLURM_JOB_ID}.err"
    echo ""
    
    # Create failure summary
    SUMMARY_DIR="$OUTPUT_DIR/summary"
    mkdir -p "$SUMMARY_DIR"
    
    cat > "$SUMMARY_DIR/job_info.txt" << EOF
SLURM Job Summary
=================
Job ID:       $SLURM_JOB_ID
Job Name:     $SLURM_JOB_NAME
Node:         $SLURM_NODELIST
Runtime:      ${HOURS}h ${MINUTES}m ${SECONDS}s
Exit Code:    $EXIT_CODE

Status:       FAILED

Check logs:
  - logs/pipeline.log
  - logs/slurm_${SLURM_JOB_ID}.out
  - logs/slurm_${SLURM_JOB_ID}.err
EOF
fi

echo ""
echo "========================================================================"

exit $EXIT_CODE
