#!/usr/bin/env Rscript
# =============================================================================
# Example Seurat Pipeline Workflow
# =============================================================================
#
# This script demonstrates how to use the Seurat pipeline functions
# for a complete single-cell RNA-seq analysis workflow.
#
# Author: Ji (gynecoloji)
# Date: 2025-01-03
# =============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# Source pipeline functions
source("R/utils.R")
source("R/qc_functions.R")
source("R/processing_functions.R")
source("R/plotting_functions.R")

# Set seed for reproducibility
set.seed(42)

# =============================================================================
# EXAMPLE 1: Create Test Data
# =============================================================================

cat("Creating test dataset...\n")

# Create a simple test dataset (replace with your real data)
create_test_data <- function(n_cells = 1000, n_genes = 2000) {
  # Simulate count matrix
  counts <- matrix(
    rpois(n_cells * n_genes, lambda = 5),
    nrow = n_genes,
    ncol = n_cells
  )
  
  # Add gene names
  rownames(counts) <- paste0("Gene_", 1:n_genes)
  
  # Add cell barcodes
  colnames(counts) <- paste0("Cell_", 1:n_cells)
  
  # Add some mitochondrial genes
  mt_genes <- paste0("MT-", 1:20)
  rownames(counts)[1:20] <- mt_genes
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = counts,
    project = "TestProject",
    min.cells = 0,
    min.features = 0
  )
  
  # Add batch information for testing
  seurat_obj$batch <- sample(c("Batch1", "Batch2"), n_cells, replace = TRUE)
  
  return(seurat_obj)
}

# Create test object
test_obj <- create_test_data(n_cells = 1000, n_genes = 2000)

# Save test data
dir.create("data", showWarnings = FALSE)
saveRDS(test_obj, "data/test_data.rds")
cat("Test data created and saved to: data/test_data.rds\n\n")

# =============================================================================
# EXAMPLE 2: Manual Step-by-Step Analysis
# =============================================================================

cat("=== Running step-by-step analysis ===\n\n")

# Create output directories
output_dir <- "results/example"
create_output_dir(output_dir, verbose = TRUE)
create_output_dir(file.path(output_dir, "qc"), verbose = FALSE)
create_output_dir(file.path(output_dir, "plots"), verbose = FALSE)

# -----------------------------------------------------------------------------
# Step 1: Add QC metrics
# -----------------------------------------------------------------------------

cat("Step 1: Adding QC metrics...\n")
test_obj <- add_qc_metrics(
  seurat_obj = test_obj,
  species = "human",
  verbose = TRUE
)

# -----------------------------------------------------------------------------
# Step 2: Plot QC metrics
# -----------------------------------------------------------------------------

cat("\nStep 2: Plotting QC metrics...\n")
qc_plots <- plot_qc_metrics(
  seurat_obj = test_obj,
  sample_name = "TestProject",
  output_dir = file.path(output_dir, "qc"),
  qc_thresholds = list(
    min_umi = 50,
    min_genes = 50,
    max_mito = 25,
    min_complexity = 0.5
  ),
  verbose = TRUE
)

# -----------------------------------------------------------------------------
# Step 3: Filter cells and genes
# -----------------------------------------------------------------------------

cat("\nStep 3: Filtering cells and genes...\n")
test_filtered <- filter_cells_and_genes(
  seurat_obj = test_obj,
  sample_name = "TestProject",
  output_dir = file.path(output_dir, "qc"),
  min_umi = 50,
  min_genes = 50,
  max_mito = 25,
  min_complexity = 0.5,
  min_cells_per_gene = 5,
  save_stats = TRUE,
  verbose = TRUE
)

# -----------------------------------------------------------------------------
# Step 4: Process the data
# -----------------------------------------------------------------------------

cat("\nStep 4: Processing data...\n")

# Create processing config
proc_config <- list(
  n_features = 500,
  dims_use = 10,
  resolution = 0.5,
  n_cores = 1,
  run_umap = TRUE,
  run_tsne = FALSE,
  run_jackstraw = FALSE
)

test_processed <- process_seurat_object(
  seurat_obj = test_filtered,
  config = proc_config,
  species = "human",
  verbose = TRUE
)

# -----------------------------------------------------------------------------
# Step 5: Visualize results
# -----------------------------------------------------------------------------

cat("\nStep 5: Creating visualizations...\n")

# Elbow plot
plot_elbow(
  seurat_obj = test_processed,
  ndims = 20,
  output_dir = file.path(output_dir, "plots"),
  sample_name = "TestProject",
  save_plot = TRUE,
  verbose = TRUE
)

# Cluster QC
plot_cluster_qc(
  seurat_obj = test_processed,
  output_dir = file.path(output_dir, "plots"),
  sample_name = "TestProject",
  save_plots = TRUE,
  verbose = TRUE
)

# UMAP plot
p_umap <- DimPlot(
  test_processed,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE
)

ggsave(
  filename = file.path(output_dir, "plots", "TestProject_umap.pdf"),
  plot = p_umap,
  width = 8,
  height = 7
)

# -----------------------------------------------------------------------------
# Step 6: Save results
# -----------------------------------------------------------------------------

cat("\nStep 6: Saving results...\n")

# Save processed object
save_seurat_object(
  seurat_obj = test_processed,
  file_path = file.path(output_dir, "TestProject_final.rds"),
  compress = "gzip",
  verbose = TRUE
)

# Save processing summary
summary_df <- get_processing_summary(test_processed)
write.csv(
  summary_df,
  file.path(output_dir, "TestProject_summary.csv"),
  row.names = FALSE
)

cat("Summary saved to:", file.path(output_dir, "TestProject_summary.csv"), "\n")

# =============================================================================
# EXAMPLE 3: Batch Correction Workflow
# =============================================================================

cat("\n=== Testing batch correction ===\n\n")

# Create batch correction config
batch_config <- list(
  method = "harmony",
  dims = 10,
  resolution = 0.5,
  n_cores = 1,
  run_umap = TRUE,
  run_tsne = FALSE
)

# Run batch correction
test_batch_corrected <- process_batch_correction(
  seurat_obj = test_filtered,
  batch_var = "batch",
  config = batch_config,
  species = "human",
  verbose = TRUE
)

# Plot batch correction results
batch_plots <- plot_batch_correction_results(
  seurat_obj = test_batch_corrected,
  batch_var = "batch",
  output_dir = file.path(output_dir, "plots"),
  sample_name = "TestProject_batch",
  reduction = "umap",
  save_plots = TRUE,
  verbose = TRUE
)

# =============================================================================
# EXAMPLE 4: Using Configuration Files
# =============================================================================

cat("\n=== Testing configuration-based workflow ===\n\n")

# Load configuration
config <- load_config("config/default_config.yaml", verbose = TRUE)

# Modify for test data
config$qc$min_umi <- 50
config$qc$min_genes <- 50
config$processing$n_features <- 500
config$processing$dims_use <- 10

# Save modified config
modified_config_file <- "config/test_config.yaml"
save_config(
  config = config,
  config_file = modified_config_file,
  verbose = TRUE
)

cat("Modified configuration saved to:", modified_config_file, "\n")

# =============================================================================
# COMPLETE
# =============================================================================

cat("\n=============================================================================\n")
cat("                    Example Workflow Completed Successfully!                 \n")
cat("=============================================================================\n\n")

cat("Results saved to:", output_dir, "\n")
cat("\nGenerated files:\n")
cat("  - QC plots:", file.path(output_dir, "qc"), "\n")
cat("  - Analysis plots:", file.path(output_dir, "plots"), "\n")
cat("  - Processed object:", file.path(output_dir, "TestProject_final.rds"), "\n")
cat("  - Summary statistics:", file.path(output_dir, "TestProject_summary.csv"), "\n")
cat("\nYou can now:\n")
cat("  1. Examine the QC plots\n")
cat("  2. Review the cluster assignments\n")
cat("  3. Load the processed object for further analysis\n")
cat("\n")

# Print session info
cat("Session Information:\n")
cat("===================\n")
print(sessionInfo())
