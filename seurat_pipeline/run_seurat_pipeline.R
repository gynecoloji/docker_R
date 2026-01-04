#!/usr/bin/env Rscript
# =============================================================================
# Seurat Single-Cell RNA-seq Analysis Pipeline
# =============================================================================
#
# Production-ready workflow for single-cell RNA-seq data analysis
#
# Author: Ji (gynecoloji)
# Date: 2025-01-03
#
# Usage:
#   Rscript run_seurat_pipeline.R --config config/default_config.yaml
#
# =============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
})

# =============================================================================
# COMMAND LINE ARGUMENTS
# =============================================================================

option_list <- list(
  make_option(
    c("-c", "--config"),
    type = "character",
    default = "config/default_config.yaml",
    help = "Path to configuration YAML file [default: %default]",
    metavar = "FILE"
  ),
  make_option(
    c("-i", "--input"),
    type = "character",
    default = NULL,
    help = "Path to input data (overrides config file)",
    metavar = "FILE"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = NULL,
    help = "Output directory (overrides config file)",
    metavar = "DIR"
  ),
  make_option(
    c("-n", "--cores"),
    type = "integer",
    default = NULL,
    help = "Number of cores to use (overrides config file)",
    metavar = "INT"
  ),
  make_option(
    c("--skip-qc"),
    action = "store_true",
    default = FALSE,
    help = "Skip QC filtering step"
  ),
  make_option(
    c("--skip-batch"),
    action = "store_true",
    default = FALSE,
    help = "Skip batch correction step"
  )
)

parser <- OptionParser(
  usage = "%prog [options]",
  option_list = option_list,
  description = "\nProduction-ready Seurat pipeline for single-cell RNA-seq analysis"
)

args <- parse_args(parser)

# =============================================================================
# SETUP
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("                 Seurat Single-Cell RNA-seq Pipeline                           \n")
cat("================================================================================\n")
cat("\n")

# Source all pipeline functions
source("R/utils.R")
source("R/qc_functions.R")
source("R/processing_functions.R")
source("R/plotting_functions.R")

# Load configuration
cat("Loading configuration...\n")
config <- load_config(args$config, verbose = TRUE)

# Override config with command line arguments
if (!is.null(args$input)) {
  config$paths$input_data <- args$input
}

if (!is.null(args$output)) {
  config$paths$output_dir <- args$output
}

if (!is.null(args$cores)) {
  config$processing$n_cores <- args$cores
  config$batch_correction$n_cores <- args$cores
}

# Create output directories
create_output_dir(config$paths$output_dir, verbose = TRUE)
create_output_dir(config$paths$qc_dir, verbose = TRUE)
create_output_dir(config$paths$processing_dir, verbose = TRUE)
create_output_dir(config$paths$plots_dir, verbose = TRUE)
create_output_dir(config$paths$final_objects_dir, verbose = TRUE)
create_output_dir(dirname(config$paths$log_file), verbose = TRUE)

# Create logger
logger <- create_logger(
  log_file = config$paths$log_file,
  verbose = config$general$verbose
)

logger("=== Pipeline started ===", "INFO")
logger(sprintf("Configuration file: %s", args$config), "INFO")
logger(sprintf("Input data: %s", config$paths$input_data), "INFO")
logger(sprintf("Output directory: %s", config$paths$output_dir), "INFO")

# Set random seed for reproducibility
set.seed(config$general$random_seed)

# Log system information
sys_info <- get_system_info()
logger(sprintf("R version: %s", sys_info$r_version), "INFO")
logger(sprintf("Platform: %s", sys_info$platform), "INFO")
logger(sprintf("User: %s", sys_info$user), "INFO")

# =============================================================================
# STEP 1: LOAD DATA
# =============================================================================

logger("=== STEP 1: Loading data ===", "INFO")

if (!file.exists(config$paths$input_data)) {
  logger(sprintf("Input file not found: %s", config$paths$input_data), "ERROR")
  stop("Input file not found")
}

# Load Seurat object
seurat_obj <- load_seurat_object(
  config$paths$input_data,
  verbose = config$general$verbose
)

logger(sprintf("Loaded %d cells and %d features", 
               ncol(seurat_obj), nrow(seurat_obj)), "INFO")

# =============================================================================
# STEP 2: QUALITY CONTROL
# =============================================================================

if (!args$skip_qc) {
  logger("=== STEP 2: Quality Control ===", "INFO")
  
  # Add QC metrics
  logger("Adding QC metrics...", "INFO")
  seurat_obj <- add_qc_metrics(
    seurat_obj,
    species = config$general$species,
    verbose = config$general$verbose
  )
  
  # Plot QC metrics
  if (config$qc$save_qc_plots) {
    logger("Generating QC plots...", "INFO")
    plot_qc_metrics(
      seurat_obj = seurat_obj,
      sample_name = config$general$project_name,
      output_dir = config$paths$qc_dir,
      plot_width = config$qc$plot_width,
      plot_height = config$qc$plot_height,
      qc_thresholds = list(
        min_umi = config$qc$min_umi,
        min_genes = config$qc$min_genes,
        max_mito = config$qc$max_mito,
        min_complexity = config$qc$min_complexity
      ),
      verbose = config$general$verbose
    )
  }
  
  # Filter cells and genes
  logger("Filtering cells and genes...", "INFO")
  seurat_filtered <- filter_cells_and_genes(
    seurat_obj = seurat_obj,
    sample_name = config$general$project_name,
    output_dir = config$paths$qc_dir,
    min_umi = config$qc$min_umi,
    min_genes = config$qc$min_genes,
    max_mito = config$qc$max_mito,
    min_complexity = config$qc$min_complexity,
    min_cells_per_gene = config$qc$min_cells_per_gene,
    save_stats = config$output$save_statistics,
    verbose = config$general$verbose
  )
  
  # Save QC object
  if (config$output$save_qc_object) {
    qc_file <- file.path(
      config$paths$final_objects_dir,
      paste0(config$general$project_name, "_qc_filtered.rds")
    )
    
    save_seurat_object(
      seurat_obj = seurat_filtered,
      file_path = qc_file,
      compress = config$output$compression,
      verbose = config$general$verbose
    )
  }
  
  seurat_obj <- seurat_filtered
  
} else {
  logger("Skipping QC filtering (--skip-qc flag used)", "WARNING")
}

# =============================================================================
# STEP 3: PROCESSING
# =============================================================================

logger("=== STEP 3: Processing ===", "INFO")

if (config$batch_correction$enabled && !args$skip_batch) {
  
  logger("Batch correction enabled", "INFO")
  
  # Process with batch correction
  seurat_processed <- process_batch_correction(
    seurat_obj = seurat_obj,
    batch_var = config$batch_correction$batch_var,
    config = list(
      method = config$batch_correction$method,
      dims = config$batch_correction$dims,
      resolution = config$batch_correction$resolution,
      vars_to_regress = config$batch_correction$vars_to_regress,
      harmony_theta = config$batch_correction$harmony$theta,
      harmony_lambda = config$batch_correction$harmony$lambda,
      harmony_max_iter = config$batch_correction$harmony$max_iter,
      k_anchor = config$batch_correction$integration$k_anchor,
      k_filter = config$batch_correction$integration$k_filter,
      k_score = config$batch_correction$integration$k_score,
      sct_n_genes = config$batch_correction$sct$n_genes,
      sct_vars_to_regress = config$batch_correction$sct$vars_to_regress,
      run_umap = config$batch_correction$run_umap,
      run_tsne = config$batch_correction$run_tsne,
      n_cores = config$batch_correction$n_cores,
      future_gb = config$batch_correction$future_gb
    ),
    species = config$general$species,
    verbose = config$general$verbose
  )
  
} else {
  
  if (args$skip_batch) {
    logger("Skipping batch correction (--skip-batch flag used)", "WARNING")
  } else {
    logger("Batch correction not enabled in config", "INFO")
  }
  
  # Standard processing without batch correction
  seurat_processed <- process_seurat_object(
    seurat_obj = seurat_obj,
    config = list(
      normalization_method = config$processing$normalization_method,
      scale_factor = config$processing$scale_factor,
      selection_method = config$processing$selection_method,
      n_features = config$processing$n_features,
      vars_to_regress = config$processing$vars_to_regress,
      scale_all_genes = config$processing$scale_all_genes,
      regress_cell_cycle = config$processing$regress_cell_cycle,
      n_pcs = config$processing$n_pcs,
      dims_use = config$processing$dims_use,
      resolution = config$processing$resolution,
      algorithm = config$processing$algorithm,
      run_umap = config$processing$run_umap,
      run_tsne = config$processing$run_tsne,
      umap_dims = config$processing$umap_dims,
      tsne_dims = config$processing$tsne_dims,
      run_jackstraw = config$processing$run_jackstraw,
      jackstraw_replicates = config$processing$jackstraw_replicates,
      n_cores = config$processing$n_cores,
      future_gb = config$processing$future_gb
    ),
    species = config$general$species,
    verbose = config$general$verbose
  )
}

# Save processed object
if (config$output$save_processed_object) {
  proc_file <- file.path(
    config$paths$final_objects_dir,
    paste0(config$general$project_name, "_processed.rds")
  )
  
  save_seurat_object(
    seurat_obj = seurat_processed,
    file_path = proc_file,
    compress = config$output$compression,
    verbose = config$general$verbose
  )
}

# =============================================================================
# STEP 4: VISUALIZATION
# =============================================================================

logger("=== STEP 4: Generating Plots ===", "INFO")

if (config$plotting$save_plots) {
  
  # Elbow plot
  logger("Creating elbow plot...", "INFO")
  plot_elbow(
    seurat_obj = seurat_processed,
    ndims = config$processing$n_pcs,
    output_dir = config$paths$plots_dir,
    sample_name = config$general$project_name,
    save_plot = TRUE,
    verbose = config$general$verbose
  )
  
  # Cluster QC plots
  logger("Creating cluster QC plots...", "INFO")
  plot_cluster_qc(
    seurat_obj = seurat_processed,
    output_dir = config$paths$plots_dir,
    sample_name = config$general$project_name,
    save_plots = TRUE,
    verbose = config$general$verbose
  )
  
  # UMAP/tSNE plots
  if (config$processing$run_umap || config$batch_correction$run_umap) {
    logger("Creating UMAP visualization...", "INFO")
    
    if (config$batch_correction$enabled && !args$skip_batch) {
      # Batch correction plots
      plot_batch_correction_results(
        seurat_obj = seurat_processed,
        batch_var = config$batch_correction$batch_var,
        output_dir = config$paths$plots_dir,
        sample_name = config$general$project_name,
        reduction = "umap",
        point_size = config$plotting$point_size,
        label_clusters = config$plotting$label_clusters,
        save_plots = TRUE,
        verbose = config$general$verbose
      )
    } else {
      # Standard UMAP plot
      p_umap <- DimPlot(
        seurat_processed,
        reduction = "umap",
        group.by = "seurat_clusters",
        label = config$plotting$label_clusters,
        pt.size = config$plotting$point_size
      )
      
      ggsave(
        filename = file.path(config$paths$plots_dir, 
                           paste0(config$general$project_name, "_umap.pdf")),
        plot = p_umap,
        width = 8,
        height = 7
      )
    }
  }
  
  # Feature plots for marker genes
  if (length(config$plotting$marker_genes) > 0) {
    logger("Creating marker gene expression plots...", "INFO")
    
    plot_gene_expression(
      seurat_obj = seurat_processed,
      features = config$plotting$marker_genes,
      reduction = "umap",
      output_dir = config$paths$plots_dir,
      sample_name = paste0(config$general$project_name, "_markers"),
      ncol = 2,
      point_size = config$plotting$point_size,
      save_plots = TRUE,
      verbose = config$general$verbose
    )
  }
}

# =============================================================================
# STEP 5: FINAL OUTPUT
# =============================================================================

logger("=== STEP 5: Saving Final Results ===", "INFO")

# Save final processed object
final_file <- file.path(
  config$paths$final_objects_dir,
  paste0(config$general$project_name, "_final.rds")
)

save_seurat_object(
  seurat_obj = seurat_processed,
  file_path = final_file,
  compress = config$output$compression,
  verbose = config$general$verbose
)

# Save processing summary
if (config$output$save_statistics) {
  summary_df <- get_processing_summary(seurat_processed)
  
  summary_file <- file.path(
    config$paths$final_objects_dir,
    paste0(config$general$project_name, "_summary.csv")
  )
  
  write.csv(summary_df, summary_file, row.names = FALSE)
  logger(sprintf("Processing summary saved to: %s", summary_file), "INFO")
}

# Save session info
if (config$output$save_session_info) {
  session_file <- file.path(
    config$paths$final_objects_dir,
    paste0(config$general$project_name, "_sessionInfo.txt")
  )
  
  print_session_info(session_file)
  logger(sprintf("Session info saved to: %s", session_file), "INFO")
}

# Save configuration used for this run
config_copy <- file.path(
  config$paths$final_objects_dir,
  paste0(config$general$project_name, "_config.yaml")
)

save_config(config, config_copy, verbose = config$general$verbose)

# =============================================================================
# COMPLETE
# =============================================================================

logger("=== Pipeline completed successfully ===", "SUCCESS")

cat("\n")
cat("================================================================================\n")
cat("                         Pipeline Completed Successfully!                      \n")
cat("================================================================================\n")
cat("\n")
cat(sprintf("Results saved to: %s\n", config$paths$output_dir))
cat(sprintf("Final object: %s\n", final_file))
cat(sprintf("Plots directory: %s\n", config$paths$plots_dir))
cat(sprintf("Log file: %s\n", config$paths$log_file))
cat("\n")
cat("Summary:\n")
cat(sprintf("  - Total cells: %d\n", ncol(seurat_processed)))
cat(sprintf("  - Total features: %d\n", nrow(seurat_processed)))
if ("seurat_clusters" %in% colnames(seurat_processed@meta.data)) {
  cat(sprintf("  - Number of clusters: %d\n", 
              length(levels(seurat_processed$seurat_clusters))))
}
cat("\n")
