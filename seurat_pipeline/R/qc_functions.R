#' Quality Control Functions for Single-Cell RNA-seq Analysis
#'
#' This module provides production-ready functions for quality control
#' of single-cell RNA-seq data using Seurat.
#'
#' @author Ji (gynecoloji)
#' @date 2025-01-03

# Required packages
#' @import Seurat
#' @import ggplot2
#' @import dplyr
#' @import Matrix

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Validate Seurat Object
#'
#' @param seurat_obj Object to validate
#' @param require_counts Logical, whether to require counts data
#' @return NULL (throws error if invalid)
#' @keywords internal
validate_seurat_object <- function(seurat_obj, require_counts = TRUE) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object. Got: ", class(seurat_obj)[1])
  }
  
  if (ncol(seurat_obj) == 0) {
    stop("Seurat object contains zero cells")
  }
  
  if (nrow(seurat_obj) == 0) {
    stop("Seurat object contains zero features")
  }
  
  if (require_counts) {
    if (!"RNA" %in% names(seurat_obj@assays)) {
      stop("Seurat object must contain 'RNA' assay")
    }
  }
  
  invisible(NULL)
}

#' Create Output Directory
#'
#' @param path Directory path to create
#' @param verbose Logical, print messages
#' @return Character, the created path
#' @keywords internal
create_output_dir <- function(path, verbose = TRUE) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (verbose) {
      message("Created directory: ", path)
    }
  }
  return(normalizePath(path))
}

#' Log Message
#'
#' @param msg Message to log
#' @param level Log level: "INFO", "WARNING", "ERROR"
#' @param verbose Logical, whether to print
#' @keywords internal
log_message <- function(msg, level = "INFO", verbose = TRUE) {
  if (!verbose) return(invisible(NULL))
  
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  formatted_msg <- sprintf("[%s] %s: %s", timestamp, level, msg)
  
  if (level == "ERROR") {
    message(formatted_msg)
  } else if (level == "WARNING") {
    message(formatted_msg)
  } else {
    message(formatted_msg)
  }
  
  invisible(NULL)
}

# =============================================================================
# QC METRICS CALCULATION
# =============================================================================

#' Add QC Metrics to Seurat Object
#'
#' Calculates and adds quality control metrics to a Seurat object's metadata.
#' Metrics include mitochondrial ratio, ribosomal ratio, and gene complexity.
#'
#' @param seurat_obj A Seurat object
#' @param mito_pattern Regex pattern for mitochondrial genes (default: "^MT-")
#' @param ribo_pattern Regex pattern for ribosomal genes (default: "^RP[SL]")
#' @param species Character, species for gene patterns ("human" or "mouse")
#' @param verbose Logical, print progress messages
#'
#' @return Seurat object with additional metadata columns:
#'   \itemize{
#'     \item mitoRatio: Percentage of mitochondrial reads
#'     \item riboRatio: Percentage of ribosomal reads  
#'     \item log10GenesPerUMI: Gene complexity metric
#'   }
#'
#' @examples
#' \dontrun{
#' seurat_obj <- add_qc_metrics(seurat_obj, species = "human")
#' }
#'
#' @export
add_qc_metrics <- function(seurat_obj, 
                           mito_pattern = NULL,
                           ribo_pattern = NULL,
                           species = "human",
                           verbose = TRUE) {
  
  # Validate input
  validate_seurat_object(seurat_obj, require_counts = TRUE)
  
  # Set species-specific patterns
  if (is.null(mito_pattern)) {
    mito_pattern <- if (species == "human") "^MT-" else "^mt-"
  }
  
  if (is.null(ribo_pattern)) {
    ribo_pattern <- "^RP[SL]"
  }
  
  log_message("Adding QC metrics to Seurat object", verbose = verbose)
  
  # Calculate mitochondrial percentage
  tryCatch({
    seurat_obj$mitoRatio <- PercentageFeatureSet(
      object = seurat_obj, 
      pattern = mito_pattern
    )
    log_message(
      sprintf("Mitochondrial genes: %.2f%% average", 
              mean(seurat_obj$mitoRatio)),
      verbose = verbose
    )
  }, error = function(e) {
    warning("Could not calculate mitochondrial ratio: ", e$message)
    seurat_obj$mitoRatio <- 0
  })
  
  # Calculate ribosomal percentage
  tryCatch({
    seurat_obj$riboRatio <- PercentageFeatureSet(
      object = seurat_obj,
      pattern = ribo_pattern
    )
    log_message(
      sprintf("Ribosomal genes: %.2f%% average", 
              mean(seurat_obj$riboRatio)),
      verbose = verbose
    )
  }, error = function(e) {
    warning("Could not calculate ribosomal ratio: ", e$message)
    seurat_obj$riboRatio <- 0
  })
  
  # Calculate gene complexity (genes per UMI)
  seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / 
                                  log10(seurat_obj$nCount_RNA)
  
  # Handle infinite/NaN values
  seurat_obj$log10GenesPerUMI[is.infinite(seurat_obj$log10GenesPerUMI)] <- NA
  seurat_obj$log10GenesPerUMI[is.nan(seurat_obj$log10GenesPerUMI)] <- NA
  
  log_message(
    sprintf("Gene complexity: %.3f average", 
            mean(seurat_obj$log10GenesPerUMI, na.rm = TRUE)),
    verbose = verbose
  )
  
  log_message("QC metrics added successfully", verbose = verbose)
  
  return(seurat_obj)
}

# =============================================================================
# QC VISUALIZATION
# =============================================================================

#' Plot Quality Control Metrics
#'
#' Creates comprehensive QC plots for single-cell data quality assessment.
#' Generates 5 PDF plots: UMI distribution, gene count distribution,
#' mitochondrial ratio, correlation plots, and gene complexity.
#'
#' @param seurat_obj Seurat object with QC metrics (run add_qc_metrics first)
#' @param sample_name Name of sample for plot titles and file naming
#' @param output_dir Directory to save plots
#' @param plot_width Plot width in inches (default: 8)
#' @param plot_height Plot height in inches (default: 7)
#' @param qc_thresholds Named list of QC thresholds for visualization:
#'   \itemize{
#'     \item min_umi: Minimum UMI count (default: 500)
#'     \item min_genes: Minimum gene count (default: 300)
#'     \item max_mito: Maximum mitochondrial % (default: 25)
#'     \item min_complexity: Minimum log10GenesPerUMI (default: 0.8)
#'   }
#' @param verbose Logical, print progress messages
#'
#' @return List of ggplot objects
#'
#' @examples
#' \dontrun{
#' plots <- plot_qc_metrics(
#'   seurat_obj, 
#'   sample_name = "sample1",
#'   output_dir = "results/qc"
#' )
#' }
#'
#' @export
plot_qc_metrics <- function(seurat_obj,
                            sample_name,
                            output_dir,
                            plot_width = 8,
                            plot_height = 7,
                            qc_thresholds = list(
                              min_umi = 500,
                              min_genes = 300,
                              max_mito = 25,
                              min_complexity = 0.8
                            ),
                            verbose = TRUE) {
  
  # Validate inputs
  validate_seurat_object(seurat_obj)
  
  if (missing(sample_name) || is.null(sample_name) || sample_name == "") {
    stop("sample_name is required and cannot be empty")
  }
  
  if (missing(output_dir) || is.null(output_dir)) {
    stop("output_dir is required")
  }
  
  # Check required metadata columns
  required_cols <- c("nCount_RNA", "nFeature_RNA", "mitoRatio", "log10GenesPerUMI")
  missing_cols <- setdiff(required_cols, colnames(seurat_obj@meta.data))
  
  if (length(missing_cols) > 0) {
    stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "),
         "\nRun add_qc_metrics() first")
  }
  
  # Create output directory
  sample_dir <- file.path(output_dir, sample_name)
  create_output_dir(sample_dir, verbose = verbose)
  
  log_message(
    sprintf("Generating QC plots for %s (%d cells)", 
            sample_name, ncol(seurat_obj)),
    verbose = verbose
  )
  
  # Prepare data
  df <- seurat_obj@meta.data %>%
    dplyr::rename(nUMI = nCount_RNA, nGene = nFeature_RNA) %>%
    dplyr::mutate(sample = sample_name)
  
  # Store plots
  plots <- list()
  
  # Plot 1: UMI count distribution
  log_message("Creating UMI distribution plot", verbose = verbose)
  plots$umi_density <- ggplot(df, aes(x = nUMI, fill = sample, color = sample)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10(labels = scales::comma) + 
    theme_classic() +
    theme(legend.position = "top") +
    labs(
      title = paste("UMI Count Distribution -", sample_name),
      x = "UMI Count (log10 scale)",
      y = "Cell Density"
    ) +
    geom_vline(xintercept = qc_thresholds$min_umi, 
               linetype = "dashed", color = "red", size = 0.5) +
    annotate("text", x = qc_thresholds$min_umi, y = Inf, 
             label = paste("Min UMI:", qc_thresholds$min_umi),
             hjust = -0.1, vjust = 2, size = 3)
  
  ggsave(
    filename = file.path(sample_dir, "QC_Cell_Density_UMI.pdf"),
    plot = plots$umi_density,
    width = plot_width,
    height = plot_height
  )
  
  # Plot 2: Gene count distribution
  log_message("Creating gene count distribution plot", verbose = verbose)
  plots$gene_density <- ggplot(df, aes(x = nGene, fill = sample, color = sample)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10(labels = scales::comma) + 
    theme_classic() +
    theme(legend.position = "top") +
    labs(
      title = paste("Gene Count Distribution -", sample_name),
      x = "Gene Count (log10 scale)",
      y = "Cell Density"
    ) +
    geom_vline(xintercept = qc_thresholds$min_genes, 
               linetype = "dashed", color = "red", size = 0.5) +
    annotate("text", x = qc_thresholds$min_genes, y = Inf,
             label = paste("Min Genes:", qc_thresholds$min_genes),
             hjust = -0.1, vjust = 2, size = 3)
  
  ggsave(
    filename = file.path(sample_dir, "QC_Cell_Density_genes.pdf"),
    plot = plots$gene_density,
    width = plot_width,
    height = plot_height
  )
  
  # Plot 3: Mitochondrial ratio distribution
  log_message("Creating mitochondrial ratio plot", verbose = verbose)
  plots$mito_density <- ggplot(df, aes(x = mitoRatio, fill = sample, color = sample)) + 
    geom_density(alpha = 0.2) +
    theme_classic() +
    theme(legend.position = "top") +
    labs(
      title = paste("Mitochondrial Ratio Distribution -", sample_name),
      x = "Mitochondrial Ratio (%)",
      y = "Cell Density"
    ) +
    geom_vline(xintercept = qc_thresholds$max_mito, 
               linetype = "dashed", color = "red", size = 0.5) +
    annotate("text", x = qc_thresholds$max_mito, y = Inf,
             label = paste("Max Mito:", qc_thresholds$max_mito, "%"),
             hjust = -0.1, vjust = 2, size = 3)
  
  ggsave(
    filename = file.path(sample_dir, "QC_Cell_Density_mitoRatio.pdf"),
    plot = plots$mito_density,
    width = plot_width,
    height = plot_height
  )
  
  # Plot 4: Scatter plot - UMI vs Genes colored by mito ratio
  log_message("Creating correlation plot", verbose = verbose)
  plots$correlation <- ggplot(df, aes(x = nUMI, y = nGene, color = mitoRatio)) + 
    geom_point(aes(size = log10GenesPerUMI), alpha = 0.5) + 
    scale_colour_gradient(low = "gray90", high = "black") +
    geom_smooth(method = "lm", color = "blue", se = TRUE) +
    scale_x_log10(labels = scales::comma) + 
    scale_y_log10(labels = scales::comma) + 
    theme_classic() +
    theme(legend.position = "right") +
    labs(
      title = paste("UMI vs Gene Count -", sample_name),
      x = "UMI Count (log10 scale)",
      y = "Gene Count (log10 scale)",
      color = "Mito %",
      size = "Complexity"
    ) +
    geom_vline(xintercept = qc_thresholds$min_umi, 
               linetype = "dashed", color = "red", size = 0.3) +
    geom_hline(yintercept = qc_thresholds$min_genes, 
               linetype = "dashed", color = "red", size = 0.3)
  
  ggsave(
    filename = file.path(sample_dir, "QC_Cell_correlations.pdf"),
    plot = plots$correlation,
    width = plot_width,
    height = plot_height
  )
  
  # Plot 5: Gene complexity distribution
  log_message("Creating gene complexity plot", verbose = verbose)
  plots$complexity <- ggplot(df, aes(x = log10GenesPerUMI, fill = sample, color = sample)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    theme(legend.position = "top") +
    labs(
      title = paste("Gene Complexity Distribution -", sample_name),
      x = "log10(Genes per UMI)",
      y = "Cell Density"
    ) +
    geom_vline(xintercept = qc_thresholds$min_complexity, 
               linetype = "dashed", color = "red", size = 0.5) +
    annotate("text", x = qc_thresholds$min_complexity, y = Inf,
             label = paste("Min Complexity:", qc_thresholds$min_complexity),
             hjust = -0.1, vjust = 2, size = 3)
  
  ggsave(
    filename = file.path(sample_dir, "QC_Cell_complexity.pdf"),
    plot = plots$complexity,
    width = plot_width,
    height = plot_height
  )
  
  log_message(
    sprintf("QC plots saved to: %s", sample_dir),
    verbose = verbose
  )
  
  return(invisible(plots))
}

# =============================================================================
# QC FILTERING
# =============================================================================

#' Filter Cells Based on QC Metrics
#'
#' Filters cells based on quality control thresholds. Provides detailed
#' statistics about filtering results.
#'
#' @param seurat_obj Seurat object with QC metrics
#' @param sample_name Sample name for logging
#' @param output_dir Directory to save filtering statistics and plots
#' @param min_umi Minimum UMI count threshold (default: 500)
#' @param min_genes Minimum gene count threshold (default: 300)
#' @param max_mito Maximum mitochondrial ratio threshold (default: 25)
#' @param min_complexity Minimum log10GenesPerUMI threshold (default: 0.8)
#' @param min_cells_per_gene Minimum cells a gene must be expressed in (default: 10)
#' @param save_stats Logical, save filtering statistics to file
#' @param verbose Logical, print progress messages
#'
#' @return Filtered Seurat object
#'
#' @examples
#' \dontrun{
#' seurat_filtered <- filter_cells_and_genes(
#'   seurat_obj,
#'   sample_name = "sample1",
#'   output_dir = "results/qc",
#'   min_umi = 1000,
#'   max_mito = 20
#' )
#' }
#'
#' @export
filter_cells_and_genes <- function(seurat_obj,
                                   sample_name,
                                   output_dir,
                                   min_umi = 500,
                                   min_genes = 300,
                                   max_mito = 25,
                                   min_complexity = 0.8,
                                   min_cells_per_gene = 10,
                                   save_stats = TRUE,
                                   verbose = TRUE) {
  
  # Validate inputs
  validate_seurat_object(seurat_obj)
  
  if (missing(sample_name) || missing(output_dir)) {
    stop("sample_name and output_dir are required")
  }
  
  # Check required columns
  required_cols <- c("nCount_RNA", "nFeature_RNA", "mitoRatio", "log10GenesPerUMI")
  missing_cols <- setdiff(required_cols, colnames(seurat_obj@meta.data))
  
  if (length(missing_cols) > 0) {
    stop("Missing required metadata: ", paste(missing_cols, collapse = ", "))
  }
  
  # Record initial counts
  n_cells_before <- ncol(seurat_obj)
  n_genes_before <- nrow(seurat_obj)
  
  log_message(
    sprintf("Starting filtering for %s: %d cells, %d genes",
            sample_name, n_cells_before, n_genes_before),
    verbose = verbose
  )
  
  # ===================
  # CELL FILTERING
  # ===================
  
  log_message("Filtering cells...", verbose = verbose)
  
  # Apply filters
  seurat_filtered <- subset(
    x = seurat_obj,
    subset = (nCount_RNA >= min_umi) &
             (nFeature_RNA >= min_genes) &
             (log10GenesPerUMI > min_complexity) &
             (mitoRatio < max_mito)
  )
  
  n_cells_after <- ncol(seurat_filtered)
  cells_removed <- n_cells_before - n_cells_after
  cells_kept_pct <- round((n_cells_after / n_cells_before) * 100, 2)
  
  log_message(
    sprintf("Cells: %d -> %d (removed %d, kept %.2f%%)",
            n_cells_before, n_cells_after, cells_removed, cells_kept_pct),
    verbose = verbose
  )
  
  if (n_cells_after == 0) {
    stop("All cells were filtered out. Consider relaxing QC thresholds.")
  }
  
  # ===================
  # GENE FILTERING
  # ===================
  
  log_message("Filtering genes...", verbose = verbose)
  
  counts <- GetAssayData(object = seurat_filtered, slot = "counts")
  nonzero <- counts > 0
  
  # Calculate genes expressed per cell
  genes_per_cell <- Matrix::colSums(nonzero)
  
  # Plot gene expression distribution
  sample_dir <- file.path(output_dir, sample_name)
  create_output_dir(sample_dir, verbose = FALSE)
  
  gene_df <- data.frame(
    Gene_in_nCells = Matrix::rowSums(nonzero),
    sample = sample_name
  )
  
  p_gene <- ggplot(gene_df, aes(x = Gene_in_nCells, fill = sample, color = sample)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10(labels = scales::comma) +
    labs(
      title = paste("Gene Expression Distribution -", sample_name),
      x = "Number of Cells Expressing Gene (log10)",
      y = "Density"
    ) +
    geom_vline(xintercept = min_cells_per_gene, 
               linetype = "dashed", color = "red") +
    annotate("text", x = min_cells_per_gene, y = Inf,
             label = paste("Min cells:", min_cells_per_gene),
             hjust = -0.1, vjust = 2, size = 3)
  
  ggsave(
    filename = file.path(sample_dir, "QC_genes_density.pdf"),
    plot = p_gene,
    width = 8,
    height = 7
  )
  
  # Filter genes
  keep_genes <- Matrix::rowSums(nonzero) >= min_cells_per_gene
  filtered_counts <- counts[keep_genes, ]
  
  n_genes_after <- nrow(filtered_counts)
  genes_removed <- n_genes_before - n_genes_after
  genes_kept_pct <- round((n_genes_after / n_genes_before) * 100, 2)
  
  log_message(
    sprintf("Genes: %d -> %d (removed %d, kept %.2f%%)",
            n_genes_before, n_genes_after, genes_removed, genes_kept_pct),
    verbose = verbose
  )
  
  # Create new Seurat object with filtered data
  seurat_filtered <- CreateSeuratObject(
    counts = filtered_counts,
    meta.data = seurat_filtered@meta.data
  )
  
  # ===================
  # SAVE STATISTICS
  # ===================
  
  if (save_stats) {
    stats <- data.frame(
      sample = sample_name,
      cells_before = n_cells_before,
      cells_after = n_cells_after,
      cells_removed = cells_removed,
      cells_kept_pct = cells_kept_pct,
      genes_before = n_genes_before,
      genes_after = n_genes_after,
      genes_removed = genes_removed,
      genes_kept_pct = genes_kept_pct,
      min_umi = min_umi,
      min_genes = min_genes,
      max_mito = max_mito,
      min_complexity = min_complexity,
      min_cells_per_gene = min_cells_per_gene,
      timestamp = Sys.time()
    )
    
    stats_file <- file.path(sample_dir, "filtering_statistics.csv")
    write.csv(stats, stats_file, row.names = FALSE)
    log_message(sprintf("Statistics saved to: %s", stats_file), verbose = verbose)
  }
  
  log_message("Filtering completed successfully", verbose = verbose)
  
  return(seurat_filtered)
}
