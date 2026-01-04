#' Plotting Functions for Seurat Analysis Results
#'
#' Production-ready visualization functions for single-cell RNA-seq analysis
#'
#' @author Ji (gynecoloji)
#' @date 2025-01-03

# Required packages
#' @import ggplot2
#' @import patchwork
#' @import RColorBrewer
#' @import Seurat

# =============================================================================
# BATCH CORRECTION VISUALIZATION
# =============================================================================

#' Plot Batch Correction Results
#'
#' Generate comprehensive diagnostic plots to evaluate batch correction effectiveness.
#' Creates side-by-side comparisons and diagnostic panels.
#'
#' @param seurat_obj Processed Seurat object with batch correction
#' @param batch_var Batch variable column name
#' @param output_dir Directory to save plots (optional)
#' @param sample_name Sample name for file naming (optional)
#' @param reduction Reduction to plot ("umap" or "tsne", default: "umap")
#' @param plot_width Plot width in inches (default: 12)
#' @param plot_height Plot height in inches (default: 5)
#' @param point_size Point size for scatter plots (default: 0.5)
#' @param label_clusters Logical, add cluster labels (default: TRUE)
#' @param save_plots Logical, save plots to files (default: TRUE)
#' @param verbose Logical, print messages
#'
#' @return List of ggplot objects
#'
#' @examples
#' \dontrun{
#' plots <- plot_batch_correction_results(
#'   seurat_obj,
#'   batch_var = "orig.ident",
#'   output_dir = "results/plots"
#' )
#' }
#'
#' @export
plot_batch_correction_results <- function(seurat_obj,
                                          batch_var,
                                          output_dir = NULL,
                                          sample_name = "batch_correction",
                                          reduction = "umap",
                                          plot_width = 12,
                                          plot_height = 5,
                                          point_size = 0.5,
                                          label_clusters = TRUE,
                                          save_plots = TRUE,
                                          verbose = TRUE) {
  
  # Validate inputs
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object")
  }
  
  if (!batch_var %in% colnames(seurat_obj@meta.data)) {
    stop(sprintf("Batch variable '%s' not found in metadata", batch_var))
  }
  
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop(sprintf("Reduction '%s' not found. Available: %s",
                 reduction, paste(names(seurat_obj@reductions), collapse = ", ")))
  }
  
  log_message("Generating batch correction plots", verbose = verbose)
  
  plots <- list()
  
  # ========================================
  # Plot 1: Side-by-side comparison
  # ========================================
  
  log_message("Creating batch/cluster comparison plot", verbose = verbose)
  
  # Plot colored by batch
  p_batch <- DimPlot(
    object = seurat_obj,
    reduction = reduction,
    group.by = batch_var,
    pt.size = point_size,
    raster = FALSE
  ) +
    ggtitle(paste("Colored by", batch_var)) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )
  
  # Plot colored by clusters
  p_clusters <- DimPlot(
    object = seurat_obj,
    reduction = reduction,
    group.by = "seurat_clusters",
    pt.size = point_size,
    label = label_clusters,
    label.size = 3,
    raster = FALSE
  ) +
    ggtitle("Colored by Clusters") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )
  
  # Combine plots
  plots$comparison <- p_batch + p_clusters +
    plot_annotation(
      title = paste(toupper(reduction), "Visualization"),
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    )
  
  # Save if requested
  if (save_plots && !is.null(output_dir)) {
    create_output_dir(output_dir, verbose = FALSE)
    
    filename <- file.path(
      output_dir,
      paste0(sample_name, "_", reduction, "_comparison.pdf")
    )
    
    ggsave(
      filename = filename,
      plot = plots$comparison,
      width = plot_width,
      height = plot_height
    )
    
    log_message(sprintf("Saved: %s", filename), verbose = verbose)
  }
  
  # ========================================
  # Plot 2: Batch distribution per cluster
  # ========================================
  
  log_message("Creating batch distribution plot", verbose = verbose)
  
  # Calculate proportions
  batch_cluster_df <- seurat_obj@meta.data %>%
    dplyr::group_by(seurat_clusters, .data[[batch_var]]) %>%
    dplyr::summarise(n = n(), .groups = "drop") %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::mutate(proportion = n / sum(n))
  
  plots$batch_distribution <- ggplot(
    batch_cluster_df,
    aes(x = seurat_clusters, y = proportion, fill = .data[[batch_var]])
  ) +
    geom_bar(stat = "identity", position = "fill") +
    scale_y_continuous(labels = scales::percent) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    ) +
    labs(
      title = "Batch Distribution Across Clusters",
      x = "Cluster",
      y = "Proportion",
      fill = batch_var
    )
  
  # Save if requested
  if (save_plots && !is.null(output_dir)) {
    filename <- file.path(
      output_dir,
      paste0(sample_name, "_batch_distribution.pdf")
    )
    
    ggsave(
      filename = filename,
      plot = plots$batch_distribution,
      width = 8,
      height = 6
    )
    
    log_message(sprintf("Saved: %s", filename), verbose = verbose)
  }
  
  # ========================================
  # Plot 3: Cells per batch and cluster
  # ========================================
  
  log_message("Creating cells per cluster plot", verbose = verbose)
  
  plots$cells_per_cluster <- ggplot(
    batch_cluster_df,
    aes(x = seurat_clusters, y = n, fill = .data[[batch_var]])
  ) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    ) +
    labs(
      title = "Cell Counts per Cluster and Batch",
      x = "Cluster",
      y = "Number of Cells",
      fill = batch_var
    )
  
  # Save if requested
  if (save_plots && !is.null(output_dir)) {
    filename <- file.path(
      output_dir,
      paste0(sample_name, "_cells_per_cluster.pdf")
    )
    
    ggsave(
      filename = filename,
      plot = plots$cells_per_cluster,
      width = 8,
      height = 6
    )
    
    log_message(sprintf("Saved: %s", filename), verbose = verbose)
  }
  
  log_message("Batch correction plots completed", verbose = verbose)
  
  return(invisible(plots))
}

# =============================================================================
# FEATURE PLOTS
# =============================================================================

#' Plot Gene Expression on Dimensionality Reduction
#'
#' Visualize expression of genes across cells in reduced dimensions
#'
#' @param seurat_obj Seurat object
#' @param features Character vector of gene names to plot
#' @param reduction Reduction to use (default: "umap")
#' @param output_dir Directory to save plots (optional)
#' @param sample_name Sample name for file naming
#' @param ncol Number of columns for multi-gene plots (default: 2)
#' @param point_size Point size (default: 0.5)
#' @param save_plots Logical, save plots to files
#' @param verbose Logical, print messages
#'
#' @return List of ggplot objects
#' @export
plot_gene_expression <- function(seurat_obj,
                                 features,
                                 reduction = "umap",
                                 output_dir = NULL,
                                 sample_name = "gene_expression",
                                 ncol = 2,
                                 point_size = 0.5,
                                 save_plots = TRUE,
                                 verbose = TRUE) {
  
  # Validate inputs
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object")
  }
  
  # Check which features are available
  available_features <- features[features %in% rownames(seurat_obj)]
  missing_features <- setdiff(features, available_features)
  
  if (length(missing_features) > 0) {
    warning(sprintf(
      "Features not found: %s",
      paste(missing_features, collapse = ", ")
    ))
  }
  
  if (length(available_features) == 0) {
    stop("None of the requested features were found in the object")
  }
  
  log_message(
    sprintf("Plotting %d features", length(available_features)),
    verbose = verbose
  )
  
  # Create feature plot
  p <- FeaturePlot(
    object = seurat_obj,
    features = available_features,
    reduction = reduction,
    pt.size = point_size,
    ncol = ncol,
    raster = FALSE
  )
  
  # Save if requested
  if (save_plots && !is.null(output_dir)) {
    create_output_dir(output_dir, verbose = FALSE)
    
    filename <- file.path(
      output_dir,
      paste0(sample_name, "_feature_plot.pdf")
    )
    
    # Calculate height based on number of rows
    n_rows <- ceiling(length(available_features) / ncol)
    plot_height <- max(5, n_rows * 3)
    
    ggsave(
      filename = filename,
      plot = p,
      width = ncol * 5,
      height = plot_height
    )
    
    log_message(sprintf("Saved: %s", filename), verbose = verbose)
  }
  
  return(invisible(p))
}

# =============================================================================
# ELBOW PLOT
# =============================================================================

#' Plot PCA Elbow Plot
#'
#' Visualize variance explained by principal components
#'
#' @param seurat_obj Seurat object with PCA computed
#' @param ndims Number of dimensions to show (default: 50)
#' @param output_dir Directory to save plot (optional)
#' @param sample_name Sample name for file naming
#' @param save_plot Logical, save plot to file
#' @param verbose Logical, print messages
#'
#' @return ggplot object
#' @export
plot_elbow <- function(seurat_obj,
                       ndims = 50,
                       output_dir = NULL,
                       sample_name = "elbow",
                       save_plot = TRUE,
                       verbose = TRUE) {
  
  # Validate input
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object")
  }
  
  if (!"pca" %in% names(seurat_obj@reductions)) {
    stop("PCA not found. Run RunPCA first.")
  }
  
  # Get standard deviations
  pca_stdev <- seurat_obj@reductions$pca@stdev
  
  if (ndims > length(pca_stdev)) {
    ndims <- length(pca_stdev)
    warning(sprintf("Only %d PCs available, plotting all", ndims))
  }
  
  # Create data frame
  elbow_df <- data.frame(
    PC = 1:ndims,
    stdev = pca_stdev[1:ndims]
  )
  
  # Create plot
  p <- ggplot(elbow_df, aes(x = PC, y = stdev)) +
    geom_point(size = 2, color = "steelblue") +
    geom_line(color = "steelblue", alpha = 0.5) +
    theme_classic() +
    labs(
      title = "PCA Elbow Plot",
      x = "Principal Component",
      y = "Standard Deviation"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  # Save if requested
  if (save_plot && !is.null(output_dir)) {
    create_output_dir(output_dir, verbose = FALSE)
    
    filename <- file.path(
      output_dir,
      paste0(sample_name, "_elbow_plot.pdf")
    )
    
    ggsave(
      filename = filename,
      plot = p,
      width = 8,
      height = 6
    )
    
    log_message(sprintf("Saved: %s", filename), verbose = verbose)
  }
  
  return(invisible(p))
}

# =============================================================================
# CLUSTER QC PLOTS
# =============================================================================

#' Plot Cluster QC Metrics
#'
#' Visualize QC metrics (UMI, genes, mito%) across clusters
#'
#' @param seurat_obj Seurat object
#' @param output_dir Directory to save plots (optional)
#' @param sample_name Sample name for file naming
#' @param save_plots Logical, save plots to files
#' @param verbose Logical, print messages
#'
#' @return List of ggplot objects
#' @export
plot_cluster_qc <- function(seurat_obj,
                            output_dir = NULL,
                            sample_name = "cluster_qc",
                            save_plots = TRUE,
                            verbose = TRUE) {
  
  # Validate input
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object")
  }
  
  if (!"seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    stop("Clusters not found. Run FindClusters first.")
  }
  
  log_message("Creating cluster QC plots", verbose = verbose)
  
  plots <- list()
  
  # Prepare data
  qc_df <- seurat_obj@meta.data %>%
    dplyr::select(
      seurat_clusters,
      nCount_RNA,
      nFeature_RNA,
      dplyr::any_of(c("percent.mt", "mitoRatio"))
    )
  
  # Rename mitoRatio to percent.mt if needed
  if ("mitoRatio" %in% colnames(qc_df) && !"percent.mt" %in% colnames(qc_df)) {
    qc_df$percent.mt <- qc_df$mitoRatio
  }
  
  # Plot 1: UMI counts per cluster
  plots$umi <- ggplot(qc_df, aes(x = seurat_clusters, y = nCount_RNA, fill = seurat_clusters)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    scale_y_log10(labels = scales::comma) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      title = "UMI Counts per Cluster",
      x = "Cluster",
      y = "UMI Count (log10)"
    )
  
  # Plot 2: Gene counts per cluster
  plots$genes <- ggplot(qc_df, aes(x = seurat_clusters, y = nFeature_RNA, fill = seurat_clusters)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    scale_y_log10(labels = scales::comma) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      title = "Gene Counts per Cluster",
      x = "Cluster",
      y = "Gene Count (log10)"
    )
  
  # Plot 3: Mitochondrial percentage per cluster (if available)
  if ("percent.mt" %in% colnames(qc_df)) {
    plots$mito <- ggplot(qc_df, aes(x = seurat_clusters, y = percent.mt, fill = seurat_clusters)) +
      geom_violin(alpha = 0.7) +
      geom_boxplot(width = 0.1, outlier.shape = NA) +
      theme_classic() +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      labs(
        title = "Mitochondrial % per Cluster",
        x = "Cluster",
        y = "Mitochondrial %"
      )
  }
  
  # Combine plots
  if ("percent.mt" %in% colnames(qc_df)) {
    plots$combined <- plots$umi + plots$genes + plots$mito +
      plot_layout(ncol = 3)
  } else {
    plots$combined <- plots$umi + plots$genes +
      plot_layout(ncol = 2)
  }
  
  # Save if requested
  if (save_plots && !is.null(output_dir)) {
    create_output_dir(output_dir, verbose = FALSE)
    
    filename <- file.path(
      output_dir,
      paste0(sample_name, "_cluster_qc.pdf")
    )
    
    ggsave(
      filename = filename,
      plot = plots$combined,
      width = 12,
      height = 4
    )
    
    log_message(sprintf("Saved: %s", filename), verbose = verbose)
  }
  
  log_message("Cluster QC plots completed", verbose = verbose)
  
  return(invisible(plots))
}
