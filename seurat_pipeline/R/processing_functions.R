#' Seurat Processing and Batch Correction Functions
#'
#' Production-ready functions for single-cell RNA-seq data processing,
#' including normalization, scaling, dimensionality reduction, clustering,
#' and batch correction.
#'
#' @author Ji (gynecoloji)
#' @date 2025-01-03

# Required packages
#' @import Seurat
#' @import dplyr
#' @import future
#' @import harmony

# =============================================================================
# CONFIGURATION MANAGEMENT
# =============================================================================

#' Get Default Processing Configuration
#'
#' Returns default parameters for Seurat processing pipeline.
#' Can be overridden by user-provided values.
#'
#' @return Named list of default processing parameters
#' @export
get_default_processing_config <- function() {
  list(
    # Normalization
    normalization_method = "LogNormalize",
    scale_factor = 10000,
    
    # Feature selection
    selection_method = "vst",
    n_features = 3000,
    
    # Scaling
    vars_to_regress = c("nCount_RNA", "percent.mt"),
    scale_all_genes = TRUE,
    
    # Cell cycle
    regress_cell_cycle = TRUE,
    
    # PCA
    n_pcs = 50,
    
    # Clustering
    dims_use = 30,
    resolution = 0.5,
    algorithm = 1,  # 1 = Louvain, 2 = Louvain with multilevel refinement
    
    # UMAP/tSNE
    run_umap = TRUE,
    run_tsne = TRUE,
    umap_dims = 30,
    tsne_dims = 30,
    
    # JackStraw
    run_jackstraw = FALSE,
    jackstraw_replicates = 100,
    
    # Parallel processing
    n_cores = 1,
    future_gb = 10
  )
}

#' Merge User Config with Defaults
#'
#' @param user_config User-provided configuration list
#' @return Complete configuration with defaults filled in
#' @keywords internal
merge_config <- function(user_config = NULL) {
  default_config <- get_default_processing_config()
  
  if (is.null(user_config)) {
    return(default_config)
  }
  
  # Merge user config with defaults
  config <- default_config
  for (name in names(user_config)) {
    config[[name]] <- user_config[[name]]
  }
  
  return(config)
}

# =============================================================================
# PARALLEL PROCESSING SETUP
# =============================================================================

#' Setup Parallel Processing
#'
#' Configures the future package for parallel processing in Seurat.
#'
#' @param n_cores Number of cores to use (default: 1, no parallelization)
#' @param memory_gb Memory limit per worker in GB (default: 10)
#' @param verbose Logical, print setup messages
#' @keywords internal
setup_parallel <- function(n_cores = 1, memory_gb = 10, verbose = TRUE) {
  
  if (n_cores > 1) {
    if (!requireNamespace("future", quietly = TRUE)) {
      warning("Package 'future' not available. Running sequentially.")
      return(invisible(NULL))
    }
    
    # Get available cores
    max_cores <- future::availableCores()
    
    if (n_cores > max_cores) {
      warning(sprintf(
        "Requested %d cores but only %d available. Using %d cores.",
        n_cores, max_cores, max_cores
      ))
      n_cores <- max_cores
    }
    
    # Setup multicore processing
    future::plan("multicore", workers = n_cores)
    
    # Set memory limits
    options(future.globals.maxSize = memory_gb * 1024^3)
    
    if (verbose) {
      message(sprintf(
        "Parallel processing enabled: %d cores, %d GB memory per worker",
        n_cores, memory_gb
      ))
    }
  } else {
    future::plan("sequential")
    if (verbose) {
      message("Sequential processing enabled")
    }
  }
  
  invisible(NULL)
}

# =============================================================================
# CORE PROCESSING PIPELINE
# =============================================================================

#' Process Seurat Object - Production Pipeline
#'
#' Comprehensive processing pipeline for single-cell RNA-seq data including
#' normalization, scaling, dimensionality reduction, and clustering.
#' Supports parallel processing and flexible configuration.
#'
#' @param seurat_obj Seurat object (pre-filtered)
#' @param config Named list of processing parameters (see get_default_processing_config)
#' @param species Species for cell cycle genes ("human" or "mouse")
#' @param verbose Logical, print progress messages
#'
#' @return Processed Seurat object with normalized data, scaled data, PCA, clusters, and UMAP/tSNE
#'
#' @details
#' The pipeline performs the following steps:
#' \enumerate{
#'   \item Quality metrics calculation (if not present)
#'   \item Normalization
#'   \item Feature selection
#'   \item Cell cycle scoring (optional)
#'   \item Data scaling with regression
#'   \item PCA
#'   \item JackStraw analysis (optional)
#'   \item Clustering
#'   \item UMAP/tSNE
#' }
#'
#' @examples
#' \dontrun{
#' # Using defaults
#' seurat_obj <- process_seurat_object(seurat_obj)
#'
#' # Custom configuration
#' config <- list(
#'   n_features = 2000,
#'   dims_use = 20,
#'   resolution = 0.8,
#'   n_cores = 4
#' )
#' seurat_obj <- process_seurat_object(seurat_obj, config = config)
#' }
#'
#' @export
process_seurat_object <- function(seurat_obj,
                                  config = NULL,
                                  species = "human",
                                  verbose = TRUE) {
  
  # Validate input
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object")
  }
  
  if (ncol(seurat_obj) == 0) {
    stop("Seurat object contains zero cells")
  }
  
  # Merge configuration
  config <- merge_config(config)
  
  # Setup parallel processing
  setup_parallel(
    n_cores = config$n_cores,
    memory_gb = config$future_gb,
    verbose = verbose
  )
  
  log_message(
    sprintf("Starting Seurat processing pipeline: %d cells, %d features",
            ncol(seurat_obj), nrow(seurat_obj)),
    verbose = verbose
  )
  
  # ============================================
  # STEP 1: Calculate QC metrics if not present
  # ============================================
  
  if (!"percent.mt" %in% colnames(seurat_obj@meta.data)) {
    log_message("Calculating mitochondrial percentage", verbose = verbose)
    mito_pattern <- if (species == "human") "^MT-" else "^mt-"
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
      object = seurat_obj,
      pattern = mito_pattern
    )
  }
  
  if (!"percent.rb" %in% colnames(seurat_obj@meta.data)) {
    log_message("Calculating ribosomal percentage", verbose = verbose)
    seurat_obj[["percent.rb"]] <- PercentageFeatureSet(
      object = seurat_obj,
      pattern = "^RP[SL]"
    )
  }
  
  # ============================================
  # STEP 2: Normalization
  # ============================================
  
  log_message(
    sprintf("Normalizing data: method = %s, scale factor = %d",
            config$normalization_method, config$scale_factor),
    verbose = verbose
  )
  
  seurat_obj <- NormalizeData(
    object = seurat_obj,
    normalization.method = config$normalization_method,
    scale.factor = config$scale_factor,
    verbose = FALSE
  )
  
  # ============================================
  # STEP 3: Feature Selection
  # ============================================
  
  log_message(
    sprintf("Finding variable features: method = %s, n = %d",
            config$selection_method, config$n_features),
    verbose = verbose
  )
  
  seurat_obj <- FindVariableFeatures(
    object = seurat_obj,
    selection.method = config$selection_method,
    nfeatures = config$n_features,
    verbose = FALSE
  )
  
  n_var_features <- length(VariableFeatures(seurat_obj))
  log_message(
    sprintf("Identified %d variable features", n_var_features),
    verbose = verbose
  )
  
  # ============================================
  # STEP 4: Cell Cycle Scoring
  # ============================================
  
  if (config$regress_cell_cycle) {
    log_message("Performing cell cycle scoring", verbose = verbose)
    
    tryCatch({
      # Get cell cycle genes
      if (species == "human") {
        s_genes <- cc.genes$s.genes
        g2m_genes <- cc.genes$g2m.genes
      } else if (species == "mouse") {
        # Convert human genes to mouse
        s_genes <- tolower(cc.genes$s.genes)
        s_genes <- paste0(toupper(substr(s_genes, 1, 1)), substr(s_genes, 2, nchar(s_genes)))
        
        g2m_genes <- tolower(cc.genes$g2m.genes)
        g2m_genes <- paste0(toupper(substr(g2m_genes, 1, 1)), substr(g2m_genes, 2, nchar(g2m_genes)))
      } else {
        stop("Species must be 'human' or 'mouse'")
      }
      
      seurat_obj <- CellCycleScoring(
        object = seurat_obj,
        s.features = s_genes,
        g2m.features = g2m_genes,
        set.ident = FALSE
      )
      
      # Add cell cycle scores to regression variables
      if (config$regress_cell_cycle) {
        config$vars_to_regress <- unique(c(
          config$vars_to_regress,
          "S.Score",
          "G2M.Score"
        ))
      }
      
      log_message("Cell cycle scoring completed", verbose = verbose)
      
    }, error = function(e) {
      warning("Cell cycle scoring failed: ", e$message)
      log_message("Skipping cell cycle scoring", level = "WARNING", verbose = verbose)
    })
  }
  
  # ============================================
  # STEP 5: Scaling
  # ============================================
  
  log_message(
    sprintf("Scaling data: regressing %s",
            paste(config$vars_to_regress, collapse = ", ")),
    verbose = verbose
  )
  
  # Determine features to scale
  features_to_scale <- if (config$scale_all_genes) {
    rownames(seurat_obj)
  } else {
    VariableFeatures(seurat_obj)
  }
  
  seurat_obj <- ScaleData(
    object = seurat_obj,
    features = features_to_scale,
    vars.to.regress = config$vars_to_regress,
    verbose = FALSE
  )
  
  log_message(
    sprintf("Scaled %d features", length(features_to_scale)),
    verbose = verbose
  )
  
  # ============================================
  # STEP 6: PCA
  # ============================================
  
  log_message(
    sprintf("Running PCA: %d components", config$n_pcs),
    verbose = verbose
  )
  
  seurat_obj <- RunPCA(
    object = seurat_obj,
    features = VariableFeatures(object = seurat_obj),
    npcs = config$n_pcs,
    verbose = FALSE
  )
  
  # ============================================
  # STEP 7: JackStraw (Optional)
  # ============================================
  
  if (config$run_jackstraw) {
    log_message(
      sprintf("Running JackStraw analysis: %d replicates",
              config$jackstraw_replicates),
      verbose = verbose
    )
    
    seurat_obj <- JackStraw(
      object = seurat_obj,
      num.replicate = config$jackstraw_replicates,
      dims = config$n_pcs
    )
    
    seurat_obj <- ScoreJackStraw(
      object = seurat_obj,
      dims = 1:config$n_pcs
    )
    
    log_message("JackStraw analysis completed", verbose = verbose)
  }
  
  # ============================================
  # STEP 8: Clustering
  # ============================================
  
  log_message(
    sprintf("Finding neighbors: using %d dimensions", config$dims_use),
    verbose = verbose
  )
  
  seurat_obj <- FindNeighbors(
    object = seurat_obj,
    dims = 1:config$dims_use,
    verbose = FALSE
  )
  
  log_message(
    sprintf("Finding clusters: resolution = %.2f, algorithm = %d",
            config$resolution, config$algorithm),
    verbose = verbose
  )
  
  seurat_obj <- FindClusters(
    object = seurat_obj,
    resolution = config$resolution,
    algorithm = config$algorithm,
    verbose = FALSE
  )
  
  n_clusters <- length(levels(Idents(seurat_obj)))
  log_message(
    sprintf("Identified %d clusters", n_clusters),
    verbose = verbose
  )
  
  # ============================================
  # STEP 9: UMAP
  # ============================================
  
  if (config$run_umap) {
    log_message(
      sprintf("Running UMAP: %d dimensions", config$umap_dims),
      verbose = verbose
    )
    
    seurat_obj <- RunUMAP(
      object = seurat_obj,
      dims = 1:config$umap_dims,
      verbose = FALSE
    )
  }
  
  # ============================================
  # STEP 10: tSNE
  # ============================================
  
  if (config$run_tsne) {
    log_message(
      sprintf("Running tSNE: %d dimensions", config$tsne_dims),
      verbose = verbose
    )
    
    seurat_obj <- RunTSNE(
      object = seurat_obj,
      dims = 1:config$tsne_dims,
      verbose = FALSE
    )
  }
  
  # ============================================
  # COMPLETE
  # ============================================
  
  log_message("Seurat processing pipeline completed successfully", verbose = verbose)
  
  # Add processing metadata
  seurat_obj@misc$processing_config <- config
  seurat_obj@misc$processing_timestamp <- Sys.time()
  
  return(seurat_obj)
}

# =============================================================================
# BATCH CORRECTION
# =============================================================================

#' Get Default Batch Correction Configuration
#'
#' @return Named list of default batch correction parameters
#' @export
get_default_batch_config <- function() {
  list(
    # Method
    method = "harmony",
    
    # Dimensions
    dims = 30,
    
    # Clustering
    resolution = 0.5,
    
    # Regression variables
    vars_to_regress = c("percent.mt"),
    
    # Method-specific parameters
    harmony_theta = NULL,
    harmony_lambda = NULL,
    harmony_max_iter = 10,
    
    # CCA/RPCA parameters
    k_anchor = 5,
    k_filter = 200,
    k_score = 30,
    
    # SCTransform parameters
    sct_n_genes = 3000,
    sct_vars_to_regress = c("percent.mt"),
    
    # Output
    run_umap = TRUE,
    run_tsne = TRUE,
    
    # Parallel
    n_cores = 1,
    future_gb = 10
  )
}

#' Process Seurat Object with Batch Correction
#'
#' Production-ready batch correction pipeline supporting multiple methods:
#' Harmony, CCA, RPCA, and SCTransform.
#'
#' @param seurat_obj Seurat object (should be pre-filtered and have basic QC metrics)
#' @param batch_var Column name in metadata specifying batch assignment
#' @param config Named list of batch correction parameters (see get_default_batch_config)
#' @param species Species for gene patterns ("human" or "mouse")
#' @param verbose Logical, print progress messages
#'
#' @return Seurat object with batch correction applied
#'
#' @details
#' Supported methods:
#' \itemize{
#'   \item \strong{harmony}: Fast, recommended for most use cases
#'   \item \strong{CCA}: Canonical Correlation Analysis integration
#'   \item \strong{RPCA}: Reciprocal PCA integration
#'   \item \strong{SCTransform}: SCTransform-based normalization (not integration per se)
#' }
#'
#' @examples
#' \dontrun{
#' # Using Harmony (default)
#' seurat_corrected <- process_batch_correction(
#'   seurat_obj,
#'   batch_var = "orig.ident"
#' )
#'
#' # Using CCA with custom config
#' config <- list(
#'   method = "CCA",
#'   dims = 20,
#'   resolution = 0.8
#' )
#' seurat_corrected <- process_batch_correction(
#'   seurat_obj,
#'   batch_var = "batch",
#'   config = config
#' )
#' }
#'
#' @export
process_batch_correction <- function(seurat_obj,
                                     batch_var,
                                     config = NULL,
                                     species = "human",
                                     verbose = TRUE) {
  
  # ============================================
  # VALIDATION
  # ============================================
  
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object")
  }
  
  if (missing(batch_var) || is.null(batch_var)) {
    stop("batch_var is required")
  }
  
  if (!batch_var %in% colnames(seurat_obj@meta.data)) {
    stop(sprintf("Batch variable '%s' not found in metadata", batch_var))
  }
  
  # Check number of batches
  n_batches <- length(unique(seurat_obj@meta.data[[batch_var]]))
  
  if (n_batches < 2) {
    warning(sprintf(
      "Only %d batch found in '%s'. Batch correction may not be necessary.",
      n_batches, batch_var
    ))
  }
  
  log_message(
    sprintf("Detected %d batches in '%s'", n_batches, batch_var),
    verbose = verbose
  )
  
  # Merge configuration
  config <- modifyList(get_default_batch_config(), config %||% list())
  
  # Validate method
  valid_methods <- c("harmony", "CCA", "RPCA", "SCTransform")
  if (!config$method %in% valid_methods) {
    stop(sprintf(
      "Invalid method '%s'. Must be one of: %s",
      config$method,
      paste(valid_methods, collapse = ", ")
    ))
  }
  
  log_message(
    sprintf("Starting batch correction: method = %s", config$method),
    verbose = verbose
  )
  
  # Setup parallel processing
  setup_parallel(
    n_cores = config$n_cores,
    memory_gb = config$future_gb,
    verbose = verbose
  )
  
  # ============================================
  # PRE-PROCESSING
  # ============================================
  
  log_message("Step 1: Pre-processing", verbose = verbose)
  
  # Calculate QC metrics if not present
  if (!"percent.mt" %in% colnames(seurat_obj@meta.data)) {
    mito_pattern <- if (species == "human") "^MT-" else "^mt-"
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
      object = seurat_obj,
      pattern = mito_pattern
    )
  }
  
  if (!"percent.rb" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj[["percent.rb"]] <- PercentageFeatureSet(
      object = seurat_obj,
      pattern = "^RP[SL]"
    )
  }
  
  # ============================================
  # METHOD-SPECIFIC PROCESSING
  # ============================================
  
  if (config$method == "SCTransform") {
    
    # ========================================
    # SCTransform Method
    # ========================================
    
    log_message("Step 2: SCTransform normalization", verbose = verbose)
    
    if (!requireNamespace("sctransform", quietly = TRUE)) {
      stop("Package 'sctransform' required for SCTransform method")
    }
    
    seurat_obj <- SCTransform(
      object = seurat_obj,
      vars.to.regress = config$sct_vars_to_regress,
      variable.features.n = config$sct_n_genes,
      verbose = FALSE
    )
    
    log_message("Step 3: PCA on SCTransform data", verbose = verbose)
    
    seurat_obj <- RunPCA(
      object = seurat_obj,
      npcs = config$dims,
      verbose = FALSE
    )
    
  } else {
    
    # ========================================
    # Standard Normalization for other methods
    # ========================================
    
    log_message("Step 2: Standard normalization", verbose = verbose)
    
    seurat_obj <- NormalizeData(
      object = seurat_obj,
      verbose = FALSE
    )
    
    seurat_obj <- FindVariableFeatures(
      object = seurat_obj,
      selection.method = "vst",
      nfeatures = 2000,
      verbose = FALSE
    )
    
    log_message("Step 3: Scaling and PCA", verbose = verbose)
    
    all_genes <- rownames(seurat_obj)
    seurat_obj <- ScaleData(
      object = seurat_obj,
      features = all_genes,
      vars.to.regress = config$vars_to_regress,
      verbose = FALSE
    )
    
    seurat_obj <- RunPCA(
      object = seurat_obj,
      features = VariableFeatures(object = seurat_obj),
      npcs = config$dims,
      verbose = FALSE
    )
  }
  
  # ============================================
  # APPLY BATCH CORRECTION
  # ============================================
  
  log_message(
    sprintf("Step 4: Applying %s batch correction", config$method),
    verbose = verbose
  )
  
  if (config$method == "harmony") {
    
    # ========================================
    # Harmony Integration
    # ========================================
    
    if (!requireNamespace("harmony", quietly = TRUE)) {
      stop("Package 'harmony' required for Harmony method. Install with: install.packages('harmony')")
    }
    
    seurat_obj <- harmony::RunHarmony(
      object = seurat_obj,
      group.by.vars = batch_var,
      dims.use = 1:config$dims,
      theta = config$harmony_theta,
      lambda = config$harmony_lambda,
      max.iter.harmony = config$harmony_max_iter,
      verbose = FALSE
    )
    
    reduction_use <- "harmony"
    
  } else if (config$method %in% c("CCA", "RPCA")) {
    
    # ========================================
    # CCA or RPCA Integration
    # ========================================
    
    log_message("Splitting object by batch", verbose = verbose)
    
    obj_list <- SplitObject(seurat_obj, split.by = batch_var)
    
    log_message(
      sprintf("Split into %d batches", length(obj_list)),
      verbose = verbose
    )
    
    # Find integration anchors
    if (config$method == "CCA") {
      
      log_message("Finding CCA anchors", verbose = verbose)
      
      anchors <- FindIntegrationAnchors(
        object.list = obj_list,
        dims = 1:config$dims,
        k.anchor = config$k_anchor,
        k.filter = config$k_filter,
        k.score = config$k_score,
        verbose = FALSE
      )
      
    } else {  # RPCA
      
      log_message("Running PCA on each batch", verbose = verbose)
      
      obj_list <- lapply(obj_list, function(x) {
        x <- RunPCA(x, features = VariableFeatures(x), verbose = FALSE)
        return(x)
      })
      
      log_message("Finding RPCA anchors", verbose = verbose)
      
      anchors <- FindIntegrationAnchors(
        object.list = obj_list,
        reduction = "rpca",
        dims = 1:config$dims,
        k.anchor = config$k_anchor,
        k.filter = config$k_filter,
        k.score = config$k_score,
        verbose = FALSE
      )
    }
    
    log_message("Integrating data", verbose = verbose)
    
    seurat_obj <- IntegrateData(
      anchorset = anchors,
      dims = 1:config$dims,
      verbose = FALSE
    )
    
    # Set default assay
    DefaultAssay(seurat_obj) <- "integrated"
    
    # Scale integrated data
    log_message("Scaling integrated data", verbose = verbose)
    
    seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
    
    # Run PCA on integrated data
    seurat_obj <- RunPCA(seurat_obj, npcs = config$dims, verbose = FALSE)
    
    reduction_use <- "pca"
  }
  
  # ============================================
  # CLUSTERING
  # ============================================
  
  log_message("Step 5: Clustering", verbose = verbose)
  
  seurat_obj <- FindNeighbors(
    object = seurat_obj,
    reduction = reduction_use,
    dims = 1:config$dims,
    verbose = FALSE
  )
  
  seurat_obj <- FindClusters(
    object = seurat_obj,
    resolution = config$resolution,
    verbose = FALSE
  )
  
  n_clusters <- length(levels(Idents(seurat_obj)))
  log_message(
    sprintf("Identified %d clusters", n_clusters),
    verbose = verbose
  )
  
  # ============================================
  # DIMENSIONALITY REDUCTION
  # ============================================
  
  if (config$run_umap) {
    log_message("Step 6: Running UMAP", verbose = verbose)
    
    seurat_obj <- RunUMAP(
      object = seurat_obj,
      reduction = reduction_use,
      dims = 1:config$dims,
      verbose = FALSE
    )
  }
  
  if (config$run_tsne) {
    log_message("Step 7: Running tSNE", verbose = verbose)
    
    seurat_obj <- RunTSNE(
      object = seurat_obj,
      reduction = reduction_use,
      dims = 1:config$dims,
      verbose = FALSE
    )
  }
  
  # ============================================
  # ADD METADATA
  # ============================================
  
  seurat_obj@meta.data$batch_correction_method <- config$method
  seurat_obj@meta.data$batch_variable <- batch_var
  seurat_obj@misc$batch_correction_config <- config
  seurat_obj@misc$batch_correction_timestamp <- Sys.time()
  
  log_message("Batch correction completed successfully", verbose = verbose)
  log_message(
    sprintf("Final object: %d cells, %d features, %d clusters",
            ncol(seurat_obj), nrow(seurat_obj), n_clusters),
    verbose = verbose
  )
  
  return(seurat_obj)
}

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

#' Get Processing Summary
#'
#' Extract summary statistics from a processed Seurat object
#'
#' @param seurat_obj Processed Seurat object
#' @return Data frame with processing summary
#' @export
get_processing_summary <- function(seurat_obj) {
  
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object")
  }
  
  summary_df <- data.frame(
    n_cells = ncol(seurat_obj),
    n_features = nrow(seurat_obj),
    n_variable_features = length(VariableFeatures(seurat_obj)),
    stringsAsFactors = FALSE
  )
  
  # Add cluster information if available
  if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    summary_df$n_clusters <- length(levels(seurat_obj$seurat_clusters))
  }
  
  # Add batch correction info if available
  if ("batch_correction_method" %in% colnames(seurat_obj@meta.data)) {
    summary_df$batch_correction_method <- unique(seurat_obj$batch_correction_method)[1]
  }
  
  # Add processing timestamp if available
  if (!is.null(seurat_obj@misc$processing_timestamp)) {
    summary_df$processing_timestamp <- as.character(seurat_obj@misc$processing_timestamp)
  }
  
  return(summary_df)
}
