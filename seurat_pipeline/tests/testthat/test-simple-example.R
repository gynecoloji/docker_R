# =============================================================================
# SIMPLE TEST EXAMPLE - Start Here!
# =============================================================================

library(testthat)
library(Seurat)

# Load the functions you want to test
source("../../R/utils.R", chdir = TRUE)

# =============================================================================
# EXAMPLE 1: Test a Simple Function
# =============================================================================

# The function we're testing: validate_seurat_object()
# What it does: Checks if input is a valid Seurat object

test_that("validate_seurat_object works with valid Seurat object", {
  
  # STEP 1: Create test data
  counts <- matrix(rpois(100, lambda = 5), nrow = 10, ncol = 10)
  rownames(counts) <- paste0("Gene_", 1:10)
  colnames(counts) <- paste0("Cell_", 1:10)
  seurat_obj <- CreateSeuratObject(counts = counts)
  
  # STEP 2: Run the function
  # Should NOT error
  result <- validate_seurat_object(seurat_obj)
  
  # STEP 3: Check the result
  expect_true(result)  # Should return TRUE
})

test_that("validate_seurat_object rejects invalid input", {
  
  # STEP 1: Create invalid data (not a Seurat object)
  invalid_input <- list(a = 1, b = 2)
  
  # STEP 2 & 3: Should throw an error
  expect_error(
    validate_seurat_object(invalid_input),
    "must be a Seurat object"
  )
})

# =============================================================================
# EXAMPLE 2: Test a Function That Modifies Data
# =============================================================================

# The function we're testing: add_qc_metrics()
# What it does: Adds QC columns to Seurat object metadata

test_that("add_qc_metrics adds new columns", {
  
  # STEP 1: Create test Seurat object
  counts <- matrix(rpois(100, lambda = 5), nrow = 10, ncol = 10)
  gene_names <- paste0("Gene_", 1:10)
  gene_names[1:3] <- paste0("MT-", 1:3)  # Add some mito genes
  rownames(counts) <- gene_names
  colnames(counts) <- paste0("Cell_", 1:10)
  seurat_obj <- CreateSeuratObject(counts = counts)
  
  # STEP 2: Run the function
  seurat_obj <- add_qc_metrics(seurat_obj, verbose = FALSE)
  
  # STEP 3: Check that new columns were added
  expect_true("mitoRatio" %in% colnames(seurat_obj@meta.data))
  expect_true("riboRatio" %in% colnames(seurat_obj@meta.data))
  expect_true("log10GenesPerUMI" %in% colnames(seurat_obj@meta.data))
})

test_that("add_qc_metrics calculates reasonable values", {
  
  # STEP 1: Create test data
  counts <- matrix(rpois(100, lambda = 5), nrow = 10, ncol = 10)
  gene_names <- paste0("Gene_", 1:10)
  gene_names[1:3] <- paste0("MT-", 1:3)
  rownames(counts) <- gene_names
  colnames(counts) <- paste0("Cell_", 1:10)
  seurat_obj <- CreateSeuratObject(counts = counts)
  
  # STEP 2: Run the function
  seurat_obj <- add_qc_metrics(seurat_obj, verbose = FALSE)
  
  # STEP 3: Check values are reasonable
  expect_true(all(seurat_obj$mitoRatio >= 0))
  expect_true(all(seurat_obj$mitoRatio <= 100))
  expect_true(all(seurat_obj$log10GenesPerUMI > 0, na.rm = TRUE))
})

# =============================================================================
# EXAMPLE 3: Test a Function That Saves Files
# =============================================================================

# The function we're testing: save_seurat_object()
# What it does: Saves Seurat object to .rds file

test_that("save_seurat_object creates file", {
  
  # STEP 1: Create test data
  counts <- matrix(rpois(100, lambda = 5), nrow = 10, ncol = 10)
  rownames(counts) <- paste0("Gene_", 1:10)
  colnames(counts) <- paste0("Cell_", 1:10)
  seurat_obj <- CreateSeuratObject(counts = counts)
  
  # Create temporary file path
  temp_file <- tempfile(fileext = ".rds")
  
  # STEP 2: Run the function
  save_seurat_object(seurat_obj, temp_file, verbose = FALSE)
  
  # STEP 3: Check file was created
  expect_true(file.exists(temp_file))
  
  # CLEANUP: Always delete temporary files!
  unlink(temp_file)
})

# =============================================================================
# EXAMPLE 4: Test a Function With Configuration
# =============================================================================

# The function we're testing: get_default_processing_config()
# What it does: Returns default configuration list

test_that("get_default_processing_config returns valid config", {
  
  # STEP 1: No setup needed
  
  # STEP 2: Run the function
  config <- get_default_processing_config()
  
  # STEP 3: Check the result
  expect_type(config, "list")  # Should be a list
  expect_true("n_features" %in% names(config))  # Has required keys
  expect_true("resolution" %in% names(config))
  expect_true(config$n_features > 0)  # Values are positive
  expect_true(config$resolution >= 0)
})

