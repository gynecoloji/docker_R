#!/usr/bin/env Rscript
# =============================================================================
# Test Runner for Seurat Pipeline
# =============================================================================
#
# This script runs all unit tests for the Seurat pipeline using testthat.
#
# Author: Ji (gynecoloji)
# Date: 2025-01-03
#
# Usage:
#   Rscript tests/run_tests.R
#   # Or from R:
#   source("tests/run_tests.R")
#
# =============================================================================

# Load required packages
if (!requireNamespace("testthat", quietly = TRUE)) {
  message("Installing testthat package...")
  install.packages("testthat")
}

library(testthat)

# Set working directory to project root
if (basename(getwd()) == "tests") {
  setwd("..")
}

cat("\n")
cat("================================================================================\n")
cat("                    Seurat Pipeline - Unit Tests\n")
cat("================================================================================\n")
cat("\n")

cat("Working directory:", getwd(), "\n")
cat("Test directory:   ", file.path(getwd(), "tests"), "\n")
cat("\n")

# Source all R functions
cat("Loading pipeline functions...\n")
source("R/utils.R")
source("R/qc_functions.R")
source("R/processing_functions.R")
source("R/plotting_functions.R")

cat("Functions loaded successfully.\n")
cat("\n")

# Run all tests
cat("================================================================================\n")
cat("Running tests...\n")
cat("================================================================================\n")
cat("\n")

# Run tests and capture results
results <- test_dir(
  path = "tests/testthat",
  reporter = "progress",
  stop_on_failure = FALSE
)

# Print summary
cat("\n")
cat("================================================================================\n")
cat("Test Summary\n")
cat("================================================================================\n")
cat("\n")

# Summary statistics
n_tests <- length(results)
n_passed <- sum(sapply(results, function(x) x$passed))
n_failed <- sum(sapply(results, function(x) x$failed))
n_warnings <- sum(sapply(results, function(x) x$warning))
n_skipped <- sum(sapply(results, function(x) x$skipped))

cat(sprintf("Total tests:  %d\n", n_tests))
cat(sprintf("Passed:       %d ✓\n", n_passed))
cat(sprintf("Failed:       %d ✗\n", n_failed))
cat(sprintf("Warnings:     %d ⚠\n", n_warnings))
cat(sprintf("Skipped:      %d ○\n", n_skipped))
cat("\n")

# Exit with appropriate code
if (n_failed > 0) {
  cat("Some tests FAILED. Please review the output above.\n")
  cat("\n")
  quit(status = 1)
} else {
  cat("All tests PASSED! ✓\n")
  cat("\n")
  quit(status = 0)
}