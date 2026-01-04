#' Utility Functions for Seurat Pipeline
#'
#' Helper functions for logging, validation, file operations, and package management
#'
#' @author Ji (gynecoloji)
#' @date 2025-01-03

# =============================================================================
# PACKAGE MANAGEMENT
# =============================================================================

#' Check and Load Required Packages
#'
#' Checks if packages are installed, installs missing ones, and loads them.
#' Production-ready with proper error handling and BioConductor support.
#'
#' @param packages Character vector of package names
#' @param bioc_packages Character vector of Bioconductor package names
#' @param verbose Logical, print installation messages
#'
#' @return NULL (invisible)
#' @export
check_and_load_packages <- function(packages = NULL,
                                    bioc_packages = NULL,
                                    verbose = TRUE) {
  
  # CRAN packages
  if (!is.null(packages) && length(packages) > 0) {
    for (pkg in packages) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        if (verbose) {
          message(sprintf("Installing CRAN package: %s", pkg))
        }
        
        tryCatch({
          install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org/")
          
          if (verbose) {
            message(sprintf("Successfully installed: %s", pkg))
          }
        }, error = function(e) {
          stop(sprintf("Failed to install package '%s': %s", pkg, e$message))
        })
      }
      
      # Load package
      library(pkg, character.only = TRUE, quietly = !verbose)
    }
  }
  
  # Bioconductor packages
  if (!is.null(bioc_packages) && length(bioc_packages) > 0) {
    
    # Check if BiocManager is available
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      if (verbose) {
        message("Installing BiocManager")
      }
      install.packages("BiocManager", repos = "https://cloud.r-project.org/")
    }
    
    for (pkg in bioc_packages) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        if (verbose) {
          message(sprintf("Installing Bioconductor package: %s", pkg))
        }
        
        tryCatch({
          BiocManager::install(pkg, update = FALSE, ask = FALSE)
          
          if (verbose) {
            message(sprintf("Successfully installed: %s", pkg))
          }
        }, error = function(e) {
          stop(sprintf("Failed to install Bioconductor package '%s': %s", pkg, e$message))
        })
      }
      
      # Load package
      library(pkg, character.only = TRUE, quietly = !verbose)
    }
  }
  
  invisible(NULL)
}

#' Get Required Packages List
#'
#' Returns list of packages required for the pipeline
#'
#' @param include_optional Logical, include optional packages
#' @return List with 'cran' and 'bioconductor' package vectors
#' @export
get_required_packages <- function(include_optional = TRUE) {
  
  required_cran <- c(
    "Seurat",
    "dplyr",
    "stringr",
    "ggplot2",
    "ggrepel",
    "patchwork",
    "RColorBrewer",
    "pheatmap",
    "tidyr",
    "future",
    "yaml",
    "jsonlite"
  )
  
  required_bioc <- c(
    "SingleR"
  )
  
  optional_cran <- c(
    "harmony",
    "msigdbr",
    "grid",
    "gridExtra",
    "pals"
  )
  
  optional_bioc <- c(
    "batchelor"
  )
  
  if (include_optional) {
    required_cran <- c(required_cran, optional_cran)
    required_bioc <- c(required_bioc, optional_bioc)
  }
  
  return(list(
    cran = required_cran,
    bioconductor = required_bioc
  ))
}

# =============================================================================
# LOGGING FUNCTIONS
# =============================================================================

#' Log Message with Timestamp
#'
#' Prints formatted log message with timestamp and level
#'
#' @param msg Message to log
#' @param level Log level ("INFO", "WARNING", "ERROR", "SUCCESS")
#' @param verbose Logical, whether to print
#' @param log_file Optional file path to append log messages
#'
#' @return NULL (invisible)
#' @export
log_message <- function(msg,
                       level = "INFO",
                       verbose = TRUE,
                       log_file = NULL) {
  
  if (!verbose && is.null(log_file)) {
    return(invisible(NULL))
  }
  
  # Format message
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  formatted_msg <- sprintf("[%s] %s: %s", timestamp, level, msg)
  
  # Print to console
  if (verbose) {
    if (level == "ERROR") {
      message(crayon::red(formatted_msg))
    } else if (level == "WARNING") {
      message(crayon::yellow(formatted_msg))
    } else if (level == "SUCCESS") {
      message(crayon::green(formatted_msg))
    } else {
      message(formatted_msg)
    }
  }
  
  # Write to log file
  if (!is.null(log_file)) {
    write(formatted_msg, file = log_file, append = TRUE)
  }
  
  invisible(NULL)
}

#' Create Logger
#'
#' Creates a logger function with pre-configured settings
#'
#' @param log_file Path to log file
#' @param verbose Logical, print to console
#' @return Logger function
#' @export
create_logger <- function(log_file = NULL, verbose = TRUE) {
  
  # Initialize log file
  if (!is.null(log_file)) {
    log_dir <- dirname(log_file)
    if (!dir.exists(log_dir)) {
      dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    # Write header
    header <- sprintf(
      "=== Log started at %s ===\n",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    )
    write(header, file = log_file, append = FALSE)
  }
  
  # Return logger function
  function(msg, level = "INFO") {
    log_message(msg, level, verbose, log_file)
  }
}

# =============================================================================
# FILE OPERATIONS
# =============================================================================

#' Create Output Directory
#'
#' Creates directory if it doesn't exist, with error handling
#'
#' @param path Directory path to create
#' @param verbose Logical, print messages
#' @return Character, normalized path to created directory
#' @export
create_output_dir <- function(path, verbose = TRUE) {
  
  if (is.null(path) || path == "") {
    stop("Path cannot be NULL or empty")
  }
  
  if (!dir.exists(path)) {
    tryCatch({
      dir.create(path, recursive = TRUE, showWarnings = FALSE)
      
      if (verbose) {
        log_message(sprintf("Created directory: %s", path), "INFO", verbose)
      }
    }, error = function(e) {
      stop(sprintf("Failed to create directory '%s': %s", path, e$message))
    })
  }
  
  return(normalizePath(path, mustWork = FALSE))
}

#' Save Seurat Object
#'
#' Saves Seurat object with error handling and optional compression
#'
#' @param seurat_obj Seurat object to save
#' @param file_path Output file path (should end in .rds or .RDS)
#' @param compress Logical or character, compression type
#' @param verbose Logical, print messages
#'
#' @return NULL (invisible)
#' @export
save_seurat_object <- function(seurat_obj,
                               file_path,
                               compress = "gzip",
                               verbose = TRUE) {
  
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object")
  }
  
  if (is.null(file_path) || file_path == "") {
    stop("file_path cannot be NULL or empty")
  }
  
  # Create directory if needed
  output_dir <- dirname(file_path)
  create_output_dir(output_dir, verbose = FALSE)
  
  # Save object
  if (verbose) {
    log_message(sprintf("Saving Seurat object to: %s", file_path), "INFO", verbose)
  }
  
  tryCatch({
    saveRDS(seurat_obj, file = file_path, compress = compress)
    
    # Get file size
    file_size <- file.size(file_path)
    size_mb <- round(file_size / 1024^2, 2)
    
    if (verbose) {
      log_message(
        sprintf("Successfully saved (%s MB)", size_mb),
        "SUCCESS",
        verbose
      )
    }
  }, error = function(e) {
    stop(sprintf("Failed to save Seurat object: %s", e$message))
  })
  
  invisible(NULL)
}

#' Load Seurat Object
#'
#' Loads Seurat object with error handling
#'
#' @param file_path Path to .rds file
#' @param verbose Logical, print messages
#'
#' @return Seurat object
#' @export
load_seurat_object <- function(file_path, verbose = TRUE) {
  
  if (!file.exists(file_path)) {
    stop(sprintf("File not found: %s", file_path))
  }
  
  if (verbose) {
    log_message(sprintf("Loading Seurat object from: %s", file_path), "INFO", verbose)
  }
  
  tryCatch({
    seurat_obj <- readRDS(file_path)
    
    if (!inherits(seurat_obj, "Seurat")) {
      stop("Loaded object is not a Seurat object")
    }
    
    if (verbose) {
      log_message(
        sprintf("Successfully loaded (%d cells, %d features)",
                ncol(seurat_obj), nrow(seurat_obj)),
        "SUCCESS",
        verbose
      )
    }
    
    return(seurat_obj)
    
  }, error = function(e) {
    stop(sprintf("Failed to load Seurat object: %s", e$message))
  })
}

# =============================================================================
# VALIDATION FUNCTIONS
# =============================================================================

#' Validate Seurat Object
#'
#' Comprehensive validation of Seurat object structure and contents
#'
#' @param seurat_obj Object to validate
#' @param require_counts Logical, require counts data
#' @param require_normalized Logical, require normalized data
#' @param require_scaled Logical, require scaled data
#' @param min_cells Minimum number of cells required
#' @param min_features Minimum number of features required
#'
#' @return NULL (throws error if invalid)
#' @export
validate_seurat_object <- function(seurat_obj,
                                   require_counts = TRUE,
                                   require_normalized = FALSE,
                                   require_scaled = FALSE,
                                   min_cells = 1,
                                   min_features = 1) {
  
  # Check object type
  if (!inherits(seurat_obj, "Seurat")) {
    stop(sprintf("Input must be a Seurat object. Got: %s", class(seurat_obj)[1]))
  }
  
  # Check dimensions
  if (ncol(seurat_obj) < min_cells) {
    stop(sprintf(
      "Seurat object has %d cells, minimum required: %d",
      ncol(seurat_obj), min_cells
    ))
  }
  
  if (nrow(seurat_obj) < min_features) {
    stop(sprintf(
      "Seurat object has %d features, minimum required: %d",
      nrow(seurat_obj), min_features
    ))
  }
  
  # Check for RNA assay
  if (!"RNA" %in% names(seurat_obj@assays)) {
    stop("Seurat object must contain 'RNA' assay")
  }
  
  # Check for counts
  if (require_counts) {
    if (is.null(seurat_obj@assays$RNA@counts) || 
        length(seurat_obj@assays$RNA@counts) == 0) {
      stop("Counts data not found in RNA assay")
    }
  }
  
  # Check for normalized data
  if (require_normalized) {
    if (is.null(seurat_obj@assays$RNA@data) || 
        length(seurat_obj@assays$RNA@data) == 0) {
      stop("Normalized data not found in RNA assay")
    }
  }
  
  # Check for scaled data
  if (require_scaled) {
    if (is.null(seurat_obj@assays$RNA@scale.data) || 
        length(seurat_obj@assays$RNA@scale.data) == 0) {
      stop("Scaled data not found in RNA assay")
    }
  }
  
  invisible(NULL)
}

# =============================================================================
# CONFIGURATION MANAGEMENT
# =============================================================================

#' Load Configuration from YAML
#'
#' Loads and validates configuration from YAML file
#'
#' @param config_file Path to YAML configuration file
#' @param verbose Logical, print messages
#'
#' @return List with configuration parameters
#' @export
load_config <- function(config_file, verbose = TRUE) {
  
  if (!file.exists(config_file)) {
    stop(sprintf("Configuration file not found: %s", config_file))
  }
  
  if (verbose) {
    log_message(sprintf("Loading configuration from: %s", config_file), "INFO", verbose)
  }
  
  tryCatch({
    config <- yaml::read_yaml(config_file)
    
    if (verbose) {
      log_message("Configuration loaded successfully", "SUCCESS", verbose)
    }
    
    return(config)
    
  }, error = function(e) {
    stop(sprintf("Failed to load configuration: %s", e$message))
  })
}

#' Save Configuration to YAML
#'
#' Saves configuration list to YAML file
#'
#' @param config Configuration list
#' @param config_file Output file path
#' @param verbose Logical, print messages
#'
#' @return NULL (invisible)
#' @export
save_config <- function(config, config_file, verbose = TRUE) {
  
  # Create directory if needed
  config_dir <- dirname(config_file)
  create_output_dir(config_dir, verbose = FALSE)
  
  if (verbose) {
    log_message(sprintf("Saving configuration to: %s", config_file), "INFO", verbose)
  }
  
  tryCatch({
    yaml::write_yaml(config, config_file)
    
    if (verbose) {
      log_message("Configuration saved successfully", "SUCCESS", verbose)
    }
    
  }, error = function(e) {
    stop(sprintf("Failed to save configuration: %s", e$message))
  })
  
  invisible(NULL)
}

# =============================================================================
# SYSTEM INFORMATION
# =============================================================================

#' Get System Information
#'
#' Collects system and session information for reproducibility
#'
#' @return List with system information
#' @export
get_system_info <- function() {
  
  info <- list(
    r_version = R.version.string,
    platform = R.version$platform,
    os = Sys.info()["sysname"],
    user = Sys.info()["user"],
    hostname = Sys.info()["nodename"],
    timestamp = as.character(Sys.time()),
    wd = getwd(),
    loaded_packages = (.packages())
  )
  
  return(info)
}

#' Print Session Information
#'
#' Prints detailed session information
#'
#' @param output_file Optional file to save session info
#' @export
print_session_info <- function(output_file = NULL) {
  
  info <- sessionInfo()
  
  if (!is.null(output_file)) {
    sink(output_file)
    print(info)
    sink()
    
    message(sprintf("Session info saved to: %s", output_file))
  } else {
    print(info)
  }
  
  invisible(info)
}
