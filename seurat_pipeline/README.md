# Production-Ready Seurat Pipeline for Single-Cell RNA-seq Analysis

A professional, containerized pipeline for single-cell RNA-seq data analysis using Seurat. Designed with best practices for data engineering and bioinformatics production environments.

[![Docker](https://img.shields.io/badge/docker-ready-blue)](Dockerfile)
[![R Version](https://img.shields.io/badge/R-4.5.1-blue)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/license-MIT-green)](../LICENSE)

---

## 🎯 Features

### Production-Quality Code
- ✅ **Comprehensive error handling** and input validation
- ✅ **Modular design** with separate modules for QC, processing, and visualization
- ✅ **Configuration-driven** workflow using YAML
- ✅ **Extensive logging** with timestamps and severity levels
- ✅ **Reproducibility** - saves session info and configurations
- ✅ **Well-documented** with roxygen2 format

### Analysis Capabilities
- Quality control with customizable thresholds
- Multiple batch correction methods (Harmony, CCA, RPCA, SCTransform)
- Parallel processing support
- Comprehensive visualization suite
- Flexible clustering parameters

---

## 🚀 Quick Start

### Option 1: Using the Base R Container

```bash
# Pull the base container
singularity pull docker://gynecoloji/bioinformatic_r_4_5_1:v2

# Run pipeline
singularity exec bioinformatic_r_4_5_1_v2.sif \
  Rscript seurat_pipeline/run_seurat_pipeline.R \
  --config config/default_config.yaml
```

### Option 3: Local R Installation

```R
# Source the pipeline functions
source("R/utils.R")
source("R/qc_functions.R")
source("R/processing_functions.R")
source("R/plotting_functions.R")

# Load configuration
config <- load_config("config/default_config.yaml")

# Run pipeline steps manually (see examples/)
```

---

## ⚙️ Configuration

The pipeline is controlled via YAML configuration files. See `config/default_config.yaml` for all available options.

### Key Configuration Sections:

**Quality Control**
```yaml
qc:
  min_umi: 500
  min_genes: 300
  max_mito: 25
  min_complexity: 0.8
  min_cells_per_gene: 10
```

**Processing**
```yaml
processing:
  normalization_method: "LogNormalize"
  n_features: 3000
  dims_use: 30
  resolution: 0.5
  n_cores: 4
```

**Batch Correction**
```yaml
batch_correction:
  enabled: true
  batch_var: "orig.ident"
  method: "harmony"  # Options: harmony, CCA, RPCA, SCTransform
  dims: 30
```

---

## 📊 Usage Examples

### Example 1: Basic Analysis

```bash
Rscript run_seurat_pipeline.R \
  --config config/default_config.yaml \
  --input data/pbmc_data.rds \
  --output results/pbmc_analysis
```

### Example 2: With Batch Correction

Create custom config `config/batch_correction.yaml`:
```yaml
batch_correction:
  enabled: true
  batch_var: "sample"
  method: "harmony"
  dims: 30
```

Run pipeline:
```bash
Rscript run_seurat_pipeline.R \
  --config config/batch_correction.yaml \
  --input data/multi_sample.rds
```

### Example 3: Skip QC (Pre-filtered Data)

```bash
Rscript run_seurat_pipeline.R \
  --config config/default_config.yaml \
  --input data/filtered_data.rds \
  --skip-qc
```

### Example 4: High-Performance Computing

```bash
# Use more cores for faster processing
Rscript run_seurat_pipeline.R \
  --config config/default_config.yaml \
  --cores 20 \
  --input data/large_dataset.rds
```

---

## 📦 Input Data Format

The pipeline accepts Seurat objects saved as `.rds` files:

```R
# Create and save a Seurat object
library(Seurat)

# From 10X data
data <- Read10X(data.dir = "path/to/10x/data")
seurat_obj <- CreateSeuratObject(counts = data, project = "my_project")

# Save for pipeline
saveRDS(seurat_obj, "data/my_input.rds")
```

Alternative inputs (counts matrix, etc.) can be easily accommodated by modifying the data loading section.

---

## 📈 Output Structure

```
results/
├── qc/                          # QC plots and statistics
│   ├── project_name/
│   │   ├── QC_Cell_Density_UMI.pdf
│   │   ├── QC_Cell_Density_genes.pdf
│   │   ├── filtering_statistics.csv
│   │   └── ...
├── plots/                       # Analysis visualizations
│   ├── project_name_umap.pdf
│   ├── project_name_elbow_plot.pdf
│   ├── project_name_cluster_qc.pdf
│   └── ...
├── objects/                     # Saved Seurat objects
│   ├── project_name_qc_filtered.rds
│   ├── project_name_processed.rds
│   ├── project_name_final.rds
│   ├── project_name_summary.csv
│   ├── project_name_sessionInfo.txt
│   └── project_name_config.yaml
└── logs/
    └── pipeline.log             # Detailed execution log
```

---

## 🧪 Running Tests

```bash
# Run unit tests (coming soon)
Rscript tests/run_tests.R

# Test with example data
Rscript examples/example_workflow.R
```

---

## 📚 Function Reference

### QC Functions

- `add_qc_metrics()` - Calculate QC metrics
- `plot_qc_metrics()` - Visualize QC metrics
- `filter_cells_and_genes()` - Filter based on QC thresholds

### Processing Functions

- `process_seurat_object()` - Complete processing pipeline
- `process_batch_correction()` - Batch correction workflow
- `get_default_processing_config()` - Get default parameters

### Plotting Functions

- `plot_batch_correction_results()` - Batch correction diagnostics
- `plot_gene_expression()` - Feature plots
- `plot_elbow()` - PCA elbow plot
- `plot_cluster_qc()` - Cluster quality metrics

### Utility Functions

- `check_and_load_packages()` - Package management
- `log_message()` - Logging
- `save_seurat_object()` / `load_seurat_object()` - I/O operations
- `load_config()` / `save_config()` - Configuration management

---

## 🤝 Contributing

This is a portfolio project demonstrating production-quality bioinformatics code. Suggestions and improvements are welcome!

---

## 📄 License

MIT License - see [LICENSE](../LICENSE) file for details.

---

## 👤 Author

**Ji (gynecoloji)**
- GitHub: [@gynecoloji](https://github.com/gynecoloji)
- Email: gynecoloji@gmail.com

---

## 🙏 Acknowledgments

- Built on the [Seurat](https://satijalab.org/seurat/) package
- Extends the `bioinformatic_r_4_5_1` container
- Batch correction methods: Harmony, Seurat integration tools

---

## 📝 Citation

If you use this pipeline in your research, please cite:

```
Ji. (2025). Production-Ready Seurat Pipeline for Single-Cell RNA-seq Analysis.
GitHub: https://github.com/gynecoloji/docker_R
```

---

**Last Updated**: January 2026  
**Version**: 1.0.0
