# Singularity Guide for Seurat Pipeline

Complete guide for using your existing Docker Hub image with Singularity to run the Seurat pipeline.

---

## 🎯 Overview

**You already have:**
- ✅ Docker image on Docker Hub with all packages: `gynecoloji/bioinformatic_r_4_5_1:v2`
- ✅ Production Seurat pipeline code in `seurat_pipeline/`

**What you'll do:**
1. Pull your Docker image as a Singularity container (one-time setup)
2. Run the pipeline using that container
3. Works on both local machines and HPC clusters

**No Docker needed!** - Singularity can pull directly from Docker Hub

---

## 🚀 Quick Start (3 Steps)

### **Step 1: Build Singularity Container**

```bash
cd docker_R

# Run the build script
bash scripts/build_singularity.sh

# This creates: containers/bioinformatic_r_4.5.1.sif
```

⏱️ **Takes ~5 minutes** (downloads ~2.5 GB)

### **Step 2: Test with Example**

```bash
# Run example workflow
singularity exec \
  --bind $PWD:/workspace \
  containers/bioinformatic_r_4.5.1.sif \
  Rscript /workspace/seurat_pipeline/examples/example_workflow.R
```

### **Step 3: Run on Your Data**

```bash
# Using the helper script
bash scripts/run_singularity_local.sh \
  data/your_data.rds \
  results/your_analysis
```

**Done!** Results will be in `results/your_analysis/`

---

## 📂 Directory Structure

```
docker_R/
│
├── containers/                        # NEW - Singularity containers
│   └── bioinformatic_r_4.5.1.sif     # Your R container (pulled from Docker Hub)
│
├── scripts/                           # NEW - Helper scripts
│   ├── build_singularity.sh          # Build container
│   ├── run_singularity_local.sh      # Run locally
│   ├── run_singularity_hpc.sh        # Submit to HPC
│   └── run_singularity_interactive.sh # Interactive R
│
├── seurat_pipeline/                   # Your analysis code
│   ├── R/
│   ├── config/
│   └── run_seurat_pipeline.R
│
├── data/                              # Input data
├── results/                           # Output results
└── logs/                              # Log files
```

---

## 🛠️ Setup Instructions

### **One-Time Setup**

```bash
# 1. Navigate to your repository
cd docker_R

# 2. Create necessary directories
mkdir -p containers scripts data results logs

# 3. Place the helper scripts in scripts/
# (Download from outputs)

# 4. Make scripts executable
chmod +x scripts/*.sh

# 5. Build the Singularity container
bash scripts/build_singularity.sh
```

**What happens:**
- Pulls `gynecoloji/bioinformatic_r_4_5_1:v1` from Docker Hub
- Converts to Singularity format
- Saves as `containers/bioinformatic_r_4.5.1.sif`
- Tests that R and Seurat are working

---

## 💡 Usage Examples

### **Example 1: Run Complete Pipeline**

```bash
# Using the helper script (recommended)
bash scripts/run_singularity_local.sh \
  data/pbmc_data.rds \
  results/pbmc_analysis

# Or run directly with singularity
singularity exec \
  --bind $PWD:/workspace \
  containers/bioinformatic_r_4.5.1.sif \
  Rscript /workspace/seurat_pipeline/run_seurat_pipeline.R \
    --config /workspace/seurat_pipeline/config/default_config.yaml \
    --input /workspace/data/pbmc_data.rds \
    --output /workspace/results/pbmc_analysis \
    --cores 4
```

### **Example 2: Interactive R Session**

```bash
# Using the helper script
bash scripts/run_singularity_interactive.sh

# Inside R:
> source('/workspace/seurat_pipeline/R/utils.R')
> source('/workspace/seurat_pipeline/R/qc_functions.R')
> obj <- readRDS('/workspace/data/pbmc_data.rds')
> obj <- add_qc_metrics(obj, species='human')
> # ... continue working interactively
```

### **Example 3: Custom Configuration**

```bash
# 1. Edit config
nano seurat_pipeline/config/default_config.yaml

# 2. Run with custom settings
bash scripts/run_singularity_local.sh \
  data/my_data.rds \
  results/custom_analysis \
  seurat_pipeline/config/default_config.yaml \
  8  # number of cores
```

### **Example 4: HPC SLURM Job**

```bash
# 1. Edit the HPC script
nano scripts/run_singularity_hpc.sh

# Update these lines:
#   INPUT_FILE="data/my_dataset.rds"
#   OUTPUT_DIR="results/my_analysis_${SLURM_JOB_ID}"

# 2. Submit job
sbatch scripts/run_singularity_hpc.sh

# 3. Check status
squeue -u $USER

# 4. Check logs
tail -f logs/slurm_*.out
```

### **Example 5: Skip QC or Batch Correction**

```bash
# Skip QC filtering (data already filtered)
singularity exec \
  --bind $PWD:/workspace \
  containers/bioinformatic_r_4.5.1.sif \
  Rscript /workspace/seurat_pipeline/run_seurat_pipeline.R \
    --config /workspace/seurat_pipeline/config/default_config.yaml \
    --input /workspace/data/filtered_data.rds \
    --skip-qc

# Skip batch correction
singularity exec \
  --bind $PWD:/workspace \
  containers/bioinformatic_r_4.5.1.sif \
  Rscript /workspace/seurat_pipeline/run_seurat_pipeline.R \
    --config /workspace/seurat_pipeline/config/default_config.yaml \
    --input /workspace/data/my_data.rds \
    --skip-batch
```

---

## 🖥️ HPC-Specific Instructions

### **For Your Cluster (easley015)**

**1. Load Singularity Module:**
```bash
# Check available modules
module avail

# Load singularity (adjust name for your system)
module load singularity
# or
module load apptainer
```

**2. Build Container on Login Node:**
```bash
# On login node
cd /path/to/SnakeMake_RNAseq
bash scripts/build_singularity.sh
```

**3. Submit SLURM Job:**
```bash
# Edit job parameters in the script
nano scripts/run_singularity_hpc.sh

# Adjust these SBATCH parameters:
#SBATCH --time=24:00:00        # Max runtime
#SBATCH --cpus-per-task=20     # Number of CPUs
#SBATCH --mem=64G              # Memory needed
#SBATCH --partition=general    # Your partition name

# Submit
sbatch scripts/run_singularity_hpc.sh
```

**4. Monitor Job:**
```bash
# Check queue
squeue -u $USER

# Check output (while running)
tail -f logs/slurm_*.out

# After completion
cat logs/slurm_*.out
```

---

## 🔧 Understanding Singularity Commands

### **Basic Singularity Commands**

```bash
# Pull from Docker Hub
singularity pull my_container.sif docker://username/image:tag

# Execute a command
singularity exec container.sif command args

# Interactive shell
singularity shell container.sif

# Interactive R
singularity exec container.sif R

# Check what's inside
singularity inspect container.sif
```

### **Bind Mounts (Important!)**

Singularity **automatically mounts** your home directory, but you need to explicitly bind other directories:

```bash
# Bind current directory to /workspace inside container
singularity exec --bind $PWD:/workspace container.sif command

# Bind multiple directories
singularity exec \
  --bind /data:/workspace/data \
  --bind /results:/workspace/results \
  container.sif command

# Using environment variable
export SINGULARITY_BIND="/data:/workspace/data,/results:/workspace/results"
singularity exec container.sif command
```

**Path Translation:**
- **Outside container:** `/home/user/SnakeMake_RNAseq/data/file.rds`
- **Inside container:** `/workspace/data/file.rds`

---

## 📊 What You Get

After running the pipeline, you'll have:

```
results/
└── my_analysis/
    ├── qc/
    │   └── project_name/
    │       ├── QC_Cell_Density_UMI.pdf
    │       ├── QC_Cell_Density_genes.pdf
    │       ├── QC_Cell_correlations.pdf
    │       └── filtering_statistics.csv
    │
    ├── plots/
    │   ├── project_name_umap.pdf
    │   ├── project_name_elbow_plot.pdf
    │   └── project_name_cluster_qc.pdf
    │
    └── objects/
        ├── project_name_qc_filtered.rds
        ├── project_name_processed.rds
        ├── project_name_final.rds
        ├── project_name_summary.csv
        └── project_name_sessionInfo.txt

logs/
├── pipeline.log
└── slurm_*.out  # If running on HPC
```

---

## ⚡ Performance Tips

### **1. Use Local Scratch Space on HPC**

```bash
# Copy data to node's local storage
cp data/large_file.rds $TMPDIR/
INPUT_FILE="$TMPDIR/large_file.rds"

# Run pipeline
singularity exec ...

# Copy results back
cp -r $TMPDIR/results/* results/
```

### **2. Optimize CPU Usage**

```bash
# Use all available cores on compute node
CORES=$(nproc)

singularity exec ... \
  Rscript ... --cores $CORES
```

### **3. Pre-cache Container**

```bash
# Build once, use many times
# Store in a shared location on HPC
SIF_FILE="/shared/containers/bioinformatic_r_4.5.1.sif"
```

---

## 📋 Command Reference

### **Build Container**
```bash
bash scripts/build_singularity.sh
```

### **Run Pipeline (Local)**
```bash
bash scripts/run_singularity_local.sh INPUT OUTPUT [CONFIG] [CORES]
```

### **Run Pipeline (HPC)**
```bash
sbatch scripts/run_singularity_hpc.sh
```

### **Interactive R**
```bash
bash scripts/run_singularity_interactive.sh
```

### **Direct Singularity Execution**
```bash
singularity exec \
  --bind $PWD:/workspace \
  containers/bioinformatic_r_4.5.1.sif \
  Rscript /workspace/seurat_pipeline/run_seurat_pipeline.R \
    --config /workspace/seurat_pipeline/config/default_config.yaml \
    --input /workspace/data/file.rds
```

---

## ✅ Verification Checklist

- [ ] Container built successfully
- [ ] Example workflow runs without errors
- [ ] Results appear in `results/example/`
- [ ] Can run with your own data
- [ ] Interactive R session works
- [ ] HPC job submission works (if applicable)

---

## 🎯 Advantages of This Approach

✅ **Simple** - Use existing Docker image, no new image to build  
✅ **Fast** - One-time download, then instant execution  
✅ **Portable** - Same container works everywhere  
✅ **HPC-friendly** - Singularity is standard on HPC  
✅ **No sudo** - Works without admin privileges  
✅ **Reproducible** - Fixed environment every time  

---

## 🔄 Workflow Summary

```
1. Build Container (once)
   bash scripts/build_singularity.sh
   
2. Run Analysis (many times)
   bash scripts/run_singularity_local.sh data/file.rds results/output
   
3. Check Results
   ls results/output/
```

---

**Questions? Issues? Check the logs:**
- Pipeline execution: `logs/pipeline.log`
- HPC jobs: `logs/slurm_*.out`

**Ready to start!** 🚀
