# R Environment Management with renv

This directory contains files for managing the R environment.

## Quick Start

1. **Install renv** (if not already installed):
   ```r
   install.packages("renv")
   ```

2. **Initialize renv in this project**:
   ```r
   renv::init()
   ```

3. **Restore packages from lockfile**:
   ```r
   renv::restore()
   ```

4. **After installing new packages, snapshot the state**:
   ```r
   renv::snapshot()
   ```

## Required Packages

Based on the analysis scripts, the following packages are required:

- **Seurat** (4.4.0) - Main single-cell analysis framework
- **scran** - Normalization
- **scater** - QC and visualization
- **scDblFinder** - Doublet detection
- **TSCAN** - Trajectory inference
- **monocle3** - Trajectory analysis
- **slingshot** - Lineage inference
- **tradeSeq** - Gene expression dynamics
- **clustree** - Cluster stability
- **ggplot2** - Visualization
- **pheatmap** - Heatmaps
- **DT** - Interactive tables
- **SeuratWrappers** - Seurat utilities
- **loupeR** - Loupe Browser integration

## Notes

- The `renv.lock` file will be generated when you first run `renv::init()`
- Make sure to run `renv::snapshot()` after installing any new packages
- Commit `renv.lock` to version control to ensure reproducibility
