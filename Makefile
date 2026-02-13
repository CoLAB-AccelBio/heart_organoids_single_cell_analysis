.PHONY: help setup preprocess cluster de trajectory visualize all clean

help:
	@echo "Heart Organoids Single-Cell Analysis Pipeline"
	@echo ""
	@echo "Available targets:"
	@echo "  setup       - Set up R environment with renv"
	@echo "  preprocess  - Run preprocessing (QC, filtering, normalization)"
	@echo "  cluster     - Run clustering analysis"
	@echo "  de          - Run differential expression analysis"
	@echo "  trajectory  - Run trajectory inference"
	@echo "  visualize   - Generate figures"
	@echo "  all         - Run full pipeline"
	@echo "  clean       - Remove generated results"

# Set up R environment
setup:
	@echo "Setting up R environment..."
	Rscript -e "if (!require('renv')) install.packages('renv'); renv::init()"

# Run preprocessing
preprocess:
	@echo "Running preprocessing..."
	Rscript scripts/01_preprocessing/preprocessing.R

# Run clustering
cluster:
	@echo "Running clustering..."
	Rscript scripts/02_clustering/clustering.R

# Run differential expression
de:
	@echo "Running differential expression analysis..."
	Rscript scripts/03_differential_expression/DEA_heatmap_volcano_plots.R

# Run trajectory inference
trajectory:
	@echo "Running trajectory inference..."
	Rscript scripts/04_trajectory/trajectory_inference.R

# Generate visualizations
visualize:
	@echo "Generating figures..."
	Rscript scripts/05_visualization/Figures.R

# Run full pipeline
all: preprocess cluster de trajectory visualize

# Clean results
clean:
	@echo "Cleaning results..."
	rm -rf results/figures/*
	rm -rf results/tables/*
	rm -rf results/reports/*
	rm -rf data/processed/*.rds
