#!/bin/bash

set -e

echo "=========================================="
echo "Heart Organoids Single-Cell Analysis"
echo "=========================================="

# Default values
STEPS="all"
SAMPLES="VEGF,PDGFBB"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --steps)
            STEPS="$2"
            shift 2
            ;;
        --samples)
            SAMPLES="$2"
            shift 2
            ;;
        --help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --steps STR    Pipeline steps to run (default: all)"
            echo "                 Options: preprocess, cluster, de, trajectory, visualize, all"
            echo "  --samples STR  Comma-separated sample names (default: VEGF,PDGFBB)"
            echo "  --help         Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Activate conda environment
eval "$(conda shell.bash activate)"
conda activate heart_organoids_sc

# Create results directories
mkdir -p results/figures results/tables results/reports

# Run pipeline steps
run_step() {
    local step=$1
    echo ""
    echo ">>> Running: $step"
    case $step in
        preprocess)
            Rscript scripts/01_preprocessing/preprocessing.R --samples "$SAMPLES"
            ;;
        cluster)
            Rscript scripts/02_clustering/clustering.R --samples "$SAMPLES"
            ;;
        de)
            Rscript scripts/03_differential_expression/DEA_heatmap_volcano_plots.R --samples "$SAMPLES"
            ;;
        trajectory)
            Rscript scripts/04_trajectory/trajectory_inference.R --samples "$SAMPLES"
            ;;
        visualize)
            Rscript scripts/05_visualization/Figures.R --samples "$SAMPLES"
            ;;
        all)
            run_step "preprocess"
            run_step "cluster"
            run_step "de"
            run_step "trajectory"
            run_step "visualize"
            ;;
        *)
            echo "Unknown step: $step"
            exit 1
            ;;
    esac
}

if [[ "$STEPS" == "all" ]]; then
    run_step "all"
else
    IFS=',' read -ra STEP_ARRAY <<< "$STEPS"
    for s in "${STEP_ARRAY[@]}"; do
        run_step "$s"
    done
fi

echo ""
echo "=========================================="
echo "Pipeline completed successfully!"
echo "=========================================="
