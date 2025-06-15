#!/bin/bash

# Create directories
mkdir -p data docs figures results tests notebooks scripts

# Create empty or placeholder files
touch CHANGELOG.md
touch LICENSE
touch README.md
touch environment.yml
touch notebooks/plot_vaf_metrics.ipynb

# Scripts
touch scripts/downsample.sh
touch scripts/extract_flagstat.sh
touch scripts/extract_idxstats.sh
touch scripts/compare_metrics.py
touch scripts/validate_pipeline.sh

# Confirm creation
echo "Project structure initialized successfully:"
tree -L 2
