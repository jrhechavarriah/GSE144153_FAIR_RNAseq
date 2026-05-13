#!/bin/bash
set -e

echo "========================================"
echo "Running GSE144153 reproducibility pipeline"
echo "========================================"

python scripts/01_download_extract_geo.py
python scripts/02_build_expression_matrix.py
python scripts/03_normalization_and_qc.py
python scripts/04_differential_expression.py
python scripts/05_generate_figures.py
python scripts/06_validate_repository.py

echo "========================================"
echo "Pipeline completed successfully"
echo "========================================"
