from pathlib import Path

# ============================================================
# Repository validation script for
# GSE144153_FAIR_RNAseq
#
# This script validates the expected derived resources after
# running:
#   bash run_all.sh
# ============================================================

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT.parent / "data"

required_files = [
    # Processed matrices and tables
    "processed/raw_counts_matrix.csv",
    "processed/cpm_normalized_matrix.csv",
    "processed/cpm_filtered_matrix.csv",
    "processed/preliminary_fold_change_summary.csv",
    "processed/volcano_plot_table.csv",
    "processed/pca_coordinates.csv",
    "processed/library_size_qc.csv",
    "processed/top50_variable_genes_log2cpm.csv",
    "processed/top50_variable_genes_zscore.csv",
    "processed/top50_variable_genes_zscore_labeled.csv",
    "processed/figure_manifest.csv",

    # Figures
    "figures/workflow_diagram.png",
    "figures/pca_qc.png",
    "figures/library_size_qc_barplot.png",
    "figures/top50_heatmap_zscore.png",
    "figures/top50_sample_gene_dendrograms.png",
    "figures/volcano_plot.png",

    # Metadata
    "metadata/sample_annotation.csv",
    "metadata/sample_label_mapping.csv",
    "metadata/gene_label_mapping.csv",

    # Root-level documentation
    "README.md",
    "LICENSE",
    "requirements.txt",
    "CITATION.cff",
    "dataset_metadata.csv",

    # Docs
    "docs/data_dictionary.csv",
    "docs/dataset_reuse_notes.md",
    "docs/download_and_processing_protocol.md",
    "docs/environment_info.txt",
]

required_external_data = [
    DATA_DIR / "GSE144153_RAW.tar"
]

missing = []

for f in required_files:
    path = REPO_ROOT / f
    if not path.exists():
        missing.append(str(f))

for path in required_external_data:
    if not path.exists():
        missing.append(str(path))

if len(missing) == 0:
    print("Repository validation successful.")
    print(f"Repository root: {REPO_ROOT}")
    print(f"External data directory: {DATA_DIR}")
else:
    print("Missing files:")
    for m in missing:
        print(m)
    raise FileNotFoundError("Repository validation failed.")
