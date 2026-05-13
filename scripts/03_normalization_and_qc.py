from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# ============================================================
# Normalization and quality-control script for
# GSE144153_FAIR_RNAseq
#
# Outputs:
# - processed/library_size_qc.csv
# - processed/cpm_normalized_matrix.csv
# - processed/cpm_filtered_matrix.csv
# - processed/pca_coordinates.csv
#
# Note:
# This script intentionally does not generate processed/dataset_metadata.csv.
# Dataset-level metadata is maintained only at the repository root as:
# - dataset_metadata.csv
# ============================================================

REPO_ROOT = Path(__file__).resolve().parents[1]
PROCESSED_DIR = REPO_ROOT / "processed"
META_DIR = REPO_ROOT / "metadata"

raw_counts = pd.read_csv(PROCESSED_DIR / "raw_counts_matrix.csv", index_col=0)
sample_annotation = pd.read_csv(META_DIR / "sample_annotation.csv")

# ------------------------------------------------------------
# Library-size QC
# ------------------------------------------------------------
library_sizes = raw_counts.sum(axis=0)

library_qc = pd.DataFrame({
    "sample_id": library_sizes.index,
    "library_size_raw_counts": library_sizes.values
})

library_qc = library_qc.merge(
    sample_annotation[["sample_id", "condition"]],
    on="sample_id",
    how="left"
)

library_qc.to_csv(PROCESSED_DIR / "library_size_qc.csv", index=False)

# ------------------------------------------------------------
# CPM normalization
# ------------------------------------------------------------
cpm = raw_counts.divide(library_sizes, axis=1) * 1_000_000
cpm.to_csv(PROCESSED_DIR / "cpm_normalized_matrix.csv")

# ------------------------------------------------------------
# Low-expression filtering
# ------------------------------------------------------------
cpm_filtered = cpm[(cpm > 1).sum(axis=1) >= 3]
cpm_filtered.to_csv(PROCESSED_DIR / "cpm_filtered_matrix.csv")

# ------------------------------------------------------------
# PCA from top 1000 most variable genes
# ------------------------------------------------------------
log_cpm = np.log2(cpm_filtered + 1)

top_var_genes = (
    log_cpm.var(axis=1)
    .sort_values(ascending=False)
    .head(1000)
    .index
)

pca_input = log_cpm.loc[top_var_genes].T

pca = PCA(n_components=2)
pca_values = pca.fit_transform(
    StandardScaler().fit_transform(pca_input)
)

pca_df = pd.DataFrame({
    "PC1": pca_values[:, 0],
    "PC2": pca_values[:, 1],
    "sample_id": pca_input.index,
    "condition": ["sh587" if "sh587" in s else "shE" for s in pca_input.index],
    "PC1_variance_percent": pca.explained_variance_ratio_[0] * 100,
    "PC2_variance_percent": pca.explained_variance_ratio_[1] * 100
})

pca_df.to_csv(PROCESSED_DIR / "pca_coordinates.csv", index=False)

print("Normalization completed.")
print(f"Library-size QC table saved to: {PROCESSED_DIR / 'library_size_qc.csv'}")
print(f"CPM matrix saved to: {PROCESSED_DIR / 'cpm_normalized_matrix.csv'}")
print(f"Filtered CPM matrix saved to: {PROCESSED_DIR / 'cpm_filtered_matrix.csv'}")
print(f"PCA coordinates saved to: {PROCESSED_DIR / 'pca_coordinates.csv'}")
print("Dataset-level metadata is maintained only at repository root: dataset_metadata.csv")
