from pathlib import Path
import pandas as pd
import numpy as np

# ============================================================
# Exploratory fold-change summary for GSE144153
#
# This script does not perform formal statistical differential
# expression analysis. It generates reusable exploratory
# transcriptomic fold-change summaries from CPM-filtered data.
# ============================================================

REPO_ROOT = Path(__file__).resolve().parents[1]
PROCESSED_DIR = REPO_ROOT / "processed"

cpm_filtered = pd.read_csv(
    PROCESSED_DIR / "cpm_filtered_matrix.csv",
    index_col=0
)

sh587_cols = [c for c in cpm_filtered.columns if "sh587" in c]
shE_cols = [c for c in cpm_filtered.columns if "shE" in c]

mean_sh587 = cpm_filtered[sh587_cols].mean(axis=1) + 0.5
mean_shE = cpm_filtered[shE_cols].mean(axis=1) + 0.5

log2fc = np.log2(mean_sh587 / mean_shE)
base_mean = ((mean_sh587 - 0.5) + (mean_shE - 0.5)) / 2

fold_summary = pd.DataFrame({
    "Gene_Ensembl_ID": cpm_filtered.index,
    "BaseMean_CPM": base_mean.values,
    "Mean_CPM_sh587": mean_sh587.values - 0.5,
    "Mean_CPM_shE": mean_shE.values - 0.5,
    "Log2FoldChange": log2fc.values,
    "Abs_Log2FoldChange": np.abs(log2fc.values)
}).sort_values("Log2FoldChange", ascending=False)

fold_summary.to_csv(
    PROCESSED_DIR / "preliminary_fold_change_summary.csv",
    index=False
)

# Backward-compatible file name for users who may expect it
fold_summary.to_csv(
    PROCESSED_DIR / "preliminary_differential_expression.csv",
    index=False
)

volcano_table = fold_summary.copy()
volcano_table["log10_BaseMean_CPM_plus1"] = np.log10(
    volcano_table["BaseMean_CPM"] + 1
)

volcano_table["Regulation_class"] = "Stable_or_low_change"
volcano_table.loc[
    volcano_table["Log2FoldChange"] >= 1.5,
    "Regulation_class"
] = "Higher_in_sh587"

volcano_table.loc[
    volcano_table["Log2FoldChange"] <= -1.5,
    "Regulation_class"
] = "Lower_in_sh587"

volcano_table.to_csv(
    PROCESSED_DIR / "volcano_plot_table.csv",
    index=False
)

print("Exploratory fold-change summary completed.")
print(f"Fold-change summary saved to: {PROCESSED_DIR / 'preliminary_fold_change_summary.csv'}")
print(f"Backward-compatible copy saved to: {PROCESSED_DIR / 'preliminary_differential_expression.csv'}")
print(f"Volcano table saved to: {PROCESSED_DIR / 'volcano_plot_table.csv'}")
