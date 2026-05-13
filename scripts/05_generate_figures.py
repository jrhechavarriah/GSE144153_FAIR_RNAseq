from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist

# ============================================================
# Final publication-oriented figure generation script for
# GSE144153_FAIR_RNAseq
#
# Uses consistent figure labels:
# Samples: S1–S6
# Genes: VG1–VG50
#
# Generates:
# 1. workflow_diagram.png
# 2. pca_qc.png
# 3. library_size_qc_barplot.png
# 4. top50_heatmap_zscore.png
# 5. top50_sample_gene_dendrograms.png
# 6. volcano_plot.png
# 7. figure_manifest.csv
# 8. sample_label_mapping.csv
# 9. gene_label_mapping.csv
# ============================================================

REPO_ROOT = Path(__file__).resolve().parents[1]
PROCESSED_DIR = REPO_ROOT / "processed"
FIG_DIR = REPO_ROOT / "figures"
META_DIR = REPO_ROOT / "metadata"

FIG_DIR.mkdir(parents=True, exist_ok=True)
META_DIR.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
def short_sample_labels(sample_names):
    ordered = list(sample_names)
    mapping = {}
    for i, sample in enumerate(ordered, start=1):
        condition = "sh587" if "sh587" in sample else "shE"
        mapping[sample] = {
            "Sample_Label": f"S{i}",
            "Original_Sample_ID": sample,
            "Condition": condition
        }
    return mapping

def apply_sample_labels(df, mapping):
    return df.rename(columns={k: v["Sample_Label"] for k, v in mapping.items()})

# ------------------------------------------------------------
# Load core matrices and tables
# ------------------------------------------------------------
cpm_file = PROCESSED_DIR / "cpm_filtered_matrix.csv"

if not cpm_file.exists():
    raise FileNotFoundError(
        f"Expected file not found: {cpm_file}\n"
        "Run scripts/03_normalization_and_qc.py first."
    )

cpm_filtered = pd.read_csv(cpm_file, index_col=0)

sample_mapping = short_sample_labels(cpm_filtered.columns)

sample_mapping_df = pd.DataFrame(sample_mapping.values())
sample_mapping_df.to_csv(META_DIR / "sample_label_mapping.csv", index=False)
sample_mapping_df.to_csv(PROCESSED_DIR / "sample_label_mapping.csv", index=False)

sample_label_dict = {
    original: item["Sample_Label"]
    for original, item in sample_mapping.items()
}

sample_condition_dict = {
    item["Sample_Label"]: item["Condition"]
    for item in sample_mapping.values()
}

cpm_labeled = apply_sample_labels(cpm_filtered, sample_mapping)

# ------------------------------------------------------------
# 1. Vertical workflow diagram
# ------------------------------------------------------------
workflow_steps = [
    ("1", "GEO source", "GSE144153 public RNA-seq dataset"),
    ("2", "Download", "Retrieve GSE144153_RAW.tar from NCBI GEO"),
    ("3", "Extraction", "Extract sample-level pileupGENE count files"),
    ("4", "Matrix construction", "Build raw gene count matrix"),
    ("5", "Normalization", "Apply CPM normalization"),
    ("6", "Filtering", "Retain genes with CPM > 1 in at least 3 samples"),
    ("7", "Quality control", "Generate PCA and library-size summaries"),
    ("8", "Exploratory analysis", "Generate heatmap, dendrograms, and volcano-style summaries"),
    ("9", "FAIR package", "Export reusable tables, figures, metadata, and scripts")
]

fig, ax = plt.subplots(figsize=(8.5, 11))
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")

y_positions = np.linspace(0.92, 0.10, len(workflow_steps))

for i, ((num, title, desc), y) in enumerate(zip(workflow_steps, y_positions)):
    ax.text(
        0.12, y, num,
        ha="center", va="center",
        fontsize=12, fontweight="bold",
        bbox=dict(boxstyle="circle,pad=0.35", linewidth=1.2, facecolor="white")
    )
    ax.text(
        0.22, y + 0.018, title,
        ha="left", va="center",
        fontsize=11, fontweight="bold"
    )
    ax.text(
        0.22, y - 0.020, desc,
        ha="left", va="center",
        fontsize=9
    )

    if i < len(workflow_steps) - 1:
        ax.annotate(
            "",
            xy=(0.12, y_positions[i + 1] + 0.035),
            xytext=(0.12, y - 0.035),
            arrowprops=dict(arrowstyle="->", lw=1.3)
        )

ax.set_title(
    "Reproducible workflow for the GSE144153 curated transcriptomic dataset",
    fontsize=13,
    fontweight="bold",
    pad=18
)

plt.tight_layout()
plt.savefig(FIG_DIR / "workflow_diagram.png", dpi=300, bbox_inches="tight")
plt.close()

# ------------------------------------------------------------
# 2. PCA plot with S-labels
# ------------------------------------------------------------
pca_file = PROCESSED_DIR / "pca_coordinates.csv"

if not pca_file.exists():
    raise FileNotFoundError(
        f"Expected file not found: {pca_file}\n"
        "Run scripts/03_normalization_and_qc.py first."
    )

pca = pd.read_csv(pca_file)

pca["Sample_Label"] = pca["sample_id"].map(sample_label_dict)
pca["Condition"] = pca["Sample_Label"].map(sample_condition_dict)

plt.figure(figsize=(8, 6))

for group in ["sh587", "shE"]:
    subset = pca[pca["Condition"] == group]
    plt.scatter(
        subset["PC1"],
        subset["PC2"],
        label=group,
        s=90
    )

for _, row in pca.iterrows():
    plt.text(
        row["PC1"] + 0.6,
        row["PC2"] + 0.4,
        row["Sample_Label"],
        ha="left",
        va="bottom",
        fontsize=10
    )

if "PC1_variance_percent" in pca.columns and "PC2_variance_percent" in pca.columns:
    pc1 = pca["PC1_variance_percent"].iloc[0]
    pc2 = pca["PC2_variance_percent"].iloc[0]
    plt.xlabel(f"PC1 ({pc1:.2f}%)")
    plt.ylabel(f"PC2 ({pc2:.2f}%)")
else:
    plt.xlabel("PC1")
    plt.ylabel("PC2")

plt.title("PCA of GSE144153 RNA-seq samples")
plt.legend(
    title="Condition")#,
#    fontsize=12,
#    title_fontsize=13,
#    labelspacing=1.6,
#    borderpad=1.0,
#)
plt.tight_layout()
plt.savefig(FIG_DIR / "pca_qc.png", dpi=300, bbox_inches="tight")
plt.close()

# ------------------------------------------------------------
# 3. Library-size QC barplot with S-labels
# ------------------------------------------------------------
library_qc_file = PROCESSED_DIR / "library_size_qc.csv"

if not library_qc_file.exists():
    raise FileNotFoundError(
        f"Expected file not found: {library_qc_file}\n"
        "Run scripts/03_normalization_and_qc.py first."
    )

library_qc = pd.read_csv(library_qc_file)

library_qc["Sample_Label"] = library_qc["sample_id"].map(sample_label_dict)
library_qc["Condition"] = library_qc["Sample_Label"].map(sample_condition_dict)

plt.figure(figsize=(8, 5))

for condition in ["sh587", "shE"]:
    subset = library_qc[library_qc["Condition"] == condition]
    plt.bar(
        subset["Sample_Label"],
        subset["library_size_raw_counts"],
        label=condition
    )

plt.ylabel("Total raw gene counts")
plt.xlabel("Samples")
plt.title("Library-size QC across GSE144153 samples")
plt.legend(title="Condition")
plt.tight_layout()
plt.savefig(FIG_DIR / "library_size_qc_barplot.png", dpi=300, bbox_inches="tight")
plt.close()

# ------------------------------------------------------------
# 4. Prepare top 50 variable genes matrices
# ------------------------------------------------------------
log_cpm = np.log2(cpm_labeled + 1)

top_50_genes = (
    log_cpm.var(axis=1)
    .sort_values(ascending=False)
    .head(50)
    .index
)

heatmap_log2cpm = log_cpm.loc[top_50_genes]

heatmap_zscore = heatmap_log2cpm.sub(
    heatmap_log2cpm.mean(axis=1),
    axis=0
).div(
    heatmap_log2cpm.std(axis=1).replace(0, np.nan),
    axis=0
).fillna(0)

gene_label_mapping = pd.DataFrame({
    "Gene_Label": [f"VG{i+1}" for i in range(len(heatmap_zscore.index))],
    "Gene_Ensembl_ID": heatmap_zscore.index,
    "Mean_log2_CPM_plus1": heatmap_log2cpm.mean(axis=1).values,
    "Variance_log2_CPM_plus1": heatmap_log2cpm.var(axis=1).values
})

gene_label_mapping.to_csv(PROCESSED_DIR / "gene_label_mapping.csv", index=False)
gene_label_mapping.to_csv(META_DIR / "gene_label_mapping.csv", index=False)

gene_label_dict = dict(
    zip(gene_label_mapping["Gene_Ensembl_ID"], gene_label_mapping["Gene_Label"])
)

heatmap_zscore_labeled = heatmap_zscore.copy()
heatmap_zscore_labeled.index = [
    gene_label_dict[g] for g in heatmap_zscore.index
]

heatmap_log2cpm.to_csv(PROCESSED_DIR / "top50_variable_genes_log2cpm.csv")
heatmap_zscore.to_csv(PROCESSED_DIR / "top50_variable_genes_zscore.csv")
heatmap_zscore_labeled.to_csv(PROCESSED_DIR / "top50_variable_genes_zscore_labeled.csv")

# ------------------------------------------------------------
# 5. Simple Z-score heatmap with S-labels and VG labels
# ------------------------------------------------------------
plt.figure(figsize=(10, 14))

sns.heatmap(
    heatmap_zscore_labeled,
    cmap="coolwarm",
    center=0,
    vmin=-2.5,
    vmax=2.5,
    linewidths=0.0,
    cbar_kws={"label": "Row Z-score"}
)

plt.title(
    "Top 50 Most Variable Genes - Row Z-score",
    fontsize=14,
    fontweight="bold"
)
plt.xlabel("Samples")
plt.ylabel("Variable Gene Label")
plt.xticks(rotation=0)
plt.yticks(fontsize=8)
plt.tight_layout()
plt.savefig(FIG_DIR / "top50_heatmap_zscore.png", dpi=300, bbox_inches="tight")
plt.close()

# Remove old combined clustered heatmap if present
old_combined = FIG_DIR / "top50_clustered_heatmap_zscore.png"
if old_combined.exists():
    old_combined.unlink()

# ------------------------------------------------------------
# 6. Sample and gene dendrograms with S/VG labels
# ------------------------------------------------------------
sample_linkage = linkage(
    pdist(heatmap_zscore_labeled.T),
    method="average"
)

gene_linkage = linkage(
    pdist(heatmap_zscore_labeled),
    method="average"
)

fig = plt.figure(figsize=(12, 9))

ax1 = fig.add_axes([0.08, 0.60, 0.84, 0.30])
dendrogram(
    sample_linkage,
    labels=heatmap_zscore_labeled.columns.tolist(),
    leaf_rotation=0,
    leaf_font_size=10,
    ax=ax1
)
ax1.set_title(
    "Hierarchical clustering of RNA-seq samples",
    fontsize=13,
    fontweight="bold"
)
ax1.set_ylabel("Distance")

ax2 = fig.add_axes([0.08, 0.08, 0.84, 0.36])
dendrogram(
    gene_linkage,
    labels=heatmap_zscore_labeled.index.tolist(),
    leaf_rotation=90,
    leaf_font_size=7,
    ax=ax2
)
ax2.set_title(
    "Hierarchical clustering of the top 50 most variable genes",
    fontsize=13,
    fontweight="bold"
)
ax2.set_ylabel("Distance")
ax2.set_xlabel("Variable Gene Label")

plt.savefig(FIG_DIR / "top50_sample_gene_dendrograms.png", dpi=300, bbox_inches="tight")
plt.close()

# ------------------------------------------------------------
# 7. Volcano plot
# ------------------------------------------------------------
volcano_file = PROCESSED_DIR / "volcano_plot_table.csv"

if not volcano_file.exists():
    raise FileNotFoundError(
        f"Expected file not found: {volcano_file}\n"
        "Run scripts/04_differential_expression.py first."
    )

volcano = pd.read_csv(volcano_file)

plt.figure(figsize=(9, 6))

class_order = [
    "Stable_or_low_change",
    "Higher_in_sh587",
    "Lower_in_sh587"
]

for cls in class_order:
    subset = volcano[volcano["Regulation_class"] == cls]
    plt.scatter(
        subset["Log2FoldChange"],
        subset["log10_BaseMean_CPM_plus1"],
        s=12,
        alpha=0.65,
        label=cls
    )

plt.axvline(1.5, linestyle="--", linewidth=1)
plt.axvline(-1.5, linestyle="--", linewidth=1)

plt.xlabel("Log2 fold change (sh587 / shE)")
plt.ylabel("log10(BaseMean CPM + 1)")
plt.title("Exploratory volcano plot for GSE144153")
plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig(FIG_DIR / "volcano_plot.png", dpi=300, bbox_inches="tight")
plt.close()

# ------------------------------------------------------------
# 8. Figure manifest table
# ------------------------------------------------------------
figure_manifest = pd.DataFrame([
    {
        "figure_file": "workflow_diagram.png",
        "description": "Vertical workflow diagram summarizing the reproducible data-processing pipeline."
    },
    {
        "figure_file": "pca_qc.png",
        "description": "PCA quality-control plot using compact sample labels S1–S6."
    },
    {
        "figure_file": "library_size_qc_barplot.png",
        "description": "Library-size QC barplot using compact sample labels S1–S6."
    },
    {
        "figure_file": "top50_heatmap_zscore.png",
        "description": "Heatmap of the top 50 most variable genes using row-wise Z-score values, VG1–VG50 labels, and S1–S6 sample labels."
    },
    {
        "figure_file": "top50_sample_gene_dendrograms.png",
        "description": "Dendrogram figure showing hierarchical clustering of samples and the top 50 most variable genes using S and VG labels."
    },
    {
        "figure_file": "volcano_plot.png",
        "description": "Exploratory volcano-style plot based on log2 fold change and average CPM expression."
    }
])

figure_manifest.to_csv(PROCESSED_DIR / "figure_manifest.csv", index=False)

print("Figures generated.")
print(f"Workflow figure saved to: {FIG_DIR / 'workflow_diagram.png'}")
print(f"PCA figure saved to: {FIG_DIR / 'pca_qc.png'}")
print(f"QC barplot saved to: {FIG_DIR / 'library_size_qc_barplot.png'}")
print(f"Simple heatmap saved to: {FIG_DIR / 'top50_heatmap_zscore.png'}")
print(f"Sample/gene dendrograms saved to: {FIG_DIR / 'top50_sample_gene_dendrograms.png'}")
print(f"Volcano plot saved to: {FIG_DIR / 'volcano_plot.png'}")
print(f"Sample label mapping saved to: {META_DIR / 'sample_label_mapping.csv'}")
print(f"Gene label mapping saved to: {META_DIR / 'gene_label_mapping.csv'}")
print(f"Figure manifest saved to: {PROCESSED_DIR / 'figure_manifest.csv'}")
