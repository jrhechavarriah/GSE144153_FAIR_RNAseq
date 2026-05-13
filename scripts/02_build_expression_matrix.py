from pathlib import Path
import pandas as pd

# ============================================================
# Build expression matrix from extracted GEO pileupGENE files
#
# Input:
#   ../data/*pileupGENE*.txt.gz
#
# Outputs:
#   metadata/sample_annotation.csv
#   processed/raw_counts_matrix.csv
# ============================================================

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT.parent / "data"

PROCESSED_DIR = REPO_ROOT / "processed"
META_DIR = REPO_ROOT / "metadata"

PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
META_DIR.mkdir(parents=True, exist_ok=True)

gene_files = sorted([p.name for p in DATA_DIR.glob("*pileupGENE*.txt.gz")])

if not gene_files:
    raise FileNotFoundError(
        f"No pileupGENE files found in {DATA_DIR}. "
        "Run scripts/01_download_extract_geo.py first."
    )

sample_records = []

for f in gene_files:
    gsm = f.split("_")[0]
    condition = "sh587" if "sh587" in f else "shE"
    sample_id = f"{gsm}_{condition}"

    sample_records.append({
        "sample_id": sample_id,
        "gsm_id": gsm,
        "condition": condition,
        "file_name": f
    })

sample_annotation = pd.DataFrame(sample_records)
sample_annotation.to_csv(META_DIR / "sample_annotation.csv", index=False)

count_tables = []

for _, row in sample_annotation.iterrows():
    file_path = DATA_DIR / row["file_name"]

    temp = pd.read_csv(
        file_path,
        sep="\t",
        header=None,
        names=["Gene_Ensembl_ID", row["sample_id"]],
        compression="gzip"
    )

    temp = temp.drop_duplicates(subset="Gene_Ensembl_ID")
    temp = temp.set_index("Gene_Ensembl_ID")

    count_tables.append(temp)

raw_counts = pd.concat(count_tables, axis=1).fillna(0).astype(int)
raw_counts.to_csv(PROCESSED_DIR / "raw_counts_matrix.csv")

print(f"Sample annotation saved to: {META_DIR / 'sample_annotation.csv'}")
print(f"Raw count matrix saved to: {PROCESSED_DIR / 'raw_counts_matrix.csv'}")
print(f"Raw count matrix shape: {raw_counts.shape}")
