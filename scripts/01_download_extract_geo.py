from pathlib import Path
import tarfile
import urllib.request

# ============================================================
# Download and extract GSE144153 original GEO archive
#
# The RAW archive and extracted files are stored directly in:
#   ../data/
#
# No data/extracted/ subfolder is created.
# ============================================================

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT.parent / "data"

BASE_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144153/suppl/"
RAW_FILE = "GSE144153_RAW.tar"
RAW_URL = BASE_URL + RAW_FILE

DATA_DIR.mkdir(parents=True, exist_ok=True)

raw_tar_path = DATA_DIR / RAW_FILE

if not raw_tar_path.exists() or raw_tar_path.stat().st_size < 1_000_000:
    print("Downloading GEO dataset...")
    urllib.request.urlretrieve(RAW_URL, raw_tar_path)
else:
    print("GEO dataset already exists. Skipping download.")

print(f"RAW archive location: {raw_tar_path}")
print(f"RAW archive size MB: {raw_tar_path.stat().st_size / 1024 / 1024:.2f}")

print("Extracting archive directly into data/ ...")
with tarfile.open(raw_tar_path, "r") as tar:
    try:
        tar.extractall(DATA_DIR, filter="data")
    except TypeError:
        tar.extractall(DATA_DIR)

gene_files = sorted([p.name for p in DATA_DIR.glob("*pileupGENE*.txt.gz")])
te_files = sorted([p.name for p in DATA_DIR.glob("*pileupTE*.txt.gz")])

print(f"Extracted gene count files detected: {len(gene_files)}")
print(f"Extracted TE count files detected: {len(te_files)}")
print(f"Data directory: {DATA_DIR}")
print("Download and extraction completed.")
