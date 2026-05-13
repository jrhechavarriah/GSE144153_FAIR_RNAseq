#!/usr/bin/env python3
"""
Create the final ZIP package for GSE144153_FAIR_RNAseq.

This script excludes the external ../data/ folder, cache files, checkpoints,
temporary files, and previous ZIP packages. If executed in Google Colab,
it triggers automatic download.

Usage:
    python scripts/07_package_zip.py
"""

from pathlib import Path
import zipfile


REPO_ROOT = Path(__file__).resolve().parents[1]
ZIP_NAME = f"{REPO_ROOT.name}.zip"
ZIP_PATH = REPO_ROOT.parent / ZIP_NAME


def should_exclude(path: Path) -> bool:
    excluded_dirs = {
        "__pycache__",
        ".ipynb_checkpoints",
        ".git",
        ".github",
        "data"
    }

    excluded_files = {
        ".DS_Store",
        "Thumbs.db"
    }

    excluded_suffixes = {
        ".pyc",
        ".pyo",
        ".tmp",
        ".log"
    }

    if any(part in excluded_dirs for part in path.parts):
        return True

    if path.name in excluded_files:
        return True

    if path.suffix in excluded_suffixes:
        return True

    if path.name.startswith("~$"):
        return True

    if path.name == ZIP_NAME:
        return True

    return False


def create_zip() -> Path:
    if ZIP_PATH.exists():
        ZIP_PATH.unlink()

    included_files = 0

    with zipfile.ZipFile(ZIP_PATH, "w", zipfile.ZIP_DEFLATED) as z:
        for file_path in REPO_ROOT.rglob("*"):
            if file_path.is_file() and not should_exclude(file_path):
                arcname = file_path.relative_to(REPO_ROOT.parent)
                z.write(file_path, arcname)
                included_files += 1

    print("[done] ZIP package created:")
    print(ZIP_PATH)
    print(f"[info] Included files: {included_files}")
    print(f"[info] ZIP size MB: {ZIP_PATH.stat().st_size / 1024 / 1024:.2f}")

    return ZIP_PATH


def download_if_colab(zip_path: Path):
    try:
        from google.colab import files
        print("[colab] Starting automatic download...")
        files.download(str(zip_path))
    except Exception as e:
        print("[info] Automatic download unavailable.")
        print("[info] ZIP file location:")
        print(zip_path)
        print("[debug]", e)


def main():
    zip_path = create_zip()
    download_if_colab(zip_path)


if __name__ == "__main__":
    main()
