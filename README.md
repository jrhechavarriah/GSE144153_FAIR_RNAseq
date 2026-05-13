# Curated RNA-seq Expression Dataset and Reproducible Python Workflow Derived from GSE144153 for FAIR Transcriptomic Analysis of ZNF587-Silenced Human H1 Cells

## Description

This repository contains a curated and reproducible derived transcriptomic dataset generated from publicly available RNA-seq data associated with ZNF587-silenced human H1 cells from GEO accession GSE144153.

The repository distributes only derived computational resources generated from publicly available sequencing count data and does not redistribute the original raw sequencing repositories.

## Original source dataset

NCBI GEO accession: GSE144153

Original study attribution:

Turelli et al., *Science Advances* (2020).  
DOI: 10.1126/sciadv.aba3200

## Repository structure

- `data/`
  - External folder generated outside the repository during execution.
  - Stores the downloaded GEO RAW archive and extracted files.
  - This folder is not included in the repository ZIP.

- `processed/`
  - Raw count matrix
  - CPM-normalized matrix
  - Filtered CPM matrix
  - PCA coordinates
  - Exploratory fold-change summary
  - Volcano plot table
  - Gene and sample label mappings
  - Figure manifest

- `figures/`
  - Workflow diagram
  - PCA visualization
  - Library-size QC barplot
  - Heatmap
  - Hierarchical clustering dendrograms
  - Exploratory volcano-style plot

- `scripts/`
  - Reproducible computational workflows
  - Validation script
  - ZIP packaging script

- `docs/`
  - Data dictionary
  - Dataset reuse notes
  - Download and processing protocol
  - Environment information

- `notebooks/`
  - Interactive Colab/Jupyter reproducibility notebook

## Reproducibility

Recommended execution:

```bash
bash run_all.sh
```

Alternative modular execution:

```bash
python scripts/01_download_extract_geo.py
python scripts/02_build_expression_matrix.py
python scripts/03_normalization_and_qc.py
python scripts/04_differential_expression.py
python scripts/05_generate_figures.py
python scripts/06_validate_repository.py
python scripts/07_package_zip.py
```

## Label conventions used in figures

Compact labels are used to improve readability and maintain publication-quality visualizations.

### Sample labels

| Label | GEO Sample ID | Condition |
|---|---|---|
| S1 | GSM4282198_sh587 | sh587 |
| S2 | GSM4282199_shE | shE |
| S3 | GSM4282200_shE | shE |
| S4 | GSM4282201_sh587 | sh587 |
| S5 | GSM4282202_shE | shE |
| S6 | GSM4282203_sh587 | sh587 |

The complete mapping is available in:

```text
metadata/sample_label_mapping.csv
```

### Variable gene labels

The top 50 most variable genes were relabeled as VG1–VG50.

The complete mapping between VG labels and Ensembl identifiers is available in:

```text
metadata/gene_label_mapping.csv
```

## Reuse potential

This dataset may support:

- transcriptomic exploratory analysis
- RNA-seq normalization benchmarking
- PCA and dimensionality reduction studies
- transcriptomic clustering analysis
- heatmap and dendrogram visualization benchmarking
- FAIR bioinformatics workflows
- transcriptomic machine learning preprocessing pipelines
- reproducible computational genomics education

## Important methodological note

The repository provides exploratory fold-change summaries and visualization resources. It does not claim formal statistical differential expression inference because the workflow does not implement DESeq2, edgeR, limma-voom, p-value estimation, or multiple-testing correction.

## Citation

If you use this dataset, please cite the associated repository and publication.

## Authors

Jesus Rafael Hechavarria-Hernandez a,*  
Cristian Vacacela Gomez b,*

## Affiliations

a Universidad Ecotec, Research Department, Samborondón, Guayas, EC092302, Ecuador.

b University of Calabria, Department of Environmental Engineering (DIAm), Rende, 87036, Italy

## Corresponding authors

jhechavarria@ecotec.edu.ec  
cristianisaac.vacacelagomez@fis.unical.it

## ORCID

J.R. Hechavarria-Hernandez: 0000-0002-9013-8665

C. Vacacela Gomez: 0000-0002-9248-9944

## License

The derived computational products are distributed under CC BY 4.0.

The original sequencing data, biological samples, and associated primary transcriptomic resources remain attributable to the original study and public repositories.
