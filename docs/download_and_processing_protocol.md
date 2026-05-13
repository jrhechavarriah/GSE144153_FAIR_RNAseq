# Download and Processing Protocol

## Original data source

NCBI GEO accession:
GSE144153

Original study attribution:

Turelli et al., *Science Advances* (2020).  
DOI: 10.1126/sciadv.aba3200

## Local organization

The repository stores derived resources in:

```text
GSE144153_FAIR_RNAseq/
```

The original GEO archive and extracted count files are downloaded outside the repository into:

```text
../data/
```

This design avoids redistributing original raw data in the repository ZIP.

## Transcriptomic workflow

The workflow includes:

1. GEO download
2. Direct extraction of GEO supplementary files into `../data/`
3. pileupGENE file detection
4. raw count matrix construction
5. CPM normalization
6. low-expression filtering
7. PCA generation
8. library-size quality-control summary
9. top-variable-gene heatmap generation
10. hierarchical clustering dendrogram generation
11. exploratory fold-change summary generation
12. volcano-style exploratory visualization
13. FAIR metadata and documentation generation
14. repository validation
15. ZIP packaging

## Generated derived resources

The workflow generates:

- raw count matrix
- CPM-normalized matrix
- filtered CPM matrix
- PCA coordinates
- library-size QC table
- exploratory fold-change summary table
- volcano plot table
- heatmap source matrices
- sample label mapping table
- gene label mapping table
- figure manifest
- workflow and QC figures

## Important note

This repository distributes only derived computational resources generated from publicly available transcriptomic data.

The original sequencing repositories remain attributable to the original GEO dataset, the original source article, and associated public repositories.
