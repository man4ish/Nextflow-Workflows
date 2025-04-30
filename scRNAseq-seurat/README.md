# Seurat-based Single Cell RNA-seq Pipeline (Nextflow)

## Overview

This pipeline performs single-cell RNA-seq analysis using:
- `Cell Ranger` for alignment and quantification
- `Seurat` for QC and filtering in R

## Requirements

- Docker or Singularity
- Nextflow
- nf-test (for testing)

## Usage

```bash
nextflow run main.nf --samplesheet samplesheet.csv --outdir results
```
