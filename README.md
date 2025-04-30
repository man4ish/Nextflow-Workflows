# 🧬 Monorepo: Nextflow Pipelines Collection

This repository contains multiple bioinformatics workflows implemented in **Nextflow**, each designed for a specific omics data type or computational purpose. It includes:

- 💠 Exome sequencing pipeline
- 🌍 Whole genome sequencing pipeline
- 🎧 RNA-seq pipeline
- 🧬 Methylation (Bisulfite) sequencing pipeline
- 🔬 Single-cell RNA-seq pipeline using **Seurat**
- ☁️ AWS Batch support
- 🧪 `nf-test` unit tests
- 🚀 GitHub Actions CI/CD
- 🔄 Kubernetes deployment example

---

## 📁 Directory Structure


```
├── exome_pipeline/
├── methylation_pipeline/
├── rnaseq_pipeline/
├── scRNAseq-seurat/
├── nextflow_awsbatch/
├── nextflow-monitoring-optimization/
├── nextflow-parallelization-example/
├── nextflow-pipeline-testing-demo/
├── nextflow-rna-seq-kubernetes/
├── .github/workflows/          # CI/CD for all pipelines
└── .gitmodules                 # Submodules (e.g., Kubernetes)

```

---

## 🚀 Getting Started

### 🔧 Prerequisites

- [Nextflow](https://www.nextflow.io/)
- Docker or Singularity
- Git
- [`nf-test`](https://github.com/nextflow-io/nf-test) (for testing pipelines)

### 🧪 Install nf-test

```bash
pip install nf-test
```
▶️ Running a Pipeline
Each pipeline is self-contained. Navigate to the desired pipeline folder and run it:

```
cd rnaseq_pipeline
nextflow run main.nf --samplesheet samplesheet.csv --outdir results
```
📘 All pipelines accept standard inputs like --samplesheet and output to --outdir.

🧪 Running Unit Tests
Tests are located in each pipeline’s test/ folder and can be executed using:

```
nf-test test test/
```
You can also test individual processes:

```
nf-test test test/test_star_align.nf
```
🔁 CI/CD with GitHub Actions
Each pipeline is automatically tested on push or PR via GitHub Actions (.github/workflows/ci.yml).

CI steps include:

Linting and validating pipeline structure

Running unit tests with nf-test

Reporting test results

# Pipeline Summary

## Pipeline	Description
- exome_pipeline/	Variant calling from WES data (BWA, GATK, etc.)
- rnaseq_pipeline/	Bulk RNA-seq analysis (STAR, Salmon, DESeq2)
- scRNAseq-seurat/	Single-cell RNA-seq analysis using Seurat (R)
- methylation_pipeline/	Bisulfite-seq using Bismark and methylKit
- nextflow_awsbatch/	AWS Batch compute environment for Nextflow
- nextflow-parallelization-example/	Example of parallel execution in Nextflow
- nextflow-monitoring-optimization/	Integration with monitoring tools like Tower
  
## Sample Commands
```
nextflow run exome_pipeline/main.nf --samplesheet exome_samples.csv --outdir results
nextflow run scRNAseq-seurat/main.nf --samplesheet cells.csv --outdir output/
```
## Contributing
- Fork this repository.

- Add or improve a pipeline inside its dedicated folder.

- Add unit tests using nf-test under test/.

- Update .github/workflows/ if needed.

- Submit a pull request!

