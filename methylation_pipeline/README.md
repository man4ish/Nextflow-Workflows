# Methylation Analysis Pipeline

This is a **Nextflow** pipeline designed for **methylation analysis** using popular bioinformatics tools. The pipeline integrates **FastQC** for quality control, **Bismark** for read alignment, and **MethylKit** in R for methylation analysis. The pipeline runs each process inside Docker containers to ensure a consistent and reproducible environment.

## Requirements

- **Nextflow**: A workflow management system used to orchestrate the pipeline.
- **Docker**: A platform for developing, shipping, and running applications in containers. Each step in the pipeline runs within its own Docker container.
- **R** (with **MethylKit**): For methylation analysis in the `methylkit_analysis` process.

## Getting Started

### 1. Install Dependencies

Before running the pipeline, make sure you have **Nextflow** and **Docker** installed.

#### Install Nextflow:

Follow the installation instructions for your system:  
[Nextflow Installation Guide](https://www.nextflow.io/docs/latest/getstarted.html)

```bash
curl -s https://get.nextflow.io | bash
```
### Then build the Docker image using:

```bash
docker build -t my_methylation_pipeline .
```

### Run the Pipeline
```bash
nextflow run main.nf -process.executor docker
```

