# Nextflow Pipeline for Exome Analysis

This repository contains a Nextflow pipeline for bioinformatics analysis, including quality control (FastQC), alignment (BWA), duplicate marking (Picard), base recalibration (GATK), and variant calling (GATK). The pipeline supports containerized execution via Docker for easy reproducibility and deployment.

---

## Requirements

- [Nextflow](https://www.nextflow.io/)
- Docker (for containerized execution)
- Git (for version control)
- (Optional) Nextflow Tower account for easy pipeline management and monitoring

---

## Setup

## 1. Clone the repository

```bash
git clone https://github.com/your-username/nextflow-bioinformatics-pipeline.git
cd nextflow-bioinformatics-pipeline
```
## 2. Build Docker image
The Dockerfile includes all necessary bioinformatics tools (bwa, samtools, picard, gatk, fastqc).

```bash
docker build -t exome_pipeline .
```
Optionally, push the Docker image to a container registry (required for Nextflow Tower):

```bash
docker tag exome_pipeline your_dockerhub_username/exome_pipeline:latest
docker push your_dockerhub_username/exome_pipeline:latest
```
## 3. Configure input files
The pipeline requires input and reference files organized as below:

### Input files
- input/fastq/: Paired-end FASTQ files, e.g. sample1_R1.fastq.gz, sample1_R2.fastq.gz

### Reference files
- input/reference/:

- assembly38.fasta: Reference genome FASTA file (e.g., Homo_sapiens_assembly38.fasta)

- Homo_sapiens_assembly38.dbsnp138.vcf: dbSNP VCF file for recalibration

Example directory structure:

```bash

project/
├── main.nf
├── data/
│   └── sample1_R1.fastq.gz
│   └── sample1_R2.fastq.gz
├── ref/
│   └── genome.fa
│   └── dbsnp.vcf
```

## 4. Running the Pipeline
### Option 1: Run locally with Nextflow (Docker recommended)
```
bash
nextflow run main.nf \
  --input './data' \
  --reference './ref/genome.fa' \
  --dbSNP './ref/dbsnp.vcf' \
  -with-docker your_dockerhub_username/exome_pipeline:latest

Or without Docker (if all tools installed locally):

```bash
nextflow run main.nf \
  --input './data' \
  --reference './ref/genome.fa' \
  --dbSNP './ref/dbsnp.vcf'
```
### Option 2: Run with Docker directly

```bash
docker run --rm -v $(pwd):/mnt exome_pipeline run /mnt/main.nf \
  --input '/mnt/data' \
  --reference '/mnt/ref/genome.fa' \
  --dbSNP '/mnt/ref/dbsnp.vcf'
```

### Option 3: Run with Nextflow Tower
Nextflow Tower allows you to launch, monitor, and manage your pipelines easily via a web UI.

### Prerequisites
Nextflow Tower account: https://tower.nf/

Push your code to a Git repository (GitHub, GitLab, Bitbucket, etc.)

Push your Docker image to a container registry (e.g., Docker Hub)

Steps
Push Docker image to registry

```bash
docker tag nextflow-analysis your_dockerhub_username/nextflow-analysis:latest
docker push your_dockerhub_username/nextflow-analysis:latest
```
### Push your pipeline code

```bash
git add .
git commit -m "Add pipeline with Docker support for Tower"
git push origin main
```

## Configure pipeline in Tower

- Log in to Nextflow Tower

- Create a new pipeline and connect your Git repo URL

- Set pipeline file to main.nf and branch (e.g., main)

- Create or select an execution profile

- Enable Docker usage and specify your Docker image: your_dockerhub_username/nextflow-analysis:latest

- Set pipeline parameters (input paths, reference files) as needed

- Launch

- Choose your compute environment (local, cloud, Kubernetes, etc.)

- Launch the pipeline and monitor it via Tower UI

## 5. Output
The pipeline generates output in the following structure:

```bash

output/
├── fastqc/                # FastQC reports
├── aligned_bam/           # Sorted BAM files from alignment
├── bamfiles/
│   └── renamed_bam/       # BAM files after duplicate marking and base recalibration
└── gvcf/                  # GVCF files from HaplotypeCaller
```
## 6. Pipeline Parameters (can be adjusted in main.nf or Tower UI)
```
params.input = './data'
params.output = './output'
params.reference = './ref/genome.fa'
params.dbSNP = './ref/dbsnp.vcf'
```
Make sure these paths match your data locations.

## Notes
When running on cloud or HPC with Tower, ensure that input and reference files are accessible from the compute environment (e.g., shared file systems, cloud storage buckets).

The pipeline uses 8 threads by default for multi-threaded tools; adjust based on your compute resources.

