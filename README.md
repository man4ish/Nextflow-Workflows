# Nextflow Pipeline for Bioinformatics Analysis

This repository contains a Nextflow pipeline for bioinformatics analysis, including quality control (FastQC), alignment (BWA), duplicate marking (Picard), base recalibration (GATK), and variant calling (GATK). The pipeline is encapsulated within a Docker container for easy reproducibility and deployment.

## Requirements

- [Nextflow](https://www.nextflow.io/)
- Docker (for containerized execution)
- Git (for version control)

# Setup

## 1. Clone the repository

```bash
git clone https://github.com/your-username/nextflow-bioinformatics-pipeline.git
cd nextflow-bioinformatics-pipeline
```

## 2. Build Docker image
The Dockerfile is provided for building the container with all necessary tools (bwa, samtools, picard, gatk, fastqc).

```bash
docker build -t nextflow-analysis .
```
## 3. Configure the input files
The pipeline requires several input files, which should be placed in specific directories.

### Input files:
input/fastq/: Directory containing paired-end FASTQ files for your samples (e.g., sample1_R1.fastq.gz, sample1_R2.fastq.gz).

input/reference/: Directory containing reference files for alignment and variant calling:

assembly38.fasta: Reference genome FASTA file (e.g., Homo_sapiens_assembly38.fasta).

Homo_sapiens_assembly38.dbsnp138.vcf: dbSNP VCF file for variant recalibration.

Example directory structure:
```bash
input/
    fastq/
        sample1_R1.fastq.gz
        sample1_R2.fastq.gz
    reference/
        assembly38.fasta
        Homo_sapiens_assembly38.dbsnp138.vcf
```
## 4. Running the Pipeline
You can run the pipeline either with Docker or directly on your machine with Nextflow.

Run the pipeline using Docker:
```bash
docker run --rm -v $(pwd):/mnt nextflow-analysis run /mnt/main.nf
```
Run the pipeline with Nextflow (if installed locally):
```bash
nextflow run main.nf
```
5. Output
The pipeline will generate output in the following directory:

```bash
output/
    fastqc/
    aligned_bam/
    bamfiles/
        renamed_bam/
fastqc/: FastQC results



aligned_bam/: Aligned BAM files

bamfiles/renamed_bam/: Processed BAM files with duplicates marked and base recalibrated

out/g.vcf/: GVCF files from HaplotypeCaller
```

## Configuration
The pipeline's input and reference files are configured via the following parameters in the main.nf file:

```bash
params.input = './in/non_redundant_fastq/'
params.output = './out'
params.reference = './in/reference/assembly38.fasta'
params.dbSNP = './in/reference/Homo_sapiens_assembly38.dbsnp138.vcf'
Make sure to adjust these paths if your file locations differ.
```
