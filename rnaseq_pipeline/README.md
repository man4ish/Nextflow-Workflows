# RNA-Seq Nextflow Pipeline
This Nextflow RNA-Seq pipeline performs a complete RNA-Seq analysis, including adapter trimming, transcript quantification using Kallisto, read alignment with STAR, BAM sorting and indexing, and RNA-Seq QC metrics generation using Picard. This pipeline is fully containerized using Docker to ensure reproducibility and portability.

## Features
- Adapter Trimming with Skewer

- Transcript Quantification with Kallisto

- Read Alignment with STAR

- BAM Sorting & Indexing using Samtools

- RNA-Seq QC Summary with Picard

- Containerized execution using Docker

- Scalable and portable workflow with Nextflow

## Requirements
Before running this pipeline, ensure the following are installed:

- Nextflow: For running and managing the workflow. Installation instructions.

- Docker: For containerized execution of all tasks. Installation instructions.

- Java: Required by Nextflow for execution.

## Setup
Clone this repository or download the rnaseq.nf Nextflow script.

Ensure you have the Docker image: docker.io/man4ish/rnaseq:latest.

## Inputs
The pipeline requires the following inputs:

- Input Name	Description
- fastq1	Forward FASTQ file (paired-end)
- fastq2	Reverse FASTQ file (paired-end)
- adapters	Adapter sequence file for Skewer
- skewer_threads	Number of threads for trimming
- minimum_read_length	Minimum read length to retain after trimming
- prefix	Prefix for output files
- kallisto_threads	Number of threads for Kallisto
- bootstrap_samples	Number of bootstrap samples for Kallisto
- idx	Kallisto index file
- gtf	GTF annotation file for transcriptome
- STAR_threads	Number of threads for STAR alignment
- ref_tar	Compressed reference genome bundle (.tar.gz) for STAR
- ref_flat	RefFlat annotation file for Picard
- ribosomal_interval	Ribosomal intervals file for Picard
- ref_seq	Reference FASTA file for Picard
  
## Outputs
The pipeline will produce the following outputs:

- Trimmed FASTQ files

- Kallisto quantification results (tar archive)

- STAR-aligned BAM file (sorted & indexed)

- RNA-Seq summary metrics and coverage plots from Picard

## Docker Container
All tasks are executed within a Docker container to ensure a consistent environment. The required Docker image is:

```
docker.io/man4ish/rnaseq:latest
```

## Example Execution
1. Create an inputs.config file
Create a configuration file that defines all the input parameters for the pipeline. Here's an example inputs.config:

```
params.fastq1 = "sample_R1.fastq.gz"
params.fastq2 = "sample_R2.fastq.gz"
params.adapters = "adapters.fa"
params.skewer_threads = 4
params.minimum_read_length = 50
params.prefix = "sample1"

params.kallisto_threads = 4
params.bootstrap_samples = 100
params.idx = "transcripts.idx"
params.gtf = "annotation.gtf"

params.STAR_threads = 8
params.ref_tar = "star_reference.tar.gz"

params.ref_flat = "reference.refFlat"
params.ribosomal_interval = "rRNA.interval_list"
params.ref_seq = "reference.fasta"
```

2. Run the pipeline with Nextflow
Execute the pipeline using Nextflow with the following command:
```
nextflow run rnaseq.nf -c inputs.config
```

This will run the pipeline with the configuration defined in the inputs.config file.


## Workflow Details

- The main Nextflow process includes several tasks such as:

- trim: Trims adapters from the raw FASTQ files using Skewer.

- quantification: Quantifies transcript abundances using Kallisto.

- align: Aligns reads using STAR.

- sort_index: Sorts and indexes the BAM file using Samtools.

- gen_summary: Generates RNA-Seq metrics and plots using Picard.