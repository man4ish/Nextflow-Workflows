# Nextflow Pipeline with AWS Batch

This is a Nextflow pipeline designed to run bioinformatics processes using AWS Batch. The pipeline includes multiple steps such as **FastQC**, **Trimmomatic**, **BWA**, and **GATK**, each running in its own Docker container.

## Requirements

- **Nextflow**: You need to have Nextflow installed on your system. You can install it via the following command:

  ```bash
  curl -s https://get.nextflow.io | bash
  ```

- AWS Account: An AWS account with the required AWS Batch resources configured (job queue, compute environment, IAM roles, etc.).

- Docker: Docker must be installed to use containerized processes. Make sure you have Docker installed and running on your system.

## Pipeline Overview
The pipeline consists of four bioinformatics steps, each executed in a Docker container:

- FastQC: Quality control for high throughput sequencing data.

- Trimmomatic: A tool for trimming Illumina sequencing data.

- BWA: A fast and accurate alignment program for mapping sequencing reads to a reference genome.

- GATK: The Genome Analysis Toolkit for variant discovery.

## The pipeline is configured to run on AWS Batch, using Docker containers for each process.

### Files
- nextflow.config: Configuration file for the Nextflow pipeline, defining AWS Batch settings, job parameters, and containers for each process.

- pipeline.nf: The main Nextflow script that defines the workflow steps and process dependencies.

### README.md: This file.

### Setup
- Install Nextflow: Make sure Nextflow is installed. If not, follow the instructions above.

- AWS Batch Setup: You must have an AWS Batch environment set up with the following:

- A job queue

### A compute environment

- IAM roles with appropriate permissions

- A shared S3 bucket for storing input/output files (if needed)

- Configure AWS CLI: Ensure that the AWS CLI is configured on your system with appropriate credentials and region settings.

### Edit Configuration: Update nextflow.config with your specific AWS Batch settings:

- Replace your-job-queue-name with the name of your AWS Batch job queue.

- Replace your-compute-environment-name, your-iam-role-name, and your-job-role-name with your AWS Batch environment and IAM role names.

- Set the correct AWS region and S3 work directory.

### Input Files: 
Prepare your input files, such as FASTQ files and reference genome. Update the mypipeline.nf file to point to the correct paths for these files.

### Running the Pipeline
To run the pipeline, use the following command:

``` bash

nextflow run mypipeline.nf -profile awsbatch
```
This will execute the pipeline using AWS Batch with the configuration set in nextflow.config.

## Monitoring with Nextflow Tower
### 1. Create a Nextflow Tower Account:
Sign up for Nextflow Tower if you don’t already have an account. Go to Nextflow Tower and create a new account.

### 2. Connect Nextflow to Tower:
Once you have a Nextflow Tower account, you will need to connect your local Nextflow setup to Tower. Use the following command to log in to Tower from your terminal:

```
nextflow tower login
```

This will open a web browser window asking you to authenticate with your Nextflow Tower account. Once logged in, it will give you a token to enter in your terminal to authenticate the connection.

### 3. Run the Pipeline with Tower Integration:
To submit your pipeline to AWS Batch and track it using Nextflow Tower, run the pipeline with the -with-tower flag:

```
nextflow run mypipeline.nf -profile awsbatch -with-tower
```

This command will execute the pipeline and automatically send execution logs, status, and results to your Nextflow Tower account. You will be able to monitor the status of each task in the UI.

### 4. Monitor Execution on Tower:
- After running the pipeline with Tower integration, go to your Nextflow Tower dashboard.

- You will see the pipeline listed there with detailed information on task execution, logs, errors, and progress.

- The UI allows you to monitor individual task performance and manage the execution remotely if needed.
