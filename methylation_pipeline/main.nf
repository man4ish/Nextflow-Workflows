#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define the process for quality control with FastQC
process fastqc {
    input:
    path reads

    output:
    path "fastqc_report"

    // Specify the Docker image for FastQC
    container 'quay.io/biocontainers/fastqc:v0.11.9--0'

    script:
    """
    fastqc $reads --outdir=fastqc_report
    """
}

// Define the process for aligning reads with Bismark
process bismark_align {
    input:
    path fastq_files

    output:
    path "aligned_bam"

    // Specify the Docker image for Bismark
    container 'quay.io/biocontainers/bismark:v0.23.0--h7e3b6e0_0'

    script:
    """
    bismark --genome /path/to/genome --output_dir ./ $fastq_files
    """
}

// Define the process for methylation analysis with MethylKit in R
process methylkit_analysis {
    input:
    path bam_files

    output:
    path "methylation_results.csv"

    // Specify the Docker image for MethylKit in R
    container 'quay.io/biocontainers/r-methylkit:v1.9.0--r40hfc51f62_0'

    script:
    """
    Rscript methylkit_analysis.R $bam_files
    """
}

// Define the workflow
workflow {
    // Define the input FASTQ files (you would replace this with your actual data path)
    reads = file('input/*.fastq')

    // Run the processes in the pipeline
    fastqc(reads)
    bismark_align(reads)
    methylkit_analysis(bam_files)
}
