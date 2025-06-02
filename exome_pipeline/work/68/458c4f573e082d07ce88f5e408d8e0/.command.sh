#!/bin/bash -ue
# Create the output directory for this sample's FastQC results
mkdir -p results/fastqc/sample1

# Run FastQC
# --nogroup: Prevents grouping of identical sequences, which can be useful.
# -t 8: Use the configured number of threads.
# sample1_R1.fastq.gz sample1_R2.fastq.gz: Access the individual read files from the 'reads' list.
# --outdir: Specify the output directory for FastQC reports.
fastqc --nogroup -t 8 sample1_R1.fastq.gz sample1_R2.fastq.gz --outdir results/fastqc/sample1
