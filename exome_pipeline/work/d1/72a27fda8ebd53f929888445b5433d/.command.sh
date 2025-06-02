#!/bin/bash -ue
mkdir -p results/fastqc/sample1
fastqc --nogroup -t 8 sample1_R1.fastq.gz sample1_R2.fastq.gz --outdir results/fastqc/sample1
