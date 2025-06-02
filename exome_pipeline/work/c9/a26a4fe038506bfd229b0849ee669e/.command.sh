#!/bin/bash -ue
mkdir -p out/fastqc/sample1
fastqc -t 8 sample1_R1.fastq.gz --outdir out/fastqc/sample1
