#!/bin/bash -ue
ls -l genome.fa || echo "Reference file not found: genome.fa"
mkdir -p results/bam/sample1
bwa mem -t 8 genome.fa sample1_R1.fastq.gz sample1_R2.fastq.gz |     samtools view -bS - > results/bam/sample1/sample1.bam
