#!/bin/bash -ue
ls -l /ref/genome.fa* || echo "Reference files not found in /ref"
mkdir -p results/bam/sample1
bwa mem -t 8 /ref/genome.fa sample1_R1.fastq.gz sample1_R2.fastq.gz |     samtools view -bS - > results/bam/sample1/sample1.bam
