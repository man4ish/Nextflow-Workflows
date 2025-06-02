#!/bin/bash -ue
mkdir -p out/aligned_bam
bwa mem -t 8 ./ref/genome.fa sample1_R1.fastq.gz null |         samtools view -bS |         samtools sort -o out/aligned_bam/sample1.bam
