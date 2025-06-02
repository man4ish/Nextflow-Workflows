#!/bin/bash -ue
bwa mem -t 8 genome.fa sample1_R1.fastq.gz sample1_R2.fastq.gz |     samtools view -bS - > sample1.bam
