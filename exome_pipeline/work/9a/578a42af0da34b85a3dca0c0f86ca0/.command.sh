#!/bin/bash -ue
fastqc --nogroup -t 8 sample1_R1.fastq.gz sample1_R2.fastq.gz
