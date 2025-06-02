#!/bin/bash -ue
gatk --java-options '-Xmx4g' BaseRecalibrator         -R ref/genome.fa         -I sample1.dup.bam         --known-sites false         -O results/bam/sample1/sample1_recal_data.csv
