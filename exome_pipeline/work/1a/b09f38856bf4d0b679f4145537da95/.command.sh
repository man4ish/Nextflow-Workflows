#!/bin/bash -ue
mkdir -p results/bam/sample1
mkdir -p results/dup_log
java -jar /opt/picard.jar MarkDuplicates         INPUT=sample1.bam         OUTPUT=results/bam/sample1/sample1.dup.bam         METRICS_FILE=results/dup_log/sample1_metrics.txt         VALIDATION_STRINGENCY=LENIENT         REMOVE_DUPLICATES=false         ASSUME_SORTED=false
