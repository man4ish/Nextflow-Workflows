#!/usr/bin/env nextflow

// Define your input files or parameters
params.input_fastq = '/path/to/input.fastq'   // Path to input FastQ file (replace with actual path)
params.reference_genome = '/path/to/genome.fa'  // Path to reference genome (replace with actual path)

workflow {
  // Define the order in which processes should run
  fastqc(input_fastq)
  trimming(input_fastq)
  bwa(input_fastq, reference_genome)
  gatk(input_bam)
}
