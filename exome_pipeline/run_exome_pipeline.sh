#!/bin/sh

# Set variables
INPUT_DIR="data"
OUTPUT_DIR="results"
REFERENCE="ref/genome.fa"
DBSNP="/path/to/dbsnp.vcf.gz"
THREADS=8

# Validate input directory
if [ ! -d "$INPUT_DIR" ]; then
  echo "ERROR: Input directory '$INPUT_DIR' does not exist."
  exit 1
fi

# Validate reference genome
if [ ! -f "$REFERENCE" ]; then
  echo "ERROR: Reference genome file '$REFERENCE' not found."
  exit 1
fi

# Check if dbSNP file exists
if [ -f "$DBSNP" ]; then
  DBSNP_ARG="--dbSNP $DBSNP"
else
  echo "WARNING: dbSNP file '$DBSNP' not found. The BaseRecalibrator process may fail or be skipped."
  DBSNP_ARG=""
fi

# Run the Nextflow pipeline
nextflow run main.nf \
  --input "$INPUT_DIR" \
  --output "$OUTPUT_DIR" \
  --reference "$REFERENCE" \
  --fastqc_threads "$THREADS" \
  $DBSNP_ARG
