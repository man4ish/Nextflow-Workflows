# 🧬 Whole Genome Sequencing Pipeline

This is a modular and reproducible **Nextflow** pipeline for **whole genome sequencing (WGS)** analysis. It includes the following steps: read trimming, quality control, alignment, sorting, deduplication, variant calling, and annotation.

---

## 🚀 Features

- **Trimmomatic**: Adapter trimming and quality filtering
- **FastQC**: Quality check of raw and trimmed reads
- **BWA**: Read alignment to reference genome
- **Samtools & Picard**: Sorting, indexing, and marking duplicates
- **GATK**: Variant calling using HaplotypeCaller
- **SnpEff**: Variant annotation
- **Dockerized**: Fully containerized environment
- **CI/CD**: GitHub Actions integrated
- **Unit Testing**: Process-level testing using [`nf-test`](https://github.com/nextflow-io/nf-test)

---

## 🏗️ Project Structure

```text
.
├── main.nf                 # Nextflow pipeline
├── nextflow.config         # Pipeline configuration
├── Dockerfile              # Docker image with all tools
├── .github/workflows/ci.yml # CI/CD pipeline for testing
├── tests/
│   ├── data/               # Minimal test data
│   ├── trimmomatic.nf.test
│   ├── fastqc.nf.test
│   ├── bwa_align.nf.test
│   ├── sort_index.nf.test
│   ├── mark_duplicates.nf.test
│   ├── call_variants.nf.test
│   └── annotate_snps.nf.test
```

🐳 Docker
Build and tag the Docker image:
```
docker build -t wgs-pipeline:latest .
```
Use the image in the pipeline by adding this to your nextflow.config:

```
process.container = 'wgs-pipeline:latest'
```

🔧 Run the Pipeline
```
bash
nextflow run main.nf \
  --reads 'data/*_R{1,2}.fastq.gz' \
  --reference 'genome.fa' \
  --snpeff_db 'GRCh38.86' \
  --outdir 'results'
```

## Run all tests:

```bash
nf-test tests/
```

## Run a specific test:

```bash
nf-test tests/trimmomatic.nf.test
```
