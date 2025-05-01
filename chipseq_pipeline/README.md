# ChIP-Seq Analysis Pipeline (Nextflow + Docker)

This is a modular ChIP-Seq analysis pipeline built using [Nextflow](https://www.nextflow.io/) and Docker. It performs raw FASTQ input processing through quality control, trimming, alignment, peak calling, and motif analysis.

## Tools Used

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) – Quality control
- [Cutadapt](https://cutadapt.readthedocs.io/) – Adapter trimming
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/) – Read alignment
- [MACS2](https://github.com/macs3-project/MACS) – Peak calling
- [HOMER](http://homer.ucsd.edu/homer/) – Motif discovery

## Setup

### 1. Build Docker Image

```bash
docker build -t chipseq-pipeline .
```

Make sure Docker is installed on your system.

### 2. Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
mv nextflow /usr/local/bin/
```

## Directory Structure

```
project/
├── data/
│   ├── sample1_R1.fastq
│   └── control.bam
├── main.nf
├── nextflow.config
├── Dockerfile
└── README.md
```

## Running the Pipeline

```bash
nextflow run main.nf -c nextflow.config
```

Optional parameters can be passed as follows:

```bash
nextflow run main.nf -c nextflow.config \
  --reads "./data/*.fastq" \
  --control "./data/control.bam" \
  --reference "/path/to/bowtie2/index/prefix"
```

## Output

Pipeline outputs are stored in `./results/`, including:

- `fastqc_reports/` – FastQC HTML reports
- `trimmed/` – Trimmed FASTQ files
- `aligned/` – SAM alignment files
- `peaks/` – MACS2 peak outputs
- `motif_analysis/` – HOMER motif findings

## Notes

- Ensure your Bowtie2 index is built and accessible.
- HOMER requires genome annotations. Ensure `findMotifsGenome.pl` is set up correctly.
- The pipeline is set for local execution. It can be modified for cluster/cloud environments in `nextflow.config`.

## Example Test Run

```bash
nextflow run main.nf \
  --reads "./data/*.fastq" \
  --control "./data/control.bam" \
  --reference "/genomes/hg38/bowtie2_index"
```

## Contact

For questions, improvements, or issues, please open an issue or pull request.

