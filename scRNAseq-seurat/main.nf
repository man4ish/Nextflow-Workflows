nextflow.enable.dsl=2

params.samplesheet = "samplesheet.csv"
params.outdir = "results"

workflow {
    samples = Channel.fromPath(params.samplesheet)

    cellranger_output = cellranger_count(samples)
    seurat_qc(cellranger_output)
}

process cellranger_count {
    input:
    path samplesheet

    output:
    path "cellranger_out"

    container 'biocontainers/cellranger:6.0.1--py39'

    script:
    """
    mkdir -p cellranger_out
    while IFS=',' read -r sample fastq_dir ref_path; do
        cellranger count --id=\$sample \\
                         --transcriptome=\$ref_path \\
                         --fastqs=\$fastq_dir \\
                         --sample=\$sample \\
                         --output-dir=cellranger_out/\$sample
    done < $samplesheet
    """
}

process seurat_qc {
    input:
    path cellranger_out

    output:
    path "${params.outdir}/seurat_qc_output"

    container 'quay.io/biocontainers/r-seurat:4.3.0--r42h03ef668_0'

    script:
    """
    mkdir -p ${params.outdir}/seurat_qc_output
    Rscript scripts/seurat_analysis.R $cellranger_out ${params.outdir}/seurat_qc_output
    """
}
