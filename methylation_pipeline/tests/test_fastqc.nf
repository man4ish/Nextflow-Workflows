nextflow.enable.dsl=2

process fastqc {
    input:
    path reads

    output:
    path "fastqc_report"

    container 'quay.io/biocontainers/fastqc:v0.11.9--0'

    script:
    """
    fastqc $reads --outdir=fastqc_report
    """
}

workflow {
    reads = file('test_data/sample.fastq')

    fastqc(reads)
}
