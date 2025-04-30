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

process bismark_align {
    input:
    path fastq_files

    output:
    path "aligned_bam"

    container 'quay.io/biocontainers/bismark:v0.23.0--h7e3b6e0_0'

    script:
    """
    bismark --genome /path/to/genome --output_dir ./ $fastq_files
    """
}

process methylkit_analysis {
    input:
    path bam_files

    output:
    path "methylation_results.csv"

    container 'quay.io/biocontainers/r-methylkit:v1.9.0--r40hfc51f62_0'

    script:
    """
    Rscript methylkit_analysis.R $bam_files
    """
}

workflow {
    reads = file('test_data/sample.fastq')

    fastqc(reads)
    bismark_align(reads)
    methylkit_analysis(bam_files)
}
