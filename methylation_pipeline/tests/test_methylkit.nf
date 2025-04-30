nextflow.enable.dsl=2

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
    bam_files = file('test_data/aligned_bam.bam')

    methylkit_analysis(bam_files)
}
