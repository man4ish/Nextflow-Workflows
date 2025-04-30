nextflow.enable.dsl=2

process test_methylkit_analysis {
    input:
    path bam_files

    output:
    path "methylation_results.csv"

    script:
    """
    # Run the methylkit_analysis process
    nextflow run test_methylkit.nf --bam_files $bam_files
    """
}

workflow {
    // Test the process with mocked bam file
    bam_files = file('test_data/mock_bam.bam')

    test_methylkit_analysis(bam_files)
}
