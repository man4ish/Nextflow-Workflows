#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Test function for fastqc process
process test_fastqc {
    input:
    path reads

    output:
    path "fastqc_report"

    // Mocking fastqc process
    script:
    """
    echo 'Mocked fastqc process' > fastqc_report
    """
}

// Define the test input FASTQ file
workflow {
    reads = file('test_data/sample.fastq')

    // Run the mocked fastqc process (test function)
    test_fastqc(reads)
}
