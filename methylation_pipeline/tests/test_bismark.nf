nextflow.enable.dsl=2

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

workflow {
    fastq_files = file('test_data/sample.fastq')

    bismark_align(fastq_files)
}
