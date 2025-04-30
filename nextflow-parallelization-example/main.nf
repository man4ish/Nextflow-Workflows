#!/usr/bin/env nextflow

params.input = 'input/*.txt'
params.outdir = 'results'

Channel.fromPath(params.input).set { input_files }

process PROCESS_FILE {
    cpus 1
    memory '500 MB'
    time '30m'

    input:
    path file

    output:
    path "${file.simpleName}.processed.txt"

    script:
    """
    echo "Processed file: $file" > ${file.simpleName}.processed.txt
    """
}

workflow {
    PROCESS_FILE(input_files)
}
