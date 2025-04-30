#!/usr/bin/env nextflow

params.input = "input/sample.txt"
params.outdir = "results"

process SAY_HELLO {
    cpus 2
    memory '2 GB'
    time '1h'

    input:
    path sample_file

    output:
    path "hello.txt"

    script:
    """
    echo "Hello, this is \$(hostname) processing \${sample_file}" > hello.txt
    """
}

workflow {
    Channel.fromPath(params.input).set { input_files }
    SAY_HELLO(input_files)
}
