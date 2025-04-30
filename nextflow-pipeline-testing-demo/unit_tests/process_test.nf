process process_greet {
    input:
        path x
    output:
        path "hello.txt"
    script:
        """
        echo "Hello from \$(cat $x)" > hello.txt
        """
}

workflow {
    Channel.fromPath('test_data/input.txt') | process_greet
}
