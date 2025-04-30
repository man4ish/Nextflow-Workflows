
---

### 🔹 4. **Sample Nextflow Files**

#### `main.nf`
```nextflow
workflow {
    Channel.fromPath(params.input)
        | process_greet
}

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
