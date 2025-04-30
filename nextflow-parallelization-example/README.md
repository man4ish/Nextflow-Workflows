## ⚡ Parallelization Demo with Nextflow

This example demonstrates how to process multiple input files **in parallel** using Nextflow.

---

### 📂 Input Files

All `.txt` files inside the `input/` folder are picked up and passed to a process:

```bash
input/
├── sample1.txt
├── sample2.txt
└── sample3.txt
```

## Run the Pipeline
```
nextflow run main.nf -with-trace -resume
```

### How Parallelization Works
Each file is passed individually into the PROCESS_FILE process, which runs independently for each input. Nextflow handles:

### Task spawning

### Job queuing

Automatic parallel execution

You can adjust parallelism with:

```
process {
  maxForks = 4
  queueSize = 100
}
```

### Output
Processed files will be written to results/:

```

results/
├── sample1.processed.txt
├── sample2.processed.txt
└── sample3.processed.txt

```
