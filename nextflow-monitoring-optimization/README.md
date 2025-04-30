## 🚀 Pipeline Monitoring & Optimization Guide

This section demonstrates how performance monitoring and basic optimizations are implemented in this Nextflow pipeline.

---

### 📊 Monitoring Commands

You can run the pipeline with built-in monitoring like this:

```bash
nextflow run main.nf \
  -with-report report.html \
  -with-trace trace.txt \
  -with-timeline timeline.html \
  -with-dag dag.png \
  -resume
```
These options enable:

- report.html: Summary of pipeline run

- trace.txt: Per-process performance data

- timeline.html: Visual timeline of job execution

- dag.png: Workflow diagram

-resume: Caching to skip already completed tasks

### Optimization Practices Used
Resource control per process
Each process defines cpus, memory, and time.

#### Caching
Using -resume reuses previous results for unchanged inputs.

#### Retries
Set maxRetries = 2 and errorStrategy = 'retry' for robustness.

#### Containerization
Each process runs inside an Ubuntu 20.04 container for consistency.

#### Parallelism
If input is split across multiple files, the pipeline can run in parallel using Channel.fromPath.

