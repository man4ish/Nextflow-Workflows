# 🧪 Nextflow Pipeline Testing Demo

This repository showcases different types of testing strategies for Nextflow pipelines. It includes examples of:

- ✅ Unit Testing
- 🔁 Integration Testing
- 🔍 Functional Testing
- 💥 Edge Case & Failure Testing
- nf test
- Regression Testing
- ⏱️ Performance Testing (conceptual setup)
- 🚀 CI Testing using GitHub Actions

---

## 📁 Project Structure

| Folder              | Purpose                                  |
|---------------------|------------------------------------------|
| `main.nf`           | Core pipeline with minimal logic         |
| `unit_tests/`       | Test individual processes                |
| `integration_tests/`| Full pipeline with test data             |
| `nf-test/`          | [nf-test](https://github.com/seqeralabs/nf-test) testing framework |
| `test_data/`        | Tiny mock input datasets                 |
| `.github/workflows/`| GitHub CI to run tests automatically     |

---

## ▶️ How to Run

### 🔹 Run Unit Test
```bash
nextflow run unit_tests/process_test.nf
```
### Run Integration Test
```bash
nextflow run integration_tests/test_pipeline.nf -profile test
```
###  Run Full Pipeline
``` bash
nextflow run main.nf -profile test --input test_data/input.txt
```
## nf-test
``` bash
nf-test run nf-test/basic.nftest
```

## Regression Testing
### 1. Establish a Baseline
Before making changes to your pipeline, you need to have a set of known-good outputs from your pipeline. These are the results from a stable version of the pipeline that you can compare against after new updates.

### Example:
Run your pipeline with a sample input (e.g., input.txt) and save the output (output.txt).

Store the output in a separate directory or a versioned artifact repository (e.g., Git, cloud storage, or even locally).

### 2. Make Changes to the Pipeline
Modify your pipeline (e.g., adding new features, fixing bugs, upgrading dependencies, or changing logic).

### 3. Run Regression Tests
Run the pipeline with the same input(s) and compare the new outputs with the baseline outputs.

Example Command:
```bash
nextflow run main.nf --input tests/data/input.txt
```
### 4. Compare Results
Compare the new output against the old output from the baseline:

Manual Comparison:
Compare file contents: You can compare the new and old files using diff or cmp.

```bash
diff old_output.txt new_output.txt
```
### 5. Automated Comparison:
Use a test script to automate this comparison:
```
bash
#!/bin/bash
# Run the pipeline with new changes
nextflow run main.nf --input tests/data/input.txt

# Compare new output with baseline output
if diff -q old_output.txt new_output.txt; then
  echo "Regression test passed: Outputs are the same."
else
  echo "Regression test failed: Outputs differ."
  exit 1
fi
```
### 6. Automate with CI/CD
You can automate regression testing by integrating it into a CI/CD pipeline (using tools like GitHub Actions, GitLab CI, or Jenkins). This way, every time you push changes to your pipeline, the regression tests run automatically.

Example GitHub Action for Regression Testing:

```
name: Regression Testing

on: [push, pull_request]

jobs:
  regression_test:
    runs-on: ubuntu-latest
    steps:
      - name: Checkout code
        uses: actions/checkout@v3
      - name: Install Nextflow
        run: curl -s https://get.nextflow.io | bash
      - name: Run regression tests
        run: |
          nextflow run main.nf --input tests/data/input.txt
          diff -q old_output.txt new_output.txt || exit 1
```

## CI with GitHub Actions
```bash
nextflow run main.nf -profile test --input test_data/input.txt
```
