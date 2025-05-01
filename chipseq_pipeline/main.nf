#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define parameters for the pipeline
params.reads = "./data/*_R1.fastq"
params.reference = "hg38" // Reference genome (can be customized)
params.control = "./data/control.bam" // Control sample (input)
params.outdir = "./results" // Output directory

// Define the processes

process fastqc {
    input:
    path fastq_files

    output:
    path "${task.workDir}/fastqc_reports/*"

    script:
    """
    mkdir -p ${task.workDir}/fastqc_reports
    fastqc -o ${task.workDir}/fastqc_reports $fastq_files
    """
}

process trim_reads {
    input:
    path fastq_files

    output:
    path "trimmed/*"

    script:
    """
    mkdir -p trimmed
    cutadapt -a AGATCGGAAGAGC -o trimmed/$(basename ${fastq_files}) ${fastq_files}
    """
}

process align_bowtie2 {
    input:
    path fastq_files
    val reference_genome

    output:
    path "aligned/*"

    script:
    """
    mkdir -p aligned
    bowtie2 -x $reference_genome -U $fastq_files -S aligned/$(basename ${fastq_files}.sam)
    """
}

process peak_calling_macs {
    input:
    path aligned_sam_files
    val control_bam

    output:
    path "peaks/*"

    script:
    """
    mkdir -p peaks
    macs2 callpeak -t $aligned_sam_files -c $control_bam -f SAM -g hs -n $(basename ${aligned_sam_files}.bam) --outdir peaks
    """
}

process motif_analysis {
    input:
    path peak_files

    output:
    path "motif_analysis/*"

    script:
    """
    mkdir -p motif_analysis
    findMotifsGenome.pl $peak_files hg38 motif_analysis/
    """
}

// Define the workflow
workflow {
    // Input files and parameters
    fastq_files = file(params.reads)
    reference_genome = params.reference
    control_bam = file(params.control)

    // Execution flow
    fastqc(fastq_files)
    trimmed_fastq = trim_reads(fastq_files)
    aligned_sam = align_bowtie2(trimmed_fastq, reference_genome)
    peak_files = peak_calling_macs(aligned_sam, control_bam)
    motif_analysis(peak_files)
}
