#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define parameters with defaults
params.input           = "data"
params.output          = "results"
params.fastqc_threads  = 8
params.reference       = "/ref/genome.fa"
params.dbSNP           = false

// Reference file validation
if (!params.reference || !file(params.reference).exists()) {
    error "❌ Reference genome file '${params.reference}' does not exist or is not specified."
}

// Define channels
Channel
    .fromFilePairs("${params.input}/*_{R1,R2}.fastq.gz")
    .map { sample_id, reads -> tuple(sample_id, reads.toList()) }
    .set { read_pairs }

// FASTQC: Run quality check on paired-end reads
process FastQC {
    tag "$sample_id"
    publishDir "${params.output}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_fastqc.zip")

    script:
    """
    fastqc --nogroup -t ${params.fastqc_threads} ${reads[0]} ${reads[1]}
    """
}

// Validate reference file on host (local)
if (!params.reference || !file(params.reference).exists()) {
    error "Reference genome file '${params.reference}' does not exist or is not specified."
}

// ...

process bwa_align {
    tag "$sample_id"
    input:
        tuple val(sample_id), path(reads)
        path reference  // this is local path ref/genome.fa
    output:
        tuple val(sample_id), path("${params.output}/bam/${sample_id}/${sample_id}.bam")

    script:
    """
    ls -l /ref/genome.fa* || echo "Reference files not found in /ref"
    mkdir -p ${params.output}/bam/${sample_id}
    bwa mem -t ${params.fastqc_threads} /ref/genome.fa ${reads[0]} ${reads[1]} | \
    samtools view -bS - > ${params.output}/bam/${sample_id}/${sample_id}.bam
    """
}

process sort_bam {
    tag "sort_bam"
    input:
        tuple val(sample_id), path(bam)

    output:
        tuple val(sample_id), path("${sample_id}.sorted.bam")

    script:
    """
    samtools sort -@ 4 -o ${sample_id}.sorted.bam ${bam}
    """
}

// Process 3: MarkDuplicates
process mark_duplicates {
    tag 'mark_duplicates'
    input:
        tuple val(sample_id), path(bamfile)

    output:
        tuple val(sample_id), path("${params.output}/bam/${sample_id}/${sample_id}.dup.bam")

    script:
    """
    mkdir -p ${params.output}/bam/${sample_id}
    mkdir -p ${params.output}/dup_log
    java -jar /opt/picard.jar MarkDuplicates \
        INPUT=${bamfile} \
        OUTPUT=${params.output}/bam/${sample_id}/${sample_id}.dup.bam \
        METRICS_FILE=${params.output}/dup_log/${sample_id}_metrics.txt \
        VALIDATION_STRINGENCY=LENIENT \
        REMOVE_DUPLICATES=false \
        ASSUME_SORTED=false
    """
}

process base_recalibrator {
    tag "$sample_id"
    
    input:
        tuple val(sample_id), path(bamfile)

    output:
        tuple val(sample_id), path("${params.output}/bam/${sample_id}/${sample_id}_recal_data.csv")

    script:
    """
    gatk --java-options '-Xmx4g' BaseRecalibrator \
        -R ${params.reference} \
        -I ${bamfile} \
        --known-sites ${params.dbSNP} \
        -O ${params.output}/bam/${sample_id}/${sample_id}_recal_data.csv
    """
}

// Workflow
workflow {
    // Run FastQC
    fastqc_out = FastQC(read_pairs)

    // Run BWA
    bwa_out = bwa_align(read_pairs, file(params.reference))

    // Sort
    sorted_ch = sort_bam(bwa_out)

     // Mark Duplicates
    dup_ch = mark_duplicates(sorted_ch)


    // Base Recalibrator
    recal_ch = base_recalibrator(dup_ch)

    // Optional views for progress tracking
    fastqc_out.view { it -> "✔ FASTQC complete: ${it[0]}" }
    bwa_out.view    { it -> "✔ BWA complete: ${it[0]}" }
    dup_ch.view     { it -> "✔ Duplicates marked: ${it[0]}" }
    recal_ch.view   { it -> "✔ Base recalibration data generated: ${it[0]}" }
}
