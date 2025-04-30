nextflow.enable.dsl=2

workflow {

    params.prefix = params.prefix ?: "sample"

    trimmed = trim(
        fastq1: file(params.fastq1),
        fastq2: file(params.fastq2),
        adapters: file(params.adapters),
        skewer_threads: params.skewer_threads,
        minimum_read_length: params.minimum_read_length,
        prefix: params.prefix
    )

    kallisto = quantification(
        fastq1: trimmed.out_pair1,
        fastq2: trimmed.out_pair2,
        idx: file(params.idx),
        gtf: file(params.gtf),
        kallisto_threads: params.kallisto_threads,
        bootstrap_samples: params.bootstrap_samples
    )

    aligned = align(
        ref_tar: file(params.ref_tar),
        fastq1: trimmed.out_pair1,
        fastq2: trimmed.out_pair2,
        STAR_threads: params.STAR_threads
    )

    sorted = sort_index(
        bam_file: aligned.out_bam
    )

    gen_summary(
        bam_file: sorted.sorted_bam,
        ref_flat: file(params.ref_flat),
        ribosomal_interval: file(params.ribosomal_interval),
        ref_seq: file(params.ref_seq)
    )
}

process trim {
    input:
    file fastq1
    file fastq2
    file adapters
    val skewer_threads
    val minimum_read_length
    val prefix

    output:
    file "${prefix}-trimmed-pair1.fastq" into trimmed_out1
    file "${prefix}-trimmed-pair2.fastq" into trimmed_out2

    container 'docker.io/man4ish/rnaseq:latest'

    script:
    """
    skewer -t ${skewer_threads} -y ${adapters} -l ${minimum_read_length} ${fastq1} ${fastq2} -o ${prefix}
    """
}

process quantification {
    input:
    file fastq1
    file fastq2
    file idx
    file gtf
    val kallisto_threads
    val bootstrap_samples

    output:
    file "*.tar.gz"

    container 'docker.io/man4ish/rnaseq:latest'

    script:
    def base = fastq1.getSimpleName().tokenize('-')[0]
    """
    mkdir ${base}_output
    /software/utils/kallisto quant -t ${kallisto_threads} -b ${bootstrap_samples} \
        --rf-stranded --genomebam -i ${idx} -g ${gtf} ${fastq1} ${fastq2} -o ${base}_output
    tar -czf ${base}_kallisto_output.tar.gz ${base}_output
    """
}

process align {
    input:
    file ref_tar
    file fastq1
    file fastq2
    val STAR_threads

    output:
    file "*sample.bam" into aligned_bam

    container 'docker.io/man4ish/rnaseq:latest'

    script:
    def base = ref_tar.getSimpleName().replace(".tar.gz", "")
    """
    mkdir ref_bundle
    tar -xzf ${ref_tar} -C ref_bundle --no-same-owner
    mv ref_bundle/*/* ref_bundle
    /software/utils/STAR --genomeDir ref_bundle --runThreadN ${STAR_threads} \
        --outSAMtype BAM Unsorted --readFilesIn ${fastq1} ${fastq2} --outFileNamePrefix ${base}_sample
    cp ${base}_sampleAligned.out.bam ${base}_sample.bam
    """
}

process sort_index {
    input:
    file bam_file

    output:
    file "sorted_${bam_file.simpleName}" into sorted_bam
    file "sorted_${bam_file.simpleName}.bai"

    container 'docker.io/man4ish/rnaseq:latest'

    script:
    """
    samtools sort ${bam_file} > sorted_${bam_file.simpleName}
    samtools index sorted_${bam_file.simpleName}
    """
}

process gen_summary {
    input:
    file bam_file
    file ref_flat
    file ribosomal_interval
    file ref_seq

    output:
    file "*.summary"
    file "*.plot.pdf"

    container 'docker.io/man4ish/rnaseq:latest'

    script:
    def base = bam_file.getSimpleName()
    """
    java -jar /software/utils/picard.jar CollectRnaSeqMetrics METRIC_ACCUMULATION_LEVEL=ALL_READS \
        REF_FLAT=${ref_flat} RIBOSOMAL_INTERVALS=${ribosomal_interval} \
        STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
        CHART_OUTPUT=${base}_position.vs.coverage.plot.pdf \
        INPUT=${bam_file} OUTPUT=${base}.rna.summary \
        REFERENCE_SEQUENCE=${ref_seq} VALIDATION_STRINGENCY=LENIENT
    """
}