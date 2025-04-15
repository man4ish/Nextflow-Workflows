#!/usr/bin/env nextflow

// Define paths for input and output
params.input = './in/non_redundant_fastq/'
params.output = './out'
params.reference = './in/reference/assembly38.fasta'
params.dbSNP = './in/reference/Homo_sapiens_assembly38.dbsnp138.vcf'

// Process 1: FastQC
process fastqc {
    tag 'fastqc'
    input:
    path read1, read2

    output:
    path 'out/fastqc/*'

    script:
    """
    fastqc -t 8 \$read1 \$read2 --outdir ./out/fastqc
    """
}

// Process 2: BWA alignment
process bwa_align {
    tag 'bwa'
    input:
    path read1, read2
    val name

    output:
    path "out/aligned_bam/${name}.bam"

    script:
    """
    bwa mem -t 8 $params.reference \$read1 \$read2 | samtools view -bS | samtools sort -o out/aligned_bam/${name}.bam
    """
}

// Process 3: MarkDuplicates
process mark_duplicates {
    tag 'mark_duplicates'
    input:
    path bamfile
    val name

    output:
    path "out/bamfiles/renamed_bam/${name}.dup.bam"

    script:
    """
    java -jar picard.jar MarkDuplicates \
        INPUT=\$bamfile \
        OUTPUT=out/bamfiles/renamed_bam/${name}.dup.bam \
        METRICS_FILE=out/bamfiles/dup_log/${name}_metrics.txt \
        VALIDATION_STRINGENCY=LENIENT
    """
}

// Process 4: BaseRecalibrator
process base_recalibrator {
    tag 'base_recalibrator'
    input:
    path bamfile
    val name

    output:
    path "out/bamfiles/renamed_bam/${name}_recal_data.csv"

    script:
    """
    gatk --java-options '-Xmx4g' BaseRecalibrator \
        -R $params.reference \
        -I \$bamfile \
        --known-sites $params.dbSNP \
        -O out/bamfiles/renamed_bam/${name}_recal_data.csv
    """
}

// Process 5: ApplyBQSR
process apply_bqsr {
    tag 'apply_bqsr'
    input:
    path bamfile, recal_data
    val name

    output:
    path "out/bamfiles/renamed_bam/${name}.recal.bam"

    script:
    """
    gatk --java-options '-Xmx4g' ApplyBQSR \
        -R $params.reference \
        -I \$bamfile \
        --bqsr-recal-file \$recal_data \
        -O out/bamfiles/renamed_bam/${name}.recal.bam
    """
}

// Process 6: HaplotypeCaller
process haplotype_caller {
    tag 'haplotype_caller'
    input:
    path bamfile
    val name

    output:
    path "out/bamfiles/renamed_bam/${name}.g.vcf.gz"

    script:
    """
    gatk HaplotypeCaller \
        -R $params.reference \
        -I \$bamfile \
        -O out/bamfiles/renamed_bam/${name}.g.vcf.gz \
        -ERC GVCF
    """
}

// Main workflow
workflow {
    reads = Channel.fromPath("${params.input}/*.fastq.gz")
        .splitCsv(sep: '_')

    reads
        .set { fastqc_input }

    fastqc(fastqc_input)
        .set { fastqc_results }

    fastqc_results
        .set { bwa_input }

    bwa_align(bwa_input)
        .set { aligned_bam }

    aligned_bam
        .set { mark_duplicates_input }

    mark_duplicates(mark_duplicates_input)
        .set { dedup_bam }

    dedup_bam
        .set { recalibration_input }

    base_recalibrator(recalibration_input)
        .set { recal_data }

    recal_data
        .set { apply_bqsr_input }

    apply_bqsr(apply_bqsr_input)
        .set { recal_bam }

    recal_bam
        .set { haplotype_input }

    haplotype_caller(haplotype_input)
}
