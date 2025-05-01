// ========================
// main.nf - Whole Genome Sequencing Pipeline
// ========================

nextflow.enable.dsl=2

params.reads = "data/*_R{1,2}.fastq.gz"
params.outdir = "results"
params.reference = "genome.fa"
params.snpeff_db = "GRCh38.86"

workflow {
  Channel.fromFilePairs(params.reads, flat: true).set { read_pairs }

  trimmed = trimmomatic(read_pairs)
  qc = fastqc(trimmed)
  aligned = bwa_align(trimmed, params.reference)
  sorted_bam = sort_index(aligned)
  deduped_bam = mark_duplicates(sorted_bam)
  variants = call_variants(deduped_bam, params.reference)
  annotated = annotate_snps(variants)
}

process trimmomatic {
  tag "Trimming ${pair_id}"
  input:
    tuple val(pair_id), file(reads)
  output:
    tuple val(pair_id), file("*_paired.fq.gz")
  script:
    """
    trimmomatic PE \
      ${reads[0]} ${reads[1]} \
      ${pair_id}_1_paired.fq.gz ${pair_id}_1_unpaired.fq.gz \
      ${pair_id}_2_paired.fq.gz ${pair_id}_2_unpaired.fq.gz \
      ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50
    """
}

process fastqc {
  tag "FastQC ${pair_id}"
  input:
    tuple val(pair_id), file(reads)
  output:
    file("*fastqc.zip")
  script:
    """
    fastqc ${reads.join(' ')}
    """
}

process bwa_align {
  tag "BWA ${pair_id}"
  input:
    tuple val(pair_id), file(reads)
    path reference
  output:
    tuple val(pair_id), file("${pair_id}.bam")
  script:
    """
    bwa mem $reference ${reads.join(' ')} | samtools view -Sb - > ${pair_id}.bam
    """
}

process sort_index {
  tag "SortIndex ${pair_id}"
  input:
    tuple val(pair_id), file("*.bam")
  output:
    tuple val(pair_id), file("${pair_id}.sorted.bam")
  script:
    """
    samtools sort -o ${pair_id}.sorted.bam ${pair_id}.bam
    samtools index ${pair_id}.sorted.bam
    """
}

process mark_duplicates {
  tag "Dedup ${pair_id}"
  input:
    tuple val(pair_id), file("*.sorted.bam")
  output:
    tuple val(pair_id), file("${pair_id}.dedup.bam")
  script:
    """
    picard MarkDuplicates I=${pair_id}.sorted.bam O=${pair_id}.dedup.bam M=${pair_id}.metrics.txt REMOVE_DUPLICATES=true
    """
}

process call_variants {
  tag "Variants ${pair_id}"
  input:
    tuple val(pair_id), file("*.dedup.bam")
    path reference
  output:
    tuple val(pair_id), file("${pair_id}.vcf")
  script:
    """
    gatk HaplotypeCaller -R $reference -I ${pair_id}.dedup.bam -O ${pair_id}.vcf
    """
}

process annotate_snps {
  tag "SnpEff ${pair_id}"
  input:
    tuple val(pair_id), file("*.vcf")
  output:
    file("${pair_id}.annotated.vcf")
  script:
    """
    snpeff ${params.snpeff_db} ${pair_id}.vcf > ${pair_id}.annotated.vcf
    """
}
