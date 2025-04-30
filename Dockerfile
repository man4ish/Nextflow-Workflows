# Use a base image with Java and other dependencies
FROM openjdk:11-jre-slim

# Install necessary tools: BWA, Samtools, Picard, GATK, FastQC, etc.
RUN apt-get update && apt-get install -y \
    bwa \
    samtools \
    fastqc \
    curl \
    unzip \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install GATK
RUN curl -L https://github.com/broadinstitute/gatk/releases/download/4.2.0.0/gatk-4.2.0.0.zip -o gatk.zip \
    && unzip gatk.zip -d /opt/ \
    && rm gatk.zip

# Install Picard
RUN curl -L https://github.com/broadinstitute/picard/releases/download/2.27.1/picard.jar -o /opt/picard.jar

# Set GATK and Picard as environment variables
ENV PATH="/opt/gatk-4.2.0.0:${PATH}"
ENV PICARD="/opt/picard.jar"

# Default command to run the Nextflow pipeline
CMD ["nextflow"]
