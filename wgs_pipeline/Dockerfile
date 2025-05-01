FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8

# Basic system dependencies
RUN apt-get update && apt-get install -y \
    wget curl unzip git build-essential \
    openjdk-11-jdk python3 python3-pip zlib1g-dev \
    samtools libbz2-dev liblzma-dev locales \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Set JAVA_HOME
ENV JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64
ENV PATH=$JAVA_HOME/bin:$PATH

# -------------------------
# Install FastQC
# -------------------------
RUN wget -q https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip && \
    chmod +x FastQC/fastqc && \
    mv FastQC /opt/fastqc && \
    ln -s /opt/fastqc/fastqc /usr/local/bin/fastqc

# -------------------------
# Install Trimmomatic
# -------------------------
RUN wget -q https://github.com/usadellab/Trimmomatic/releases/download/v0.39/Trimmomatic-0.39.zip && \
    unzip Trimmomatic-0.39.zip && \
    mv Trimmomatic-0.39 /opt/trimmomatic

ENV TRIMMOMATIC_JAR=/opt/trimmomatic/trimmomatic-0.39.jar

# -------------------------
# Install BWA
# -------------------------
RUN git clone https://github.com/lh3/bwa.git /opt/bwa && \
    cd /opt/bwa && make && \
    ln -s /opt/bwa/bwa /usr/local/bin/bwa

# -------------------------
# Install Picard
# -------------------------
RUN wget -q https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar -O /opt/picard.jar

# -------------------------
# Install GATK
# -------------------------
RUN wget -q https://github.com/broadinstitute/gatk/releases/download/4.4.0.0/gatk-4.4.0.0.zip && \
    unzip gatk-4.4.0.0.zip && \
    mv gatk-4.4.0.0 /opt/gatk && \
    ln -s /opt/gatk/gatk /usr/local/bin/gatk

# -------------------------
# Install SnpEff
# -------------------------
RUN wget -q https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip && \
    unzip snpEff_latest_core.zip && \
    mv snpEff /opt/snpeff && \
    ln -s /opt/snpeff/snpEff.jar /usr/local/bin/snpEff.jar

# Add wrappers for Java tools
RUN echo '#!/bin/bash\njava -jar /opt/picard.jar "$@"' > /usr/local/bin/picard && chmod +x /usr/local/bin/picard
RUN echo '#!/bin/bash\njava -Xmx4g -jar /usr/local/bin/snpEff.jar "$@"' > /usr/local/bin/snpeff && chmod +x /usr/local/bin/snpeff
RUN echo '#!/bin/bash\njava -Xmx4g -jar $TRIMMOMATIC_JAR "$@"' > /usr/local/bin/trimmomatic && chmod +x /usr/local/bin/trimmomatic

WORKDIR /data
ENTRYPOINT ["/bin/bash"]
