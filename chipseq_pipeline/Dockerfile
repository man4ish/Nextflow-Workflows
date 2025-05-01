# Use an official bioinformatics base image as a starting point
FROM ubuntu:20.04

# Set environment variables to avoid interactive prompts during installation
ENV DEBIAN_FRONTEND=noninteractive

# Update package lists and install required dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    curl \
    git \
    python3 \
    python3-pip \
    samtools \
    bowtie2 \
    fastqc \
    bedtools \
    macs2 \
    default-jdk \
    cutadapt \
    python3-setuptools \
    python3-dev \
    unzip \
    bwa \
    sudo \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    && apt-get clean

# Install HOMER (Motif Analysis)
RUN curl -L https://github.com/Defaulty/homer/archive/refs/tags/v4.11.1.tar.gz -o homer.tar.gz \
    && tar -xzvf homer.tar.gz \
    && cd homer-4.11.1 \
    && make

# Install Nextflow (for running the pipeline)
RUN curl -s https://get.nextflow.io | bash \
    && mv nextflow /usr/local/bin/nextflow

# Install Conda (optional for managing Python dependencies)
RUN curl -O https://repo.anaconda.com/archive/Anaconda3-2023.03-1-Linux-x86_64.sh \
    && bash Anaconda3-2023.03-1-Linux-x86_64.sh -b -p /opt/conda \
    && rm Anaconda3-2023.03-1-Linux-x86_64.sh \
    && /opt/conda/bin/conda init bash \
    && echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc

# Install Python dependencies
RUN pip3 install numpy pandas matplotlib seaborn

# Set working directory
WORKDIR /workspace

# Set the entrypoint to run Nextflow
ENTRYPOINT ["nextflow"]

# Expose ports for Nextflow if running locally with web interface
EXPOSE 8080

# Optional: Configure default CMD if no command is passed
CMD ["-v"]
