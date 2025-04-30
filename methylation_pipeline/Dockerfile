# Start with a base image that includes Bioconductor (for R) and necessary dependencies
FROM rocker/r-ver:4.0.3

# Install system dependencies
RUN apt-get update && apt-get install -y \
    fastqc \
    bismark \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install R packages (MethylKit)
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('methylKit')"

# Set working directory
WORKDIR /data

# Set entrypoint
ENTRYPOINT ["bash"]
