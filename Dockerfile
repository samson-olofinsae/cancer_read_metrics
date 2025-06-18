# Dockerfile for Cancer Read Metrics Project

FROM continuumio/miniconda3

# Create and activate environment
RUN conda create -n cancer-metrics python=3.10 -y && \
    echo "conda activate cancer-metrics" >> ~/.bashrc

# Activate env, install samtools and other tools
RUN /bin/bash -c "source ~/.bashrc && \
    conda activate cancer-metrics && \
    conda install -c bioconda -c conda-forge samtools=1.18 pandas matplotlib seaborn -y"

# Set working directory inside container
WORKDIR /workspace

# Default command when container starts
CMD ["/bin/bash"]

