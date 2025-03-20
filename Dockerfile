FROM condaforge/mambaforge:latest

RUN apt-get update && apt-get install -y \
    procps \
    wget \
    curl \
    grep \
    sed \
    unzip \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

COPY env.yml /

RUN mamba env create -f /env.yml && \
    mamba clean -a -y && \
    echo "Installed packages in environment:" && \
    conda list -n intrahost_analysis

ENV PATH /opt/conda/envs/intrahost_analysis/bin:$PATH

RUN echo "source activate intrahost_analysis" >> ~/.bashrc

SHELL ["conda", "run", "--no-capture-output", "-n", "intrahost_analysis", "/bin/bash", "-c"]