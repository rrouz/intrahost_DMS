FROM golang:1.23-alpine AS builder

# Install git and build dependencies
RUN apk add --no-cache git

RUN git clone https://github.com/rrouz/gofasta.git /go/src/github.com/virus-evolution/gofasta && \
    cd /go/src/github.com/virus-evolution/gofasta && \
    git checkout 4dd35a5 && \
    git pull origin master

# Build gofasta
WORKDIR /go/src/github.com/virus-evolution/gofasta
RUN go build -o /go/bin/gofasta .

# Main image
FROM condaforge/mambaforge:latest

# Copy the compiled gofasta binary from the builder stage
COPY --from=builder /go/bin/gofasta /usr/local/bin/gofasta

# Combine apt-get operations into a single layer
RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    grep \
    sed \
    minimap2 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

COPY env.yml /

RUN mamba env create -f /env.yml && \
    mamba clean -a -y

ENV PATH /opt/conda/envs/intrahost_analysis/bin:$PATH

RUN echo "source activate intrahost_analysis" >> ~/.bashrc

SHELL ["conda", "run", "--no-capture-output", "-n", "intrahost_analysis", "/bin/bash", "-c"]