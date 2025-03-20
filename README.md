# H5N1 DMS Analysis Pipeline
A Nextflow pipeline for analyzing deep mutational scanning (DMS) data from H5N1 influenza sequences.

## Overview
This pipeline processes influenza sequence data through multiple steps:
1. Reference sequence extraction
2. Sequence alignment
3. Protein translation
4. Mutation identification
5. DMS data integration

## Requirements
- Nextflow (>=23.04.0)
- Java 11 or later
- Conda

## Data Requirements
1. Obtain sequences and variants:
```bash
wget https://github.com/andersen-lab/avian-influenza/raw/master/fasta/
wget https://github.com/andersen-lab/avian-influenza/raw/master/variants/
```

# Usage
```bash
nextflow run main.nf
```

## Docker
```bash
nextflow run main.nf -profile docker
``` 
