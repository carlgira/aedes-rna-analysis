# RNA-SEQ analysis

# Requeriments

- Docker (https://www.docker.com/)

# Build DockerFile
Build DockerFile.
```
cd docker
docker build . -t carlgira/rna-analysis:latest
alias drun="docker run --memory="28000m" -it --rm -v $(pwd):/work -w /work carlgira/rna-analysis:latest /bin/bash"
alias drun="docker run -it --rm -v $(pwd):/work -w /work carlgira/rna-analysis:latest /bin/bash"
```

## Tools
- fastqc 0.11.8
- hisat2 2.1.0
- samtools 1.9
- subread 1.6.4
- multiqc 1.7
- star 2.7.0e
- edger 3.24.1
- biomart 2.38.0
- noiseq 2.26.0

# Run pipeline
- Download samples (TODO)
- Data preparation
- Quality Control
- Alignment
- Counting and Differential Expression (R notebooks)

For re-running from the beginning
```
drun scripts/run-all.sh
```

For re-running only the Differential Expression analysis
```
drun scripts/run-all.sh DE
```

# References

- A comprehensive evaluation of normalization methods for Illumina high-throughput RNA sequencing data analysis https://academic.oup.com/bib/article/14/6/671/189645
- A survey of best practices for RNA-seq data analysis https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4728800/
- Noiseq https://bioconductor.org/packages/release/bioc/vignettes/NOISeq/inst/doc/NOISeq.pdf
- RNAseq analysis in R https://bioinformatics-core-shared-training.github.io/RNAseq-R/
- START Aligner https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html
