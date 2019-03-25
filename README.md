# RNA-SEQ analysis

# Requeriments

- MacOS or Linux
- Docker (https://www.docker.com/)

# Build DockerFile
Build DockerFile.
```
cd docker
docker build . -t carlgira/rna-analysis:latest
```

## Tools
- fastqc 0.11.8
- hisat2 2.1.0
- samtools 1.9
- subread 1.6.4

# Run pipeline
- Download samples (TODO)
- Data preparation
- Quality Control
- Alignment
- Counting and Differential Expression
- Differential Expression Stats

## Download samples (TODO)
Create a folder called "reads" and download all the fastq files from FTP. (not public url yet)

## Data preparation
- Download reference genome from aedes. (5.0v)
- Download GTF file from aedes. (5.1v)
- Build Index for reference genome
- Prepare samples for procesing
```
sh scripts/prepare-data.sh
```

## Quality Control
- Run *fastqc* to get quality metrics for all the fastq samples.
```
sh scripts/quality.sh
```

## Alignment
- Align each sample to the reference genome using *hisat2*
```
sh scripts/align.sh
```

## Counting and Differential Expression
- Counting using *featureCounts* and differential expression of samples using *edgeR*
```
sh scripts/count.sh
```

## Differential Expression Stats
- Genereate stats
```
sh scripts/stats.sh
```
