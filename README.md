# RNA-SEQ analysis

# Requeriments

- MacOS or Linux
- Docker (https://www.docker.com/)

# Build DockerFile
Build DockerFile.
```
cd docker
docker build . -t carlgira/rna-analysis:latest
alias drun="docker run -it --rm -v $(pwd):/work -w /work carlgira/rna-analysis:latest /bin/bash"
alias drun="docker run -it --rm -v $(pwd):/work -w /work carlgira/rna-analysis-u:latest /bin/bash"

alias drun="docker run --memory="28000m" -it --rm -v $(pwd):/work -w /work carlgira/rna-analysis-u:latest /bin/bash"
```

## Tools
- fastqc 0.11.8
- hisat2 2.1.0
- samtools 1.9
- subread 1.6.4
- multiqc
-

# Run pipeline
- Download samples (TODO)
- Data preparation (1 hr)
- Quality Control (2 min)
- Alignment (10 min for each sample, 2.5 hr)
- Counting and Differential Expression (5 min)
- Differential Expression Stats (1 min)

## Download samples (TODO)
Create a folder called "reads" and download all the fastq files from FTP. (not public url yet)

## Data preparation
- Download reference genome from aedes (5.0v)
- Download GTF file from aedes (5.1v)
- Build Index for reference genome
- Prepare samples for procesing
```
drun scripts/prepare-data.sh
```

## Quality Control
- Run **fastqc** to get quality metrics for all the fastq samples.
```
drun scripts/quality.sh
```
Results in *quality* folder.

## Alignment and Counting
- Align each sample to the reference genome using **hisat2** and count with **featureCounts**
```
drun scripts/align-and-count.sh c3a_vs_control setup/analysis_samples_vitro_exp1.tsv
drun scripts/align-and-count.sh c5a_vs_control setup/analysis_samples_vitro_exp2.tsv

drun scripts/align-and-count.sh normal_vs_sugarfed setup/analysis_samples_vivo_exp1.tsv
drun scripts/align-and-count.sh inactivated_vs_sugarfed setup/analysis_samples_vivo_exp2.tsv
drun scripts/align-and-count.sh inactivated_vs_normal setup/analysis_samples_vivo_exp3.tsv
```
Results in *bam* folder.

## Counting and Differential Expression
- Counting using **featureCounts** and differential expression of samples using **edgeR**
```
drun scripts/count.sh
```
Results in *count* folder.

## Differential Expression Stats
- Genereate stats
```
drun scripts/stats.sh
```
Results in *stats* folder.

# References

- A comprehensive evaluation of normalization methods for Illumina high-throughput RNA sequencing data analysis https://academic.oup.com/bib/article/14/6/671/189645
- A survey of best practices for RNA-seq data analysis https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4728800/
- Noiseq https://bioconductor.org/packages/release/bioc/vignettes/NOISeq/inst/doc/NOISeq.pdf
- RNAseq analysis in R https://bioinformatics-core-shared-training.github.io/RNAseq-R/
- START Aligner https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html


- Aedes Genome https://www.ncbi.nlm.nih.gov/assembly/GCF_002204515.2

picard CollectAlignmentSummaryMetrics -Xmx2G R=refs/Aedes-aegypti-LVP_AGWG_CHROMOSOMES_AaegL5.fa I=bam/Normal1_mappable.bam O=output.txt

picard CollectGcBiasMetrics -Xmx2G I=bam/Normal1_mappable.bam O=gc_bias_metrics.txt CHART=gc_bias_metrics.pdf S=summary_metrics.txt R=refs/Aedes-aegypti-LVP_AGWG_CHROMOSOMES_AaegL5.fa
