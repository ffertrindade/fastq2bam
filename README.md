# Snakemake workflow: `fastq2bam`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.1.0-brightgreen.svg)](https://snakemake.github.io)


A Snakemake workflow for `Processing Illumina raw short read data and mapping to a reference genome`.


It trims the reads using fastp v0.23.2 (Chen, 2023) and checks the filtered reads with FastQC. Then, it independently maps the filtered reads (paired and colapsed) using BWA mem v0.7.17 (Li and Durbin, 2009). Next, it merges these BAM files (from the same sample and across multiple lanes), removes duplicates, and indexes them using Sambamba v0.6.6 (Tarasov et al., 2015). After that, it applies quality and region filtering using SAMtools v1.5 (Danecek et al., 2021). It checks the mapping status at each step using SAMtools flagstat, deepTools v2.5.7 (Ramírez et al., 2016), and Qualimap v2.2.2a (Okonechnikov et al., 2016).



## Usage

1 - Modify the [config.yaml](https://github.com/ffertrindade/fastq2bam/blob/main/config/config.yaml) file properlly with your files and parameters.

1a - You can use the script [creatingYamlRawFastq.py](https://github.com/ffertrindade/fastq2bam/blob/main/scripts/creatingYamlRawFastq.py) to assist you in creating the config.yaml for several input samples.

2 - The simplest way of running is `snakemake --snakefile workflow/Snakefile --use-conda --cores 8`.

2a - For testing, you can use the example files by running `snakemake -d example/ --snakefile workflow/Snakefile --use-conda --cores 8`.


If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository.
