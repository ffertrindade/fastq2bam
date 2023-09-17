# Snakemake workflow: `fastq2bam`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)


A Snakemake workflow for `Processing Illumina raw short read data and mapping to a reference genome`.


It trims the reads using fastp [version used for test 0.23.4] (Chen, 2023) and checks the filtered reads with FastQC. Then, it independently maps the filtered reads (paired and colapsed) using BWA mem  [version used for test 0.7.12-r1039] (Li and Durbin, 2009). Next, it merges these BAM files (from the same sample and across multiple lanes), removes duplicates, and indexes them using Sambamba  [version used for test 0.6.6] (Tarasov et al., 2015). After that, it applies quality and region filtering using SAMtools [version used for test 0.1.20] (Danecek et al., 2021). It checks the mapping status at each step using SAMtools flagstat, deepTools  [version used for test 3.5.2] (Ramírez et al., 2016), and Qualimap  [version used for test v.2.2.2-dev] (Okonechnikov et al., 2016).



## Usage

1 - Modify the config/config.yaml file properlly with your files and parameters.

2 - Edit the "configfile:" line of the workflow/snakefile with respective config file name.

3 - The simplest way of running is just setting the number of threads <snakemake --cores 2>


If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository.
