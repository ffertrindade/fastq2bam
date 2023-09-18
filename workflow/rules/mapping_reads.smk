rule bwa_map_paired:
	message: """##### Running bwa mem for mapping paired reads of {wildcards.sample}.{wildcards.lane}... #####"""
	input:
		{config["reference"]},
		"../results/trimmed_reads/{sample}.{lane}.trimmed_R1.fastq.gz",
		"../results/trimmed_reads/{sample}.{lane}.trimmed_R2.fastq.gz"
	output:
		"../results/mapped_reads/{sample}.{lane}.paired.bam"
	log:
		"../results/logs/mapping/{sample}.{lane}.bwa_map_paired.log"
	threads:
		config["threads"]
	params:
		rg=r"@RG\tID:{sample}\tSM:{sample}",
		prefix="../results/mapped_reads/{sample}.{lane}.paired.bam"
	shell:
		"(bwa mem -R '{params.rg}' -t {threads} {input} | "
		"samtools view -Shb -F 4 - | "
		"samtools sort -@ {threads} -f - {params.prefix}) 2> {log}"

rule bwa_map_colap:
	message: """##### Running bwa mem for mapping collapsed reads of {wildcards.sample}.{wildcards.lane}... #####"""
	input:
		{config["reference"]},
		"../results/trimmed_reads/{sample}.{lane}.trimmed_colap.fastq.gz"
	output:
		"../results/mapped_reads/{sample}.{lane}.colap.bam"
	log:
		"../results/logs/mapping/{sample}.{lane}.bwa_map_colap.log"
	threads:
		config["threads"]
	params:
		rg=r"@RG\tID:{sample}\tSM:{sample}",
		prefix="../results/mapped_reads/{sample}.{lane}.colap.bam"
	shell:
		"(bwa mem -R '{params.rg}' -t {threads} {input} | "
		"samtools view -Shb -F 4 - | "
		"samtools sort -@ {threads} -f - {params.prefix}) 2> {log}"

rule merge_raw_bam:
	message: """##### Merging collapsed and paired bams for {wildcards.sample}... #####"""
	input:
		pe=lambda wildcards: expand("../results/mapped_reads/{{sample}}.{lane}.paired.bam", lane=config["samples"][wildcards.sample]),
		se=lambda wildcards: expand("../results/mapped_reads/{{sample}}.{lane}.colap.bam", lane=config["samples"][wildcards.sample])
	output:
		temp("../results/mapped_reads/{sample}.sorted.merged.bam")
	log:
		"../results/logs/mapping/{sample}.merge_raw_bam.log"
	threads:
		config["threads"]
	shell:
		"sambamba merge -p -t {threads} {output} {input.pe} {input.se} 2> {log}"
