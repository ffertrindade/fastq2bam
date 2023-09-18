rule fastp:
	message: """##### Running fastp for {wildcards.sample}.{wildcards.lane}... #####"""
	input:
		pe1=lambda wildcards: config["samples"][wildcards.sample][wildcards.lane]["pe1"],
		pe2=lambda wildcards: config["samples"][wildcards.sample][wildcards.lane]["pe2"]
	output:
		pe1="../results/trimmed_reads/{sample}.{lane}.trimmed_R1.fastq.gz",
		pe2="../results/trimmed_reads/{sample}.{lane}.trimmed_R2.fastq.gz",
		colap="../results/trimmed_reads/{sample}.{lane}.trimmed_colap.fastq.gz"
	log:
		"../results/logs/trimming_reads/{sample}.{lane}.fastp.log"
	threads:
		config["threads"]
	params:
		config["fastpPar"]
	shell:
		"fastp --in1 {input.pe1} --in2 {input.pe2} "
		"{params} --thread {threads} "
		"--out1 {output.pe1} --out2 {output.pe2} --merged_out {output.colap} 2> {log}"

rule fastqc:
	message: """##### Running fastQC for trimmed reads {wildcards.sample}.{wildcards.lane}... #####"""
	input:
		"../results/trimmed_reads/{sample}.{lane}.trimmed_{type}.fastq.gz"
	output:
		"../results/fastqc_trimmed_reads/{sample}.{lane}.trimmed_{type}_fastqc.zip",
		"../results/fastqc_trimmed_reads/{sample}.{lane}.trimmed_{type}_fastqc.html"
	log:
		"../results/logs/trimming_reads/{sample}.{lane}.{type}.fastqc.log"
	threads:
		1
	shell:
		"fastqc -t {threads} -o ../results/fastqc_trimmed_reads {input} 2> {log}"
