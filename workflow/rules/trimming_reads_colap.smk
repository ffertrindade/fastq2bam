rule fastp:
	message: """##### Running fastp for {wildcards.sample}.{wildcards.lane}... #####"""
	input:
		pe1=lambda wildcards: config["samples"][wildcards.sample][wildcards.lane]["pe1"],
		pe2=lambda wildcards: config["samples"][wildcards.sample][wildcards.lane]["pe2"]
	output:
		pe1="results/trimmed_reads/{sample}.{lane}.trimmed_R1.fastq.gz",
		pe2="results/trimmed_reads/{sample}.{lane}.trimmed_R2.fastq.gz",
		colap="results/trimmed_reads/{sample}.{lane}.trimmed_colap.fastq.gz"
	conda:
		config["environment"]
	log:
		"results/logs/trimming_reads/{sample}.{lane}.fastp.log"
	threads:
		config["threads"]
	params:
		config["fastpPar"]
	shell:
		"fastp --in1 {input.pe1} --in2 {input.pe2} "
		"{params} --merge --thread {threads} "
		"--json results/logs/trimming_reads/{wildcards.sample}.{wildcards.lane}.fastp.json "
		"--html results/logs/trimming_reads/{wildcards.sample}.{wildcards.lane}.fastp.html "
		"--out1 {output.pe1} --out2 {output.pe2} --merged_out {output.colap} 2> {log}"