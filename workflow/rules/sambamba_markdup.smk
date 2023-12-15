rule markdup:
	message: """##### Removing duplicated reads for sample {wildcards.sample}... #####"""
	input:
		"results/mapped_reads/{sample}.bam",
		"results/mapped_reads/{sample}.bam.bai"
	output:
		temp("results/mapped_reads/{sample}.markdup.bam"),
		temp("results/mapped_reads/{sample}.markdup.bam.bai")
	conda:
		config["environment"]	
	log:
		"results/logs/mapping/{sample}.markdup.log"
	threads:
		config["threads"]
	params:
		config["markdupPar"]
	shell:
		"sambamba markdup {params} -t {threads} --tmpdir results/mapped_reads/ {input[0]} {output[0]} 2> {log}"
