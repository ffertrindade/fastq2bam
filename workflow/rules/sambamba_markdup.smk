rule markdup:
	message: """##### Removing duplicated reads for sample {wildcards.sample}... #####"""
	input:
		"results/mapped_reads/{sample}.sorted.merged.bam"
	output:
		temp("results/mapped_reads/{sample}.sorted.merged.markdup.bam")
	conda:
		config["environment"]	
	log:
		"results/logs/mapping/{sample}.markdup.log"
	threads:
		config["threads"]
	params:
		config["markdupPar"]
	shell:
		"sambamba markdup {params} -t {threads} --tmpdir results/mapped_reads/ {input} {output} 2> {log}"

rule bam_index_markdup:
	message: """##### Indexing {wildcards.sample}.sorted.merged.markdup... #####"""
	input:
		"results/mapped_reads/{sample}.sorted.merged.markdup.bam"
	output:
		temp("results/mapped_reads/{sample}.sorted.merged.markdup.bam.bai")
	conda:
		config["environment"]
	log:
		"results/logs/mapping/{sample}.sorted.merged.markdup.bam_index.log"
	threads:
		config["threads"]
	shell:
		"sambamba index -p -t {threads} {input} {output} 2> {log}"