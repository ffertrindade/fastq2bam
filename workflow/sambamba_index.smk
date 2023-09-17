rule bam_index:
	message: """##### Indexing {wildcards.sample}.{wildcards.stage}... #####"""
	input:
		"../results/mapped_reads/{sample}.{stage}.bam"
	output:
		"../results/mapped_reads/{sample}.{stage}.bam.bai"
	log:
		"../results/logs/mapping/{sample}.{stage}.bam_index.log"
	threads:
		config["threads"]
	shell:
		"sambamba index -p -t {threads} {input} {output} 2> {log}"
