rule bam_filtering:
	message: """##### Filtering bam {wildcards.sample}... #####"""
	input:
		bam="results/mapped_reads/{sample}.sorted.merged.markdup.bam",
		bai="results/mapped_reads/{sample}.sorted.merged.markdup.bam.bai",
	output:
		"results/mapped_reads/{sample}.sorted.merged.markdup.filtered.bam"
	conda:
		config["environment"]
	log:
		"results/logs/mapping/{sample}.bam_filtering.log"
	threads:
		config["threads"]
	params:
		config["samtoolsFiltPar"]
	shell:
		"samtools view -@ {threads} {params} -o {output} {input.bam} 2> {log}"
	
rule bam_index_filt:
	message: """##### Indexing {wildcards.sample}.sorted.merged.markdup.filtered... #####"""
	input:
		"results/mapped_reads/{sample}.sorted.merged.markdup.filtered.bam"
	output:
		"results/mapped_reads/{sample}.sorted.merged.markdup.filtered.bam.bai"
	conda:
		config["environment"]
	log:
		"results/logs/mapping/{sample}.sorted.merged.markdup.filtered.bam_index.log"
	threads:
		config["threads"]
	shell:
		"sambamba index -p -t {threads} {input} {output} 2> {log}"