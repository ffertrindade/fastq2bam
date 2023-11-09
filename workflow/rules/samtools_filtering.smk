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
	
