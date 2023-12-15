rule bam_filtering_reg:
	message: """##### Filtering bam regions {wildcards.sample}... #####"""
	input:
		bam="results/mapped_reads/{sample}.markdup.filtered.bam",
		bai="results/mapped_reads/{sample}.markdup.filtered.bam.bai",
		reg={config['regions']}
	output:
		"results/mapped_reads/{sample}.markdup.filtered.reg.bam"
	conda:
		config["environment"]
	log:
		"results/logs/mapping/{sample}.bam_filtering_reg.log"
	threads:
		config["threads"]
	shell:
		"samtools view -hb -@ {threads} -L {input.reg} -o {output} {input.bam} 2> {log}"

rule bam_index_reg:
	message: """##### Indexing {wildcards.sample}.markdup.filtered.reg... #####"""
	input:
		"results/mapped_reads/{sample}.markdup.filtered.reg.bam"
	output:
		"results/mapped_reads/{sample}.markdup.filtered.reg.bam.bai"
	conda:
		config["environment"]
	log:
		"results/logs/mapping/{sample}.markdup.filtered.reg.bam_index.log"
	threads:
		config["threads"]
	shell:
		"sambamba index -p -t {threads} {input} {output} 2> {log}"
