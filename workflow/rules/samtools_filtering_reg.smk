rule bam_filtering_reg:
	message: """##### Filtering bam regions {wildcards.sample}... #####"""
	input:
		bam="../results/mapped_reads/{sample}.sorted.merged.markdup.filtered.bam",
		bai="../results/mapped_reads/{sample}.sorted.merged.markdup.filtered.bam.bai",
		reg={config['regions']}
	output:
		"../results/mapped_reads/{sample}.sorted.merged.markdup.filtered.reg.bam"
	log:
		"../results/logs/mapping/{sample}.bam_filtering_reg.log"
	threads:
		config["threads"]
	shell:
		"samtools view -hb -@ {threads} -L {input.reg} -o {output} {input.bam} 2> {log}"
	
