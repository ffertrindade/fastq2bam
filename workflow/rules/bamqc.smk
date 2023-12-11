rule qualimap:
	message: """##### Running qualimap for bam {wildcards.sample}.{wildcards.stage}... #####"""
	input:
		"results/mapped_reads/{sample}.{stage}.bam"
	output:
		pdf="results/qualimap/{sample}.{stage}/{sample}.{stage}.pdf",
		txt="results/qualimap/{sample}.{stage}/{sample}.{stage}.txt"
	conda:
		config["environment"]
	log:
		"results/logs/mapping/{sample}.{stage}.bam_qc.log"
	threads:
		config["threads"]
	shell:
		"""
		(qualimap bamqc -c -nt {threads} -bam {input} -outfile {wildcards.sample}.{wildcards.stage}.pdf -outdir results/qualimap/{wildcards.sample}.{wildcards.stage} &&
		mv results/qualimap/{wildcards.sample}.{wildcards.stage}/genome_results.txt {output.txt}) > {log}
		"""

rule flagstat:
	message: """##### Running samtools flagstat for bam {wildcards.sample}.{wildcards.stage}... #####"""
	input:
		bam="results/mapped_reads/{sample}.{stage}.bam",
		bai="results/mapped_reads/{sample}.{stage}.bam.bai"
	output:
		"results/mapped_reads/flagstat/{sample}.{stage}.flagstat.txt"
	conda:
		config["environment"]
	log:
		"results/logs/mapping/flagstat/{sample}.{stage}.flagstat.log"
	shell:
		"samtools flagstat {input.bam} > {output} 2> {log}"

rule plotCoverage:
	message: """##### Running plotCoverage for bams {wildcards.stage}... #####"""
	input:
		bam=expand("results/mapped_reads/{sample}.{{stage}}.bam", sample=config["samples"]),
		bai=expand("results/mapped_reads/{sample}.{{stage}}.bam.bai", sample=config["samples"])
	output:
		png="results/deeptools/plotCoverage_{stage}.png",
		txt="results/deeptools/plotCoverage_{stage}.txt"
	conda:
		config["environment"]
	log:
		"results/logs/mapping/{stage}.plotCoverage.log"
	threads:
		config["threadsPlotCoverage"]
	shell:
		"plotCoverage -b {input.bam} -p {threads} -o {output.png} > {output.txt} 2> {log}"
