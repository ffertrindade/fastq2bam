def output_merge(wildcards):
	samples_multiple_lane = []
	samples_single_lane = []

	for sample in list(config["samples"]):

		lanes = []
		for lane in list(config["samples"][sample]):
			lanes.append(lane)
		
		if len(lanes) > 1:
			samples_multiple_lane.append(sample)
		else:
			samples_single_lane.append(sample)

	return samples_multiple_lane, samples_single_lane

samples_multiple_lane, samples_single_lane = output_merge(config)

if samples_multiple_lane == []:
	wildcard_constraints: samples_multiple_lane = "^$"
elif samples_single_lane == []:
	wildcard_constraints: samples_single_lane = "^$"
else:
	wildcard_constraints: samples_multiple_lane = '|'.join([x for x in samples_multiple_lane]),
	wildcard_constraints: samples_single_lane = '|'.join([x for x in samples_single_lane])

rule bwa_map_paired:
	message: """##### Running bwa mem for mapping paired reads of {wildcards.sample}.{wildcards.lane}... #####"""
	input:
		{config["reference"]},
		"results/trimmed_reads/{sample}.{lane}.trimmed_pair1.fastq.gz",
		"results/trimmed_reads/{sample}.{lane}.trimmed_pair2.fastq.gz"
	output:
		"results/mapped_reads/individual_raw_bam/{sample}.{lane}.paired.bam"
	conda:
		config["environment"]
	log:
		"results/logs/mapping/individual_raw_bam/{sample}.{lane}.bwa_map_paired.log"
	threads:
		config["threads"]
	params:
		rg=r"@RG\tID:{sample}_{lane}\tSM:{sample}",
		prefix="results/mapped_reads/individual_raw_bam/{sample}.{lane}.paired.bam"
	shell:
		"(bwa mem -R '{params.rg}' -t {threads} {input} | "
		"samtools view -Shb -F 4 - | "
		"samtools sort -@ {threads} -O bam -T results/mapped_reads/individual_raw_bam/{wildcards.sample}.{wildcards.lane}.pe -o {output} -) 2> {log}"

rule clip:
	message: """##### Clipping raw mapping for {wildcards.sample}.{wildcards.lane}... #####"""
	input:
		"results/mapped_reads/individual_raw_bam/{sample}.{lane}.paired.bam"
	output:
		temp("results/mapped_reads/individual_raw_bam/{sample}.{lane}.paired.clipped.bam")
	conda:
		config["environment"]
	log:
		"results/logs/mapping/{sample}.{lane}.clip.log"
	shell:
		"bam clipOverlap --in {input} --out {output} --stats --noPhoneHome > {log}"

rule merge_raw_bam:
	message: """##### Merging mapped lanes for {wildcards.samples_multiple_lane}... #####"""
	input:
		lambda wildcards: expand("results/mapped_reads/individual_raw_bam/{{samples_multiple_lane}}.{lane}.paired.clipped.bam", lane=config["samples"][wildcards.samples_multiple_lane]),
	output:
		temp("results/mapped_reads/individual_raw_bam/merged_lanes/{samples_multiple_lane}.bam"),
		temp("results/mapped_reads/individual_raw_bam/merged_lanes/{samples_multiple_lane}.bam.bai")
	conda:
		config["environment"]
	log:
		"results/logs/mapping/{samples_multiple_lane}.merge_raw_bam.log"
	threads:
		config["threads"]
	shell:
		"sambamba merge -p -t {threads} {output[0]} {input} 2> {log}"

rule markdup_merged:
	message: """##### Removing duplicated reads for sample {wildcards.samples_multiple_lane}... #####"""
	input:
		"results/mapped_reads/individual_raw_bam/merged_lanes/{samples_multiple_lane}.bam",
		"results/mapped_reads/individual_raw_bam/merged_lanes/{samples_multiple_lane}.bam.bai"
	output:
		temp("results/mapped_reads/{samples_multiple_lane}.markdup.bam"),
		temp("results/mapped_reads/{samples_multiple_lane}.markdup.bam.bai")
	conda:
		config["environment"]	
	log:
		"results/logs/mapping/{samples_multiple_lane}.markdup.log"
	threads:
		config["threads"]
	params:
		config["markdupPar"]
	shell:
		"sambamba markdup {params} -t {threads} --tmpdir results/mapped_reads/ {input[0]} {output[0]} 2> {log}"

rule markdup:
	message: """##### Removing duplicated reads for sample {wildcards.samples_single_lane}... #####"""
	input:
		lambda wildcards: expand("results/mapped_reads/individual_raw_bam/{{samples_single_lane}}.{lane}.paired.clipped.bam", lane=config["samples"][wildcards.samples_single_lane])
	output:
		temp("results/mapped_reads/{samples_single_lane}.markdup.bam"),
		temp("results/mapped_reads/{samples_single_lane}.markdup.bam.bai")
	conda:
		config["environment"]
	log:
		"results/logs/mapping/{samples_single_lane}.markdup.log"
	threads:
		config["threads"]
	params:
		config["markdupPar"]
	shell:
		"sambamba markdup {params} -t {threads} --tmpdir results/mapped_reads/ {input} {output[0]} 2> {log}"

