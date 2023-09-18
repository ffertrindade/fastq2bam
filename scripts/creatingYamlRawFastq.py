#!/usr/bin/python3.6
## Build a list of fastq files in yaml format
## All fastq files must be in the same directory
## Used to create samples input files of the config file for snakemake workflow fastq2bam

# Libraries and arguments
import sys
import os
import glob
from py_linq import Enumerable
import yaml

my_args = sys.argv
if len(my_args) != 3:
        print("Usage: creatingYamlFastq.py list_id_samples paths_with_fastq")
        exit()
else:
        sample_list_file = my_args[1]
        fastq_path = bam1 = my_args[2]

        if os.path.isfile(sample_list_file):
                pass
        else:
                print("File doesnt exist!")
                exit()

# List of fastq files from path
def get_samples_list(path):
	file = open(path, "r")
	samples = [line.rstrip() for line in file]
	file.close()
	return samples

## List of fastq from target samples
def filter_files_by_sample(path, samples):
	files = glob.glob(path + '/*')
	sample_files = []
	for s in samples:
		enumerable = Enumerable(files)
		filtered_files = enumerable.where(lambda x: x.find(s) != -1)
		sample_files.extend(filtered_files.to_list())

	return sample_files

## Create a dictionary to be read into a yaml format
def create_dict_for_yaml(files):
	sample = {}
	sample['samples'] = {}
	for f in files:
		file_name = f.split("/")[-1]
		sample_name = file_name.split("_")[1]
		dataset_name = file_name.split("_")[0]

		if sample_name not in sample['samples']:
			sample['samples'][sample_name] = {}

		if dataset_name not in sample['samples'][sample_name]:
			sample['samples'][sample_name][dataset_name] = {}

		if f.endswith("_R1.fastq.gz"):
			sample['samples'][sample_name][dataset_name]["pe1"] = f
		if f.endswith("_R2.fastq.gz"):
			sample['samples'][sample_name][dataset_name]["pe2"] = f

	return sample

## Running using the given parameters
samples = get_samples_list(sample_list_file)
sample_files = filter_files_by_sample(fastq_path, samples)
samples_dict = create_dict_for_yaml(sample_files)
print(yaml.dump(samples_dict))

