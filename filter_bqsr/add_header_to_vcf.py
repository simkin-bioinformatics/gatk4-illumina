#!/usr/bin/env python3

import os
import bgzip
import gzip
import subprocess
from joblib import Parallel, delayed

header_line = "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"

input_file_directory = '/home/charlie/git/gatk4-illumina/input/PASS_vcf'
output_folder = '/home/charlie/git/gatk4-illumina/output/PASS_vcf_with_headers'
merged_vcf_name = 'merged'

os.makedirs(output_folder, exist_ok = True)

def add_header(input_file_name):
	with gzip.open(os.path.join(input_file_directory, input_file_name), 'rt') as input_file:
		with open(os.path.join(output_folder, input_file_name).strip('.gz'), 'w') as output_file:
			insert_header = 0
			samples = 0
			for line in input_file:
				if "##" in line:
					output_file.write(line)
				if "##" not in line:
					insert_header += 1
					samples += 1
				if insert_header == 1:
					output_file.write(header_line +'\n')
					insert_header += 1
				if insert_header > 1:
					output_file.write(line)
			if samples == 0:
				output_file.write(header_line +'\n')

def bgzip(input_file):
	subprocess.run(['bgzip', os.path.join(output_folder, input_file)])

def bcf_index(input_file):
	subprocess.run(['bcftools', 'index', os.path.join(output_folder, input_file)])

def vcf_merge(input_list):
	subprocess.run(
		['bcftools', 'merge'] + 
		input_list + 
		['-o', os.path.join(output_folder, f"{merged_vcf_name}.vcf.gz"), '-O', 'z']
	)

Parallel(n_jobs=-1)(delayed(add_header)(i) for i in os.listdir(input_file_directory))
Parallel(n_jobs=-1)(delayed(bgzip)(i) for i in os.listdir(output_folder))
Parallel(n_jobs=-1)(delayed(bcf_index)(i) for i in os.listdir(output_folder))

input_list = [os.path.join(output_folder, item) for item in os.listdir(output_folder) if '.csi' not in item]
vcf_merge(input_list)
subprocess.run(['bcftools', 'index', os.path.join(output_folder, f"{merged_vcf_name}.vcf.gz")])



