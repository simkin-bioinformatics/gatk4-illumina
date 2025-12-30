'''
base quality score recalibration (BQSR) needs a list of high quality variant
calls as inputs. According to Chiyun Lee at the Sanger Institute, when he runs
this tool he uses as inputs files from here:
https://pf8-release.cog.sanger.ac.uk/vcf/index.html and filters them to list
only variants that are 'pass' and only keeps the following columns:
#CHROM  POS     ID          REF  ALT    QUAL  FILTER  INFO

This script re-does this BQSR parsing. Because the input vcf files are so large,
the lines from these files are parsed one by one.
'''

url_path=snakemake.params.url_path
unfiltered_vcf=snakemake.output.unfiltered_vcf
filtered_vcf=snakemake.output.filtered_vcf

def get_header(header):
	'''
	returns the indices associated with every column in the header
	'''
	h_dict={}
	for column_number, column in enumerate(header):
		h_dict[column]=column_number
	return h_dict

import gzip
import subprocess
targets=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']

subprocess.call(['wget', '-O', unfiltered_vcf, url_path])
parsing=False
output_file=open(filtered_vcf, 'w')
for line_number, line in enumerate(gzip.open(unfiltered_vcf, mode='rt')):
	if line_number%1000==0:
		print(unfiltered_vcf, line_number)
	split_line=line.strip().split('\t')
	if line.startswith('#CHROM'):
		h_dict=get_header(split_line)
		parsing=True
	if parsing:
		if split_line[h_dict['FILTER']]=='PASS':
			output_line=[]
			for target in targets:
				output_line.append(split_line[h_dict[target]])
			output_file.write('\t'.join(output_line)+'\n')
	else:
		print('line is', line)
		output_file.write(line)
