'''
base quality score recalibration (BQSR) needs a list of high quality variant
calls as inputs. According to Chiyun Lee at the Sanger Institute, when he runs
this tool he uses as inputs files from here:
https://pf8-release.cog.sanger.ac.uk/vcf/index.html and filters them to list
only variants that are 'pass' and only keeps the following columns:
#CHROM  POS     ID          REF  ALT    QUAL  FILTER  INFO

This script re-does this BQSR parsing. Because the input vcf files are so large,
the lines from these files are parsed one by one.

This version of the script attempts to fix a parsing error where the header line
that starts #CHROM was not added to the parsed vcf file - the edited script has
not been tested because re-generating the parsed VCF files would take several
days. (attempted fix is on line 44).
'''

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
for line in open('sanger_vcf_file_links.txt'):
	html_path=line.strip()
	file_path=html_path.split('/')[-1]
#	subprocess.call(['wget', html_path])
	if file_path.endswith('.vcf.gz'):
		parsing=False
		output_file=open('parsed_files/'+file_path[:-3], 'w')
		for line_number, line in enumerate(gzip.open(file_path, mode='rt')):
			if line_number%1000==0:
				print(line_number)
			split_line=line.strip().split('\t')
			if line.startswith('#CHROM'):
				h_dict=get_header(split_line)
				parsing=True
				output_file.write('\t'.join(targets)+'\n')
			if parsing:
#				print(split_line[h_dict['FILTER']])
				if split_line[h_dict['FILTER']]=='PASS':
					output_line=[]
					for target in targets:
						output_line.append(split_line[h_dict[target]])
					output_file.write('\t'.join(output_line)+'\n')
			else:
				output_file.write(line)
				#print(split_line)