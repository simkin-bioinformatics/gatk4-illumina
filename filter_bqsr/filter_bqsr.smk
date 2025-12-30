'''
downloads unfiltered vcf files, filters on PASS variants, removes sample variant
calls, and sends output to a filtered file.
'''

def parse_filter_files(url_paths):
	files=[]
	for url_path in open(url_paths):
		url_path=url_path.strip()
		if url_path.endswith('.vcf.gz'):
			files.append(url_path.split('/')[-1][:-7])
	return files

url_paths_file='test_example.txt'
filter_files=parse_filter_files(url_paths_file)

rule all:
	input:
		final_files=expand('PASS_vcf_files/{filter_file}.vcf', filter_file=filter_files)

rule filter_bqsr:
	params:
		url_path='https://pf8-release.cog.sanger.ac.uk/vcf/{filter_file}.vcf.gz'
	output:
		unfiltered_vcf=temp('TEMP_vcf_files/{filter_file}.vcf'),
		filtered_vcf='PASS_vcf_files/{filter_file}.vcf'
	script:
		'scripts/filter_bqsr.py'