'''
downloads unfiltered vcf files, filters on PASS variants, removes sample variant
calls, and sends output to a filtered file, one per chromosome. These filtered
chromosomes are then merged into a single file of known sites.
'''

def parse_filter_files(url_paths):
	files=[]
	for url_path in open(url_paths):
		url_path=url_path.strip()
		if url_path.endswith('.vcf.gz'):
			files.append(url_path.split('/')[-1][:-12])
	return files

chroms=parse_filter_files(url_paths_file)

rule all:
	input:
		known_sites='known_sites/known_sites.vcf.gz'

rule filter_bqsr:
	'''
	downloads files one at a time, filters them, and then deletes the original
	unfiltered file. This keeps the hard drive from filling up.
	'''
	params:
		url_path='https://pf8-release.cog.sanger.ac.uk/vcf/{chrom}.filt.vcf.gz'
	output:
		unfiltered_vcf=temp('TEMP_vcf_files/{chrom}.filt.vcf'),
		filtered_vcf=temp('PASS_vcf_files/{chrom}.filt.vcf')
	script:
		'scripts/filter_bqsr.py'

rule zip_files:
	'''
	because the original filtered vcf files are unzipped raw text, the filtered
	vcf files need to be zipped up.
	'''
	input:
		pass_file='PASS_vcf_files/{chrom}.filt.vcf'
	output:
		zipped_file=temp('zipped_files/{chrom}.filt.vcf.gz')
	shell:
		'bgzip -c {input.pass_file} > {output.zipped_file}'

rule index_files:
	'''
	indexes the zipped vcf files using bcftools
	'''
	input:
		zipped_file='zipped_files/{chrom}.filt.vcf.gz'
	output:
		indexed_file=temp('zipped_files/{chrom}.filt.vcf.gz.csi')
	shell:
		'bcftools index {input.zipped_file}'

rule merge_files:
	'''
	merges the zipped vcf files together into a single multi-chromosome vcf
	file.
	'''
	input:
		indexed_files=expand('zipped_files/{chrom}.filt.vcf.gz.csi', chrom=chroms),
		zipped_files=expand('zipped_files/{chrom}.filt.vcf.gz', chrom=chroms)
	output:
		known_sites='known_sites/known_sites.vcf.gz'
	shell:
		'bcftools merge {input.zipped_files} -o {output.known_sites} -O z'