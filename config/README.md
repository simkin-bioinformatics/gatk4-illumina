# filter bqsr

This script is for obtaining a list of high quality known variable sites in the
falciparum genome. It uses falciparum data from the Sanger institute (8,000
variant called samples I think). Because these vcf files (which are separated
into chromosomes) are huge and we don't need the variant calls of every sample
(only a list of genomic positions where variants are found), the first step is
to download one file at a time, throw out any genomic positions that don't
"PASS" the criteria set by the sanger institute, and remove all of the
individual sample variant calls, then delete the original large unfiltered file.

Filtered chromosome files are then gzipped, indexed, and combined together to
create a single genome-wide file of known high quality variants.

After obtaining micromamba, the dependencies (tabix, snakemake, python, and
bcftools) can be met with conda create -f environment.yml

The environment can be activated with micromamba activate gatk_known_sites

The workflow can be executed with snakemake -s filter_bqsr.smk --cores 2

If some time passes before you need to execute this pipeline, you may need to
update the links to the original vcf files, which are currently stored in
sanger_vcf_file_links.txt