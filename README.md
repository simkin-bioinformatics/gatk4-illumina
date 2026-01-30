# gatk4-illumina

This pipeline is a general-purpose and simple pipeline for analyzing drug
resistance mutations in illumina datasets. It takes fastq files as inputs and
outputs amino acid tables that list the number of times each mutation was found
within each sample. These amino acid tables can then be input into a different
pipeline for visualizing the prevalences of drug resistance mutations:

https://github.com/simkin-bioinformatics/AA_table_visualization

## Prerequisite files
This pipeline assumes that a user has prepared the following:
1. A folder of demultiplexed paired end fastq files, with two fastq files per
   sample.
2. A fasta file containing a reference genome.
3. A gff file containing annotations of protein-coding genes within the genome.
4. A vcf file of high-confidence mutations. In our case, we created this dataset
   using this github repository:
   https://github.com/simkin-bioinformatics/sanger_known_sites_vcf.
5. A targets.tsv file of mutations that the user would like to "force" the
   program to make variant calls on. We provide an example targets.tsv file in
   this repo.

## Installation instructions

We are using a package manager called pixi for this project. You can obtain a
copy of pixi here:
https://pixi.prefix.dev/latest/installation/

After obtaining pixi, clone this github repository and change directory into the
cloned folder, then install the environment with the following command:
```
pixi install
```

## Running instructions

Open the 'config' folder and find the 'config.yaml' file. Open this file and
read the instructions carefully. Provide the paths for your input and output
files.

When you are finished editing the config.yaml file, run the pipeline with:
```
pixi run gatk4
```

## What this pipeline is doing:

### analyzing variants

The first step is variant calling with gatk, which can be roughly broken down
into the following steps:
1. bwa aligns Illumina reads against a reference genome.
2. gatk MarkDuplicates marks aligned bam files for potential duplicate reads.
3. gatk ApplyBQSR compares the observed mutations in the bam file against a list
   of known high confidence mutations, and uses this information to adjust the
   quality scores in the bam files.
4. gatk HaplotypeCaller calls mutations as "real" or "not real" and outputs them
   to a vcf file.
5. bcftools merge consolidates samples into a single vcf
6. gatk GenotypeGVCFs uses information from multiple samples to make a final
   decision about which mutations are real.

### annotating variants and generating an AA table
The second step is annotating the gatk-called variants with snpEff, which can
roughly be broken down into the following steps:
1. snpEff indexes your genome and gff files and builds a snpEff database
2. snpEff annotates your vcf file with protein-coding mutations.
3. A custom parser (written in Python) converts your vcf files into amino acid
   tables, with the number of times each mutation occurs within each sample.
