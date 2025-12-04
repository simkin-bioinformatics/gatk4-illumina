I'm trying to determine the standard process for running GATK. As a first pass,
I asked gemini AI for 'standard practice'. Gemini came up with the following
basic steps:

For each individual sample:
1. Run bwa, with a special header line added
2. sort the sam file and output in bam format
3. Run the gatk markDuplicates step
4. Index the marked duplicate bam file
5. Run a baseRecalibrator to remove false positive variants - requires a file
of known variant sites (other sites treated as possible false positives)
6. Run gatk HaplotypeCaller

Combining samples:
1. run GenomicsDBImport to create a database from the individual samples
2. run GenotypeVCFs on the database to do joint variant calling on the database.

Option 1 for BQSR: could run the code below (either on my data or on the Pf8 database) to
pull out high confidence SNPs for BQSR.
gatk SelectVariants -R ref.fasta -V raw_variants.vcf.gz --select-type-to-include SNP -O raw_snps.vcf.gz
gatk VariantFiltration -R ref.fasta -V raw_snps.vcf.gz -O bootstrap.sites.vcf.gz --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" --filter-name "bootstrap_filter"

Option 2 for BQSR: 

Of course. If you're skipping BQSR, you'll use hard-filtering. This process involves two tools:

    SelectVariants: You first split your VCF into two files, one for SNPs and one for Indels, because they have different filtering criteria.

    VariantFiltration: You run this tool on each file to add "filter" tags (like LowQD or HighFS) to variants that fail your quality thresholds.

Here are the example commands, which assume your unfiltered VCF from GenotypeGVCFs is named cohort.vcf.gz.

Step 1: Split VCF into SNPs and Indels

First, you need to separate the SNPs and Indels from your main VCF file.

A. Select SNPs
Bash

gatk SelectVariants \
  -R reference.fasta \
  -V cohort.vcf.gz \
  --select-type-to-include SNP \
  -O cohort.snps.vcf.gz

B. Select Indels
Bash

gatk SelectVariants \
  -R reference.fasta \
  -V cohort.vcf.gz \
  --select-type-to-include INDEL \
  -O cohort.indels.vcf.gz

  Step 2: Apply Filters with VariantFiltration

Now you filter each file. These commands don't remove the variants; they add information to the FILTER column of the VCF. Any variant that passes all checks will have PASS in this column.

A. Filter SNPs

This command uses the standard GATK hard-filtering recommendations. Each filter-expression is paired with a filter-name that will be applied if the expression is true (i.e., if the variant fails the check).
Bash

gatk VariantFiltration \
  -R reference.fasta \
  -V cohort.snps.vcf.gz \
  -O cohort.snps.filtered.vcf.gz \
  --filter-expression "QD < 2.0" --filter-name "LowQD" \
  --filter-expression "QUAL < 30.0" --filter-name "LowQual" \
  --filter-expression "SOR > 3.0" --filter-name "HighSOR" \
  --filter-expression "FS > 60.0" --filter-name "HighFS" \
  --filter-expression "MQ < 40.0" --filter-name "LowMQ" \
  --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum" \
  --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum"

B. Filter Indels

The recommendations for Indels are different (and fewer).
Bash

gatk VariantFiltration \
  -R reference.fasta \
  -V cohort.indels.vcf.gz \
  -O cohort.indels.filtered.vcf.gz \
  --filter-expression "QD < 2.0" --filter-name "LowQD" \
  --filter-expression "QUAL < 30.0" --filter-name "LowQual" \
  --filter-expression "FS > 200.0" --filter-name "HighFS" \
  --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum"

Step 3 (Optional but Recommended): Create Final "PASS" VCF

After filtering, you have two VCFs (.filtered.vcf.gz) where bad variants are flagged. To create a final, clean VCF containing only the variants that passed all filters, you can merge the files and select only the "PASS" records.

A. Merge the Filtered VCFs
Bash

gatk MergeVcfs \
  -I cohort.snps.filtered.vcf.gz \
  -I cohort.indels.filtered.vcf.gz \
  -O cohort.all_variants.filtered.vcf.gz

B. Select Only "PASS" Variants This is your final, high-confidence variant file.
Bash

gatk SelectVariants \
  -V cohort.all_variants.filtered.vcf.gz \
  --exclude-filtered \
  -O cohort.final_pass_variants.vcf.gz

The file cohort.final_pass_variants.vcf.gz is the one you would use for your downstream analysis.

