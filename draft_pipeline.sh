# Variables for clarity
REF_GENOME="/home/charlie/git/gatk4-illumina/input/genome/PlasmoDB-62_Pfalciparum3D7_Genome.fasta"
READ1="/home/charlie/git/gatk4-illumina/input/fastqs/Undetermined_S0_R1_001.fastq.gz"
READ2="/home/charlie/git/gatk4-illumina/input/fastqs/Undetermined_S0_R2_001.fastq.gz"
OUTPUTFOLDER="output"
OUTPUTNAME='Undetermined'
MERGEDVCF='/home/charlie/git/gatk4-illumina/output/PASS_vcf_with_headers/merged.vcf.gz'

# OUTPUT_BAM="output/Undetermined.sorted.bam"

# # Index the genome for bwa
# # creates index files in same folder as genome
# bwa index $REF_GENOME


# # The -R option adds the read group. Adjust SM (Sample) and LB (Library) accordingly.
# # We pipe (|) the output of bwa mem directly to samtools sort to save space and time.
# bwa mem -o $OUTPUTFOLDER/$OUTPUTNAME.sam -t 4 -R "@RG\tID:group1\tSM:sample1\tPL:ILLUMINA\tLB:lib1" $REF_GENOME $READ1 $READ2

# samtools sort -@ 4 -o $OUTPUTFOLDER/$OUTPUTNAME.sorted.bam $OUTPUTFOLDER/$OUTPUTNAME.sam

# # Mark Duplicates
# gatk MarkDuplicates \
#  -I $OUTPUTFOLDER/$OUTPUTNAME.sorted.bam \
#  -O $OUTPUTFOLDER/$OUTPUTNAME.marked_dups.bam \
#  -M $OUTPUTFOLDER/$OUTPUTNAME.dup_metrics.txt

# # Index the new BAM file
# samtools index $OUTPUTFOLDER/$OUTPUTNAME.marked_dups.bam

# # Create the FASTA Index
# samtools faidx $REF_GENOME

# # Create GATK sequence dictionary
# gatk CreateSequenceDictionary \
#      -R $REF_GENOME

# # create index for known-sites file
# gatk IndexFeatureFile -I $MERGEDVCF

# gatk BaseRecalibrator \
#  -I $OUTPUTFOLDER/$OUTPUTNAME.marked_dups.bam \
#  -R $REF_GENOME \
#  --known-sites $MERGEDVCF \
#  -O $OUTPUTFOLDER/sample1.recal_data.table

gatk HaplotypeCaller \
  -I $OUTPUTFOLDER/$OUTPUTNAME.marked_dups.bam \
  -R $REF_GENOME \
  -ERC GVCF \
  -O $OUTPUTFOLDER/$OUTPUTNAME.g.vcf.gz
