# Variables for clarity
REF_GENOME="/home/alfred/gatk4-illumina/Pf_3D7_genome/PlasmoDB-62_Pfalciparum3D7_Genome.fasta"
READ1=/home/alfred/gatk4-illumina/251006_miseq/fastq/Undetermined_S0_R1_001.fastq.gz
READ2=/home/alfred/gatk4-illumina/251006_miseq/fastq/Undetermined_S0_R2_001.fastq.gz
OUTPUT_BAM="Undetermined.sorted.bam"

# The -R option adds the read group. Adjust SM (Sample) and LB (Library) accordingly.
# We pipe (|) the output of bwa mem directly to samtools sort to save space and time.
#bwa mem -o Undetermined.sam -t 4 -R "@RG\tID:group1\tSM:sample1\tPL:ILLUMINA\tLB:lib1" $REF_GENOME $READ1 $READ2

#samtools sort -@ 4 -o Undetermined.sorted.bam Undetermined.sam

#gatk MarkDuplicates \
#  -I Undetermined.sorted.bam \
#  -O Undetermined.marked_dups.bam \
#  -M Undetermined.dup_metrics.txt

# Index the new BAM file
samtools index Undetermined.marked_dups.bam

#gatk BaseRecalibrator \
#  -I Undetermined.marked_dups.bam \
#  -R $REF_GENOME \
#  --known-sites known_sites.vcf.gz \
#  -O sample1.recal_data.table

gatk HaplotypeCaller \
  -I Undetermined.marked_dups.bam \
  -R $REF_GENOME \
  -ERC GVCF \
  -O Undetermined.g.vcf.gz
