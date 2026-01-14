import os


ref_genome = "gatk4-illumina_test_data/genome/PlasmoDB-62_Pfalciparum3D7_Genome.fasta"
output_directory = "gatk4-illumina_test_data/output2"
output_name = "Undetermined"

read_1 = "gatk4-illumina_test_data/fastqs/Undetermined_S0_R1_001.fastq.gz"
read_2 = "gatk4-illumina_test_data/fastqs/Undetermined_S0_R2_001.fastq.gz"

merged_vcf='gatk4-illumina_test_data/PASS_vcf_with_headers/merged.vcf.gz'

copied_ref_genome = os.path.join(output_directory, "copied_data", os.path.basename(ref_genome))
copied_vcf = os.path.join(output_directory, "copied_data", os.path.basename(merged_vcf))
print(copied_ref_genome)
rule all:
    input:
        cohort = os.path.join(output_directory, 'cohort.vcf.gz')

rule copy_genome_and_vcf:
    input:
        ref_genome = ref_genome,
        merged_vcf = merged_vcf
    output:
        copied_ref_genome = copied_ref_genome,
        copied_vcf = copied_vcf,
    shell:
        '''
        cp {input.ref_genome} {output.copied_ref_genome}
        cp {input.merged_vcf} {output.copied_vcf}
        '''
    
rule index_genome:
    input:
        copied_ref_genome = copied_ref_genome,
        copied_vcf = copied_vcf
    output:
        amb = copied_ref_genome + ".amb",
        ann = copied_ref_genome + ".ann",
        bwt = copied_ref_genome + ".bwt",
        pac = copied_ref_genome + ".pac",
        sa = copied_ref_genome + ".sa",
        fai = copied_ref_genome + ".fai",
        gatk_dict = copied_ref_genome.replace(".fasta", ".dict"),
        vcf_index = copied_vcf + ".csi"
    shell:
        '''
        bwa index {input.copied_ref_genome} && \
        samtools faidx {input.copied_ref_genome} && \
        gatk CreateSequenceDictionary -R {input.copied_ref_genome}
        bcftools index {input.copied_vcf}
        '''

rule bwa_mem:
    input:
        amb = copied_ref_genome + ".amb",
        ann = copied_ref_genome + ".ann",
        bwt = copied_ref_genome + ".bwt",
        pac = copied_ref_genome + ".pac",
        sa = copied_ref_genome + ".sa",
        fai = copied_ref_genome + ".fai",
        gatk_dict = copied_ref_genome.replace(".fasta", ".dict"),
        vcf_index = copied_vcf + ".csi",
        copied_ref_genome = copied_ref_genome,
        read_1 = read_1,
        read_2 = read_2,

    output:
        sam = os.path.join(output_directory, f"{output_name}.sam")
    shell:
        # note that snakemake replaces \t with a literal tab so the \\t was necessitated below
        '''
        bwa mem \
            -o {output.sam} \
            -t 4 \
            -R "@RG\\tID:group1\\tSM:sample1\\tPL:ILLUMINA\\tLB:lib1" \
            {input.copied_ref_genome} {input.read_1} {input.read_2} \
        '''

rule samtools_sort:
    input:
        sam = os.path.join(output_directory, f"{output_name}.sam")
    output:
        sorted = os.path.join(output_directory, f"{output_name}.sorted.bam")
    shell:
        '''
        samtools sort -@ 4 -o {output.sorted} {input.sam}
        '''

rule mark_duplicates:
    input:
        sorted = os.path.join(output_directory, f"{output_name}.sorted.bam")
    output:
        marked_dups = os.path.join(output_directory, f"{output_name}.marked_dups.bam"),
        dup_metrics = os.path.join(output_directory, f"{output_name}.dup_metrics.txt")
    shell:
        '''
        gatk MarkDuplicates \
            -I {input.sorted} \
            -O {output.marked_dups} \
            -M {output.dup_metrics}
        '''

rule index_marked_bam:
    input:
        marked_dups = os.path.join(output_directory, f"{output_name}.marked_dups.bam"),
        dup_metrics = os.path.join(output_directory, f"{output_name}.dup_metrics.txt")
    output:
        marked_dups_index = os.path.join(output_directory, f"{output_name}.marked_dups.bam.bai")
    shell:
        '''
        samtools index {input.marked_dups}
        '''

rule index_feature_file:
    input:
        merged_vcf = copied_vcf
    output:
        vcf_index = f"{copied_vcf}.tbi"
    shell:
        '''
        gatk IndexFeatureFile -I {input.merged_vcf}
        '''

rule BQSR:
    input:
        marked_dups = os.path.join(output_directory, f"{output_name}.marked_dups.bam"),
        marked_dups_index = os.path.join(output_directory, f"{output_name}.marked_dups.bam.bai"),
        copied_ref_genome = copied_ref_genome,
        gatk_dict = copied_ref_genome.replace(".fasta", ".dict"),
        copied_vcf = copied_vcf,
        vcf_index = f"{copied_vcf}.tbi",
    output:
        recal_data_table = os.path.join(output_directory, f"{output_name}.recal_data.table")
    shell:
        '''
        gatk BaseRecalibrator \
            -I {input.marked_dups} \
            -R {input.copied_ref_genome} \
            --known-sites {input.copied_vcf} \
            -O {output.recal_data_table}
        '''

rule apply_BQSR:
    input:
        marked_dups = os.path.join(output_directory, f"{output_name}.marked_dups.bam"),
        copied_ref_genome = copied_ref_genome,
        recal_data_table = os.path.join(output_directory, f"{output_name}.recal_data.table")
    output:
        analysis_ready_bam = os.path.join(output_directory, f"{output_name}.analysis_ready.bam")
    shell:
        '''
        gatk ApplyBQSR \
        -I {input.marked_dups} \
        -R {input.copied_ref_genome} \
        --bqsr-recal-file {input.recal_data_table} \
        -O {output.analysis_ready_bam}
        '''

rule index_BQSR_bam:
    input:
        analysis_ready_bam = os.path.join(output_directory, f"{output_name}.analysis_ready.bam")
    output:
        analysis_bam_index = os.path.join(output_directory, f"{output_name}.analysis_ready.bam.bai")
    shell:
        "samtools index {input.analysis_ready_bam}"

rule haplotype_caller:
    input:
        analysis_ready_bam = os.path.join(output_directory, f"{output_name}.analysis_ready.bam"),
        analysis_bam_index = os.path.join(output_directory, f"{output_name}.analysis_ready.bam.bai"),
        ref_genome = copied_ref_genome
    output:
        called_haplotypes = os.path.join(output_directory, f"{output_name}.g.vcf.gz")
    shell:
        '''
        gatk HaplotypeCaller \
            -I {input.analysis_ready_bam} \
            -R {input.ref_genome} \
            -ERC GVCF \
            -O {output.called_haplotypes}
        '''

# #### Step 4: Consolidate gVCFs
# rule consolidate_gVCFs:
# #   * **Note:** This tool works on intervals. For a basic run, you can specify your main chromosomes 
# #   (e.g., `-L chr1 -L chr2`). For whole-genome data, it's often run per chromosome or chromosome group.
#     input:
#         called_haplotypes = os.path.join(output_directory, f"{output_name}.g.vcf.gz")
#     output:
#         touch(os.path.join(output_directory, 'consolidate_done.txt'))
#     params:
#         genomicsdb_workspace_path = os.path.join(output_directory, "my_database")
#     shell:
#         '''
#         gatk GenomicsDBImport \
#             --genomicsdb-workspace-path {params.genomicsdb_workspace_path} \
#             -L Pf3D7_01_v3 \
#             -L Pf3D7_02_v3 \
#             -L Pf3D7_03_v3 \
#             -L Pf3D7_04_v3 \
#             -L Pf3D7_05_v3 \
#             -L Pf3D7_06_v3 \
#             -L Pf3D7_07_v3 \
#             -L Pf3D7_08_v3 \
#             -L Pf3D7_09_v3 \
#             -L Pf3D7_10_v3 \
#             -L Pf3D7_11_v3 \
#             -L Pf3D7_12_v3 \
#             -L Pf3D7_13_v3 \
#             -L Pf3D7_14_v3 \
#             -V {input.called_haplotypes}
#         '''
#             # -V sample2.g.vcf.gz \
#             # -V sample3.g.vcf.gz


rule joint_genotyping:
    input:
        ref_genome = copied_ref_genome,
        called_haplotypes = os.path.join(output_directory, f"{output_name}.g.vcf.gz")
    output:
        cohort = os.path.join(output_directory, 'cohort.vcf.gz')
    shell:
        '''
        gatk GenotypeGVCFs \
            -R {input.ref_genome} \
            -V {input.called_haplotypes} \
            -O {output.cohort}
        '''
#### Step 5: Joint Genotyping

# This is the final step. It takes the combined database from `GenomicsDBImport` and runs the genotyping algorithm across all samples simultaneously. This "joint" analysis uses information from all samples to make more accurate calls for each one.

#   * **Tool:** `gatk GenotypeGVCFs`
#   * **Basic Command:**
#     ```bash

#     ```

# You will now have a single file, **`cohort.vcf.gz`**, containing the variants for your entire set of samples. The next step, which is more advanced, would be to filter this VCF using `VariantRecalibrator` (VQSR), but the steps above will give you a complete, un-filtered callset.