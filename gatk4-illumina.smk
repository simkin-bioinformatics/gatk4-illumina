import os

configfile: 'gatk4-illumina.yaml'
output_directory=config['output_directory']
ref_genome=config['ref_genome']
samples_file=config['samples_file']
sample_reads_folder=config['sample_reads_folder']
R1_suffix=config['R1_suffix']
R2_suffix=config['R2_suffix']
known_sites_vcf=config['known_sites_vcf']
known_sites_copy_path=os.path.join("temp_copied_input", os.path.basename(known_sites_vcf))
ref_genome_copy_path=os.path.join("temp_copied_input", os.path.basename(ref_genome))
#### get a list of samples from the samples file
samples = []
with open(samples_file, 'r') as f:
    for line in f:
        sample = line.strip()
        samples.append(sample)

# samples = ['Undetermined_S0']


rule all:
    input:
        # merged_samples = os.path.join(output_directory, "all_samples_merged.vcf.gz")
        cohort = os.path.join(output_directory, 'cohort.vcf.gz')

rule copy_genome_and_vcf:
    '''
    copy the reference genome into a local folder so that the indices generated
    here don't populate the user's original folder. The same needs to be done
    for the vcf file of known high confidence mutations.
    '''
    input:
        ref_genome = ref_genome,
        known_sites_vcf = known_sites_vcf
    output:
        copied_ref_genome = ref_genome_copy_path,
        copied_known_sites = known_sites_copy_path
    shell:
        '''
        cp {input.ref_genome} {output.copied_ref_genome}
        cp {input.known_sites_vcf} {output.copied_known_sites}
        '''
    
rule index_genome:
    input:
        copied_ref_genome = ref_genome_copy_path,
        copied_known_sites = known_sites_copy_path
    output:
        amb = ref_genome_copy_path + ".amb",
        ann = ref_genome_copy_path + ".ann",
        bwt = ref_genome_copy_path + ".bwt",
        pac = ref_genome_copy_path + ".pac",
        sa = ref_genome_copy_path + ".sa",
        fai = ref_genome_copy_path + ".fai",
        vcf_index = known_sites_copy_path + ".csi"
    shell:
        '''
        bwa index {input.copied_ref_genome}
        samtools faidx {input.copied_ref_genome}
        bcftools index {input.copied_known_sites}
        '''

rule create_gatk_dict:
    input:
        copied_ref_genome = ref_genome_copy_path
    output:
        gatk_dict = ref_genome_copy_path.replace(".fasta", ".dict"),
    # conda:
    #     "envs/gatk.yaml"
    shell:
        "gatk CreateSequenceDictionary -R {input.copied_ref_genome}"


rule bwa_mem:
    input:
        amb = ref_genome_copy_path + ".amb",
        ann = ref_genome_copy_path + ".ann",
        bwt = ref_genome_copy_path + ".bwt",
        pac = ref_genome_copy_path + ".pac",
        sa = ref_genome_copy_path + ".sa",
        fai = ref_genome_copy_path + ".fai",
        gatk_dict = ref_genome_copy_path.replace(".fasta", ".dict"),
        vcf_index = known_sites_copy_path + ".csi",
        copied_ref_genome = ref_genome_copy_path,
        read_1 = os.path.join(sample_reads_folder, "{sample}_"+R1_suffix),
        read_2 = os.path.join(sample_reads_folder, "{sample}_"+R2_suffix),

    output:
        sam = os.path.join(output_directory, "intermediate", "{sample}.sam")
    shell:
        # note that snakemake replaces \t with a literal tab so the \\t was necessitated below
        '''
        bwa mem \
            -o {output.sam} \
            -t 4 \
            -R "@RG\\tID:group1\\tSM:{wildcards.sample}\\tPL:ILLUMINA\\tLB:lib1" \
            {input.copied_ref_genome} {input.read_1} {input.read_2} \
        '''

rule samtools_sort:
    input:
        sam = os.path.join(output_directory, "intermediate", "{sample}.sam")
    output:
        sorted = os.path.join(output_directory, "intermediate", "{sample}.sorted.bam")
    shell:
        '''
        samtools sort -@ 4 -o {output.sorted} {input.sam}
        '''

rule mark_duplicates:
    input:
        sorted = os.path.join(output_directory, "intermediate", "{sample}.sorted.bam")
    output:
        marked_dups = os.path.join(output_directory, "intermediate", "{sample}.marked_dups.bam"),
        dup_metrics = os.path.join(output_directory, "intermediate", "{sample}.dup_metrics.txt")
    # conda:
    #     "envs/gatk.yaml"
    shell:
        '''
        gatk MarkDuplicates \
            -I {input.sorted} \
            -O {output.marked_dups} \
            -M {output.dup_metrics}
        '''

rule index_marked_bam:
    input:
        marked_dups = os.path.join(output_directory, "intermediate", "{sample}.marked_dups.bam"),
        dup_metrics = os.path.join(output_directory, "intermediate", "{sample}.dup_metrics.txt")
    output:
        marked_dups_index = os.path.join(output_directory, "intermediate", "{sample}.marked_dups.bam.bai")
    shell:
        '''
        samtools index {input.marked_dups}
        '''

rule index_feature_file:
    input:
        known_sites_vcf = known_sites_copy_path
    output:
        vcf_index = f"{known_sites_copy_path}.tbi"
    # conda:
    #     "envs/gatk.yaml"
    shell:
        '''
        gatk IndexFeatureFile -I {input.known_sites_vcf}
        '''

rule BQSR:
    input:
        marked_dups = os.path.join(output_directory, "intermediate", "{sample}.marked_dups.bam"),
        marked_dups_index = os.path.join(output_directory, "intermediate", "{sample}.marked_dups.bam.bai"),
        copied_ref_genome = ref_genome_copy_path,
        gatk_dict = ref_genome_copy_path.replace(".fasta", ".dict"),
        known_sites_vcf = known_sites_copy_path,
        vcf_index = f"{known_sites_copy_path}.tbi",
    output:
        recal_data_table = os.path.join(output_directory, "intermediate", "{sample}.recal_data.table")
    # conda:
    #     "envs/gatk.yaml"    
    shell:
        '''
        gatk BaseRecalibrator \
            -I {input.marked_dups} \
            -R {input.copied_ref_genome} \
            --known-sites {input.known_sites_vcf} \
            -O {output.recal_data_table}
        '''

rule apply_BQSR:
    input:
        marked_dups = os.path.join(output_directory, "intermediate", "{sample}.marked_dups.bam"),
        copied_ref_genome = ref_genome_copy_path,
        recal_data_table = os.path.join(output_directory, "intermediate", "{sample}.recal_data.table")
    output:
        analysis_ready_bam = os.path.join(output_directory, "intermediate", "{sample}.analysis_ready.bam")
    # conda:
    #     "envs/gatk.yaml"
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
        analysis_ready_bam = os.path.join(output_directory, "intermediate", "{sample}.analysis_ready.bam")
    output:
        analysis_bam_index = os.path.join(output_directory, "intermediate", "{sample}.analysis_ready.bam.bai")
    shell:
        "samtools index {input.analysis_ready_bam}"

rule haplotype_caller:
    input:
        analysis_ready_bam = os.path.join(output_directory, "intermediate", "{sample}.analysis_ready.bam"),
        analysis_bam_index = os.path.join(output_directory, "intermediate", "{sample}.analysis_ready.bam.bai"),
        ref_genome = ref_genome_copy_path
    output:
        called_haplotypes = os.path.join(output_directory, "intermediate", "{sample}.g.vcf.gz")
    # conda:
    #     "envs/gatk.yaml"
    shell:
        '''
        gatk HaplotypeCaller \
            -I {input.analysis_ready_bam} \
            -R {input.ref_genome} \
            -ERC GVCF \
            -O {output.called_haplotypes}
        '''

rule consolidate gVCFs:
# there is an alternate method using gatk GenomicsDBImport but here we will use bcftools merge
# the other method is probably faster
    input:
        called_haplotypes = expand(os.path.join(output_directory, "intermediate", "{sample}.g.vcf.gz"), sample = samples)
    output:
        called_merged_haplotypes = os.path.join(output_directory, "intermediate", "called_merged_haplotypes.vcf.gz")
    shell:
        # the --force-samples was used because each sample in every sample file in the test data had
        # the same sample1 name
        '''
        bcftools merge {input.called_haplotypes} \
            -O z -o {output.called_merged_haplotypes} \
            --force-single
        '''

rule index_consolidated_vcf:
    input:
        called_merged_haplotypes = os.path.join(output_directory, "intermediate", "called_merged_haplotypes.vcf.gz")
    output:
        merged_vcf_index = f"{os.path.join(output_directory, "intermediate", "called_merged_haplotypes.vcf.gz")}.tbi"
    # conda:
    #     "envs/gatk.yaml"
    shell:
        '''
        gatk IndexFeatureFile -I {input.called_merged_haplotypes}
        '''


rule joint_genotyping:
    input:
        ref_genome = ref_genome_copy_path,
        called_merged_haplotypes = os.path.join(output_directory, "intermediate", "called_merged_haplotypes.vcf.gz"),
        merged_vcf_index = f"{os.path.join(output_directory, "intermediate", "called_merged_haplotypes.vcf.gz")}.tbi"

    output:
        cohort = os.path.join(output_directory, 'cohort.vcf.gz')
    # conda:
    #     "envs/gatk.yaml"
    shell:
        '''
        gatk GenotypeGVCFs \
            -R {input.ref_genome} \
            -V {input.called_merged_haplotypes} \
            -O {output.cohort}
        '''
