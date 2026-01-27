rule bwa_mem:
    input:
        ref_genome = copied_ref_genome,
        bwa_indices = multiext(copied_ref_genome, ".amb", ".ann", ".bwt", ".pac", ".sa"),
        sam_index = f"{copied_ref_genome}.fai",
        read_1 = Path(config['sample_reads_folder']) / str("{sample}"+config['R1_suffix']),
        read_2 = Path(config['sample_reads_folder']) / str("{sample}" + config['R2_suffix']),
    output:
        bam = results / "aligned_reads" / "{sample}.bam"
    shell:
        "(bwa mem "
        "-R '@RG\\tID:group1\\tSM:{wildcards.sample}\\tPL:ILLUMINA\\tLB:lib1' "
        "{input.ref_genome} {input.read_1} {input.read_2} | "
        "samtools sort -o {output.bam} -)"

rule mark_duplicates:
    input:
        bam = results / "aligned_reads" / "{sample}.bam"
    output:
        marked_dups = results / "aligned_reads" / "{sample}.marked_dups.bam",
        dup_metrics = results / "aligned_reads" / "{sample}.dup_metrics.txt",
        marked_dups_index = results / "aligned_reads" / "{sample}.marked_dups.bai"
    threads: 2
    shell:
        '''
        gatk MarkDuplicates \
            -I {input.bam} \
            -O {output.marked_dups} \
            -M {output.dup_metrics} \
            --CREATE_INDEX true
       '''

rule BQSR:
    input:
        ref_genome = copied_ref_genome,
        gatk_ref_dict = copied_ref_genome.replace(".fasta", ".dict"),
        known_sites_vcf = copied_known_sites_vcf,
        known_sites_index = f"{copied_known_sites_vcf}.tbi",
        marked_dups = results / "aligned_reads" / "{sample}.marked_dups.bam",
        marked_dups_index = results / "aligned_reads" / "{sample}.marked_dups.bai"
    output:
        recal_data_table = results / "aligned_reads" / "{sample}.recal_data.table",
    threads: 2
    shell:
        '''
        gatk BaseRecalibrator \
            -I {input.marked_dups} \
            -R {input.ref_genome} \
            --known-sites {input.known_sites_vcf} \
            -O {output.recal_data_table}
        '''

rule apply_BQSR:
    input:
        recal_data_table = results / "aligned_reads" / "{sample}.recal_data.table",
        marked_dups = results / "aligned_reads" / "{sample}.marked_dups.bam",
        ref_genome = copied_ref_genome,
    output:
        analysis_ready_bam = results / "aligned_reads" / "{sample}.analysis_ready.bam",
        analysys_ready_bam_index = results / "aligned_reads" / "{sample}.analysis_ready.bai"
    threads: 2
    shell:
        '''
        gatk ApplyBQSR \
            -I {input.marked_dups} \
            -R {input.ref_genome} \
            --bqsr-recal-file {input.recal_data_table} \
            -O {output.analysis_ready_bam} \
            --create-output-bam-index true
        '''

rule haplotype_caller:
    input:
        analysis_ready_bam = results / "aligned_reads" / "{sample}.analysis_ready.bam",
        analysys_ready_bam_index = results / "aligned_reads" / "{sample}.analysis_ready.bai",
        ref_genome = copied_ref_genome
    output:
        called_haplotypes = results / "aligned_reads" / "{sample}.g.vcf.gz"
    threads: 2
    shell:
        '''
        gatk HaplotypeCaller \
            -I {input.analysis_ready_bam} \
            -R {input.ref_genome} \
            -ERC GVCF \
            -O {output.called_haplotypes} \
            --alleles {config[targets_vcf]}\
            --output-mode EMIT_ALL_ACTIVE_SITES
        '''
