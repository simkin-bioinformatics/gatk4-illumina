rule bwa_mem:
    input:
        rules.collect_indices.input,
        ref_genome = copied_ref_genome,
        read_1 = Path(config['sample_reads_folder']) / str("{sample}"+config['R1_suffix']),
        read_2 = Path(config['sample_reads_folder']) / str("{sample}" + config['R2_suffix']),
    output:
        bam = results / "temp" / "aligned_reads" / "{sample}.bam"
    shell:
        "(bwa mem "
        "-R '@RG\\tID:group1\\tSM:{wildcards.sample}\\tPL:ILLUMINA\\tLB:lib1' "
        "{input.ref_genome} {input.read_1} {input.read_2} | "
        "samtools sort -o {output.bam} -)"

rule mark_duplicates:
    input:
        bam = results / "temp" / "aligned_reads" / "{sample}.bam"
    output:
        marked_dups = results / "temp" / "marked_duplicates" / "{sample}.marked_dups.bam",
        dup_metrics = results / "temp" / "marked_duplicates" / "{sample}.dup_metrics.txt",
        marked_dups_index = results / "temp" / "marked_duplicates" / "{sample}.marked_dups.bai",
    threads: 2
    shell:
        '''
        pixi run gatk MarkDuplicates \
            -I {input.bam} \
            -O {output.marked_dups} \
            -M {output.dup_metrics} \
            --CREATE_INDEX true
       '''
       # rules.bwa_index.output,
       # rules.samtools_index.output,
       # rules.index_feature_file.output,
       # rules.gatk_Sequence_Dictionary.output,
       # rules.create_targets_vcf.output,
rule BQSR:
    input:
        rules.collect_indices.input,
        ref_genome = copied_ref_genome,
        known_sites_vcf = copied_known_sites_vcf,
        marked_dups = results / "temp" / "marked_duplicates" / "{sample}.marked_dups.bam",
        marked_dups_index = results / "temp" / "marked_duplicates" / "{sample}.marked_dups.bai"
    output:
        recal_data_table = results / "temp" / "BQSR" / "{sample}.recal_data.table",
    threads: 2
    shell:
        '''
        pixi run gatk BaseRecalibrator \
            -I {input.marked_dups} \
            -R {input.ref_genome} \
            --known-sites {input.known_sites_vcf} \
            -O {output.recal_data_table}
        '''

rule apply_BQSR:
    input:
        rules.collect_indices.input,
        recal_data_table = results / "temp" / "BQSR" / "{sample}.recal_data.table",
        marked_dups = results / "temp" / "marked_duplicates" / "{sample}.marked_dups.bam",
        ref_genome = copied_ref_genome,
    output:
        analysis_ready_bam = results / "temp" / "BQSR" / "{sample}.analysis_ready.bam",
        analysys_ready_bam_index = results / "temp" / "BQSR" / "{sample}.analysis_ready.bai",
    threads: 2
    shell:
        '''
        pixi run gatk ApplyBQSR \
            -I {input.marked_dups} \
            -R {input.ref_genome} \
            --bqsr-recal-file {input.recal_data_table} \
            -O {output.analysis_ready_bam} \
            --create-output-bam-index true
        '''

rule haplotype_caller:
    input:
        rules.collect_indices.input,
        analysis_ready_bam = results / "temp" / "BQSR" / "{sample}.analysis_ready.bam",
        analysys_ready_bam_index = results / "temp" / "BQSR" / "{sample}.analysis_ready.bai",
        ref_genome = copied_ref_genome,
        targets_vcf = results / "temp" / "targets.vcf"
    output:
        called_haplotypes = results / "temp" / "called_haplotypes" / "{sample}.g.vcf.gz",
        called_haplotypes_index = results / "temp" / "called_haplotypes" / "{sample}.g.vcf.gz.tbi",
    threads: 2
    shell:
        '''
        pixi run gatk HaplotypeCaller \
            -I {input.analysis_ready_bam} \
            -R {input.ref_genome} \
            -ERC GVCF \
            -O {output.called_haplotypes} \
            --alleles {input.targets_vcf}\
            --output-mode EMIT_ALL_ACTIVE_SITES
        '''
