rule create_variant_arguments_file:
    input:
        called_haplotypes = expand(results / "temp" / "called_haplotypes" / "{sample}.g.vcf.gz", sample = samples),
    output:
        variant_arguments_file = results / "temp" / "consolidate_vcf_args.txt",
    run:
        with open(output.variant_arguments_file, 'w') as out:
            for variant in input.called_haplotypes:
                out.write(f"--variant {variant}\n")


rule consolidate_gVCFs:
    input:
        rules.collect_indices.input,
        called_haplotypes = expand(results / "temp" / "called_haplotypes" / "{sample}.g.vcf.gz", sample = samples),
        called_haplotypes_index = expand(results / "temp" / "called_haplotypes" / "{sample}.g.vcf.gz.tbi", sample = samples),
        variant_arguments_file = results / "temp" / "consolidate_vcf_args.txt",
        ref_genome = copied_ref_genome,
    output:
        called_merged_haplotypes = results / "temp" / "called_haplotypes_merged.vcf.gz",
        called_merged_haplotypes_index = results / "temp" / "called_haplotypes_merged.vcf.gz.tbi",
    shell:
        '''
         pixi run gatk CombineGVCFs \
           --R {input.ref_genome} \
           --arguments_file {input.variant_arguments_file} \
           -O {output.called_merged_haplotypes}
        '''

rule joint_genotyping:
    input:
        rules.collect_indices.input,
        ref_genome = copied_ref_genome,
        called_merged_haplotypes = results / "temp" / "called_haplotypes_merged.vcf.gz",
        called_merged_haplotypes_index = results / "temp" / "called_haplotypes_merged.vcf.gz.tbi",
        targets_vcf = results / "temp" / "targets.vcf"
    output:
        genotyped_cohort = results / 'genotyped_cohort.vcf.gz'
    shell:
        '''
        pixi run gatk GenotypeGVCFs \
            -R {input.ref_genome} \
            -V {input.called_merged_haplotypes} \
            -O {output.genotyped_cohort} \
            --force-output-intervals {input.targets_vcf}
        '''
