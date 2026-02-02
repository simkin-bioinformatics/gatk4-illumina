rule create_sample_map_file:
    input:
        called_haplotypes = expand(results / "temp" / "called_haplotypes" / "{sample}.g.vcf.gz", sample = samples),
    output:
        sample_map_file = results / "temp" / "cohort.sample_map",
    run:
        with open(output.sample_map_file, 'w') as out:
            for path in input.called_haplotypes:
                sample_name = Path(path).name.removesuffix('.g.vcf.gz')
                out.write(sample_name + '\t' + str(path) + "\n")

rule create_intervals_for_processing:
    input:
        rules.collect_indices.input,
        ref_genome = copied_ref_genome,
    output:
        interval_list = results / "temp" / "intervals.interval_list",
    shell:
        '''
        pixi run gatk PreprocessIntervals \
            -R {input.ref_genome} \
            --bin-length 1000 \
            --padding 0 \
            -O {output.interval_list}
        '''

rule genomics_db_import:
    input:
        rules.collect_indices.input,
        called_haplotypes = expand(results / "temp" / "called_haplotypes" / "{sample}.g.vcf.gz", sample = samples),
        called_haplotypes_index = expand(results / "temp" / "called_haplotypes" / "{sample}.g.vcf.gz.tbi", sample = samples),
        interval_list = results / "temp" / "intervals.interval_list",
        ref_genome = copied_ref_genome,
        sample_map_file = results / "temp" / "cohort.sample_map",
    output:
        database = directory(results / "temp" / "dbimport")
    shell:
        '''
        pixi run gatk GenomicsDBImport \
            --genomicsdb-workspace-path {output.database} \
            -L {input.interval_list} \
            --sample-name-map {input.sample_map_file}
        '''

rule joint_genotyping_with_dbi:
    input:
        rules.collect_indices.input,
        ref_genome = copied_ref_genome,
        database = directory(results / "temp" / "dbimport"),
        targets_vcf = results / "temp" / "targets.vcf",
    output:
        genotyped_cohort = results / 'genotyped_cohort.vcf.gz'
    shell:
        '''
        pixi run gatk GenotypeGVCFs \
            -R {input.ref_genome} \
            -V gendb://{input.database} \
            -O {output.genotyped_cohort} \
            --force-output-intervals {input.targets_vcf}
        '''
