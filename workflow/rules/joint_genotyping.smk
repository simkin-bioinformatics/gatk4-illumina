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
        sample_map_file = results / "temp" / "cohort.sample_map",
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

checkpoint split_intervals:
    input:
        rules.collect_indices.input,
        ref_genome = copied_ref_genome,
        interval_list = results / "temp" / "intervals.interval_list",
    output:
        interval_folder = directory(results / "temp" / "interval_lists")
    shell:
        '''
        pixi run gatk SplitIntervals \
           -R {input.ref_genome} \
           -L {input.interval_list} \
           --scatter-count 10 \
           -O {output.interval_folder}
        '''

rule genomics_db_import:
    input:
        rules.collect_indices.input,
        called_haplotypes = expand(results / "temp" / "called_haplotypes" / "{sample}.g.vcf.gz", sample = samples),
        called_haplotypes_index = expand(results / "temp" / "called_haplotypes" / "{sample}.g.vcf.gz.tbi", sample = samples),
        interval_folder = results / "temp" / "interval_lists",
        interval_list = results / "temp" / "interval_lists" / "{shard}.interval_list",
        ref_genome = copied_ref_genome,
        sample_map_file = results / "temp" / "cohort.sample_map",
    output:
        database = directory(results / "temp" / "{shard}_dbimport")
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
        database = results / "temp" / "{shard}_dbimport",
        targets_vcf = results / "temp" / "targets.vcf",
    output:
        genotyped_shard = results / "temp" / '{shard}_genotyped.vcf.gz'
    shell:
        '''
        pixi run gatk GenotypeGVCFs \
            -R {input.ref_genome} \
            -V gendb://{input.database} \
            -O {output.genotyped_shard} \
            --force-output-intervals {input.targets_vcf}
        '''

def do_the_genotyping(wildcards):
    checkpoint_output = checkpoints.split_intervals.get(**wildcards).output['interval_folder']
    shards = [file.name.removesuffix('.interval_list') for file in Path(checkpoint_output).glob("*.interval_list")]
    # print(shards)
    return sorted(expand(results / "temp" / "{shard}_genotyped.vcf.gz", shard = shards))

rule gather_vcfs_arguments:
    input:
        shards = do_the_genotyping
    output:
        arguments_file = results / "temp" / "gather_vcfs_arguments.txt"
    run:
        with open(output.arguments_file, 'w') as args:
            for genotyped_shard in input.shards:
                args.write("-I"+"\t"+genotyped_shard+"\n")

rule gather_vcfs:
    input:
        arguments_file = results / "temp" / "gather_vcfs_arguments.txt"
    output:
        genotyped_cohort = results / 'genotyped_cohort.vcf.gz'
    shell:
        '''
        pixi run gatk GatherVcfs \
            --arguments_file {input.arguments_file} \
            -O {output.genotyped_cohort}
        '''
