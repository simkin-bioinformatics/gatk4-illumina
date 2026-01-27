copied_ref_genome = str(Path('copied_inputs') / Path(config['ref_genome']).name)
copied_known_sites_vcf = Path('copied_inputs') / Path(config['known_sites_vcf']).name

rule copy_genome_and_vcf:
    input:
        ref_genome = Path(config['ref_genome']).resolve(),
        known_sites_vcf = Path(config['known_sites_vcf']).resolve()
    output:
        ref_genome = copied_ref_genome,
        known_sites_vcf = copied_known_sites_vcf
    shell:
        'rsync -a {input.ref_genome} {output.ref_genome} \n'
        'rsync -a {input.known_sites_vcf} {output.known_sites_vcf}'

rule bwa_index:
    input:
        ref_genome = copied_ref_genome,
    output:
        bwa_indices = multiext(copied_ref_genome, ".amb", ".ann", ".bwt", ".pac", ".sa"),
    shell:
        "bwa index {input.ref_genome}"

rule samtools_index:
    input:
        ref_genome = copied_ref_genome
    output:
        sam_index = f"{copied_ref_genome}.fai"
    shell:
        "samtools faidx {input.ref_genome}"

rule index_feature_file:
    input:
        known_sites_vcf = copied_known_sites_vcf
    output:
        known_sites_index = f"{copied_known_sites_vcf}.tbi"
    shell:
        "gatk IndexFeatureFile -I {input.known_sites_vcf}"

rule gatk_Sequence_Dictionary:
    input:
        ref_genome = copied_ref_genome
    output:
        gatk_ref_dict = copied_ref_genome.replace(".fasta", ".dict"),
    shell:
        "gatk CreateSequenceDictionary -R {input.ref_genome}"
