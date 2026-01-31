copied_ref_genome = str(Path(config['indexed_input_directory']) / Path(config['genome_fasta']).name)
copied_known_sites_vcf = Path(config['indexed_input_directory']) / Path(config['known_sites_vcf']).name

import pandas as pd
import subprocess

rule copy_genome_and_vcf:
    input:
        ref_genome = Path(config['genome_fasta']).resolve(),
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

rule create_targets_vcf:
    input:
        gatk_ref_dict = copied_ref_genome.replace(".fasta", ".dict"),
        targets_tsv = config['targets_tsv']
    output:
        targets_vcf_no_header = results / "temp" / "targets_no_header.vcf",
        targets_vcf_with_header = results / "temp" / "targets.vcf"
    run:
        df = pd.read_table(input.targets_tsv)[['CHROM', 'POS', 'ID', 'REF', 'ALT']]
        df[['QUAL', 'FILTER', 'INFO']] = ['.', 'PASS', '.']
        df = df.rename(columns = {'CHROM': "#CHROM"})

        df.to_csv(output.targets_vcf_no_header, sep='\t', index=False)
        with open(output.targets_vcf_no_header, 'r') as f:
            contents = f.readlines()
        contents.insert(0, '##fileformat=VCFv4.2' + '\n')
        with open(output.targets_vcf_no_header, 'w') as f:
            f.writelines(contents)

        subprocess.run([
            "gatk", "UpdateVCFSequenceDictionary", 
            "-V", output.targets_vcf_no_header, 
            "--source-dictionary", input.gatk_ref_dict,
            "-O", output.targets_vcf_with_header
        ])