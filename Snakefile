#!/usr/bin/env python3

import pathlib2
import pandas
import os

#############
# FUNCTIONS #
#############

def resolve_path(x):
    return(str(pathlib2.Path(x).resolve(strict=False)))

def find_db_files(peptide_dir):
#Make list of files
    path_generator = os.walk(peptide_dir, followlinks = True)
    my_files = list((dirpath, filenames)
        for (dirpath, dirname, filenames)
        in path_generator)
#Make new dictionary & populate with files
    my_peptide_files = {}
    for dirpath, filenames in my_files:
        for filename in filenames:
            if filename.endswith('.faa'):
                my_peptide_files = str(pathlib2.Path(dirpath))
    return(my_peptide_files)

def sample_name_to_fastq(wildcards):
    sample_row = sample_key[sample_key['Name'] == wildcards.sample]
    sample_all_faa = [x for x in all_faa[sample_row]]

###########
# GLOBALS #
###########

sample_key_file = 'data/sample_key.csv'
sample_key = pandas.read_csv(sample_key_file)

peptide_dir = 'data/peptide_dbs'

all_dbs = sorted(set(sample_key['Name']))

#containers
tidyverse_container = 'shub://TomHarrop/singularity-containers:r_3.5.0'
bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'

#########
# SETUP #
#########

all_faa = find_db_files(peptide_dir)

all_samples = sorted(set(sample_key['Name']))

#########
# RULES #
#########

rule target:
    input:
        expand('output/nr_blastp/{sample}_blastp.outfmt3', sample=all_samples),
        'output/mh_exonerate/exonerate_bro_scaffolds.txt',
        'output/mh_exonerate/mh_trinity_broN_exonerate.out',
        'output/mh_exonerate/mh_baculoviridae_exonerate.out'

rule mh_transcriptome_baculoviridae_exonerate:
	input:
		baculoviridae_genes = 'output/mh_exonerate/baculoviridae_genes.fasta',
		mh_genome = 'data/Mh_assembly.fa'
	output:
		exonerate_res = 'output/mh_exonerate/mh_baculoviridae_exonerate.out'
	threads:
		50
	log:
		'output/logs/mh_transcriptome_baculoviridae_exonerate.log'
	shell:
	    'bin/exonerate-2.2.0-x86_64/bin/exonerate '
        '--model est2genome '
        '--score 400 '
        '{input.baculoviridae_genes} '
        '{input.mh_genome} '
        '> {output.exonerate_res} '
        '2> {log} '

rule filter_baculoviridae_genes:
	input:
		baculoviridae_gene_ids = 'output/mh_exonerate/baculoviridae_gene_ids.txt',
		mh_transcriptome = 'data/mh_transcriptome.fasta'
	output:
		baculoviridae_genes = 'output/mh_exonerate/baculoviridae_genes.fasta'
	threads:
		50
	singularity:
		bbduk_container
	log:
		'output/logs/filter_baculoviridae_genes.log'
	shell:
	    'filterbyname.sh '
        'in={input.mh_transcriptome} '
        'include=t '
        'names={input.baculoviridae_gene_ids} '
        'out={output.baculoviridae_genes} '
        '&> {log}'

rule mh_transcriptome_broN_exonerate:
	input:
		broN_genes = 'output/mh_exonerate/bro-n_domain_trinity_genes.fasta',
		mh_genome = 'data/Mh_assembly.fa'
	output:
		exonerate_res = 'output/mh_exonerate/mh_trinity_broN_exonerate.out'
	threads:
		50
	log:
		'output/logs/mh_transcriptome_broN_exonerate.log'
	shell:
	    'bin/exonerate-2.2.0-x86_64/bin/exonerate '
        '--model est2genome '
        '--score 400 '
        '{input.broN_genes} '
        '{input.mh_genome} '
        '> {output.exonerate_res} '
        '2> {log} '

rule grep_scaffolds:
    input:
        exonerate_res = 'output/mh_exonerate/mh_exonerate.out'
    output:
        bro_scaffolds = 'output/mh_exonerate/exonerate_bro_scaffolds.txt'
    threads:
        20
    shell:
        'egrep -i "Target: " {input.exonerate_res} > {output.bro_scaffolds} '

rule mh_genome_bro_exonerate:
    input:
        bro_peptides = 'output/mh_exonerate/bro_peptides.faa',
        mh_genome = 'data/Mh_assembly.fa'
    output:
        exonerate_res = 'output/mh_exonerate/mh_exonerate.out'
    threads:
        50
    log:
        'output/logs/mh_exonerate.log'
    shell:
        'bin/exonerate-2.2.0-x86_64/bin/exonerate '
        '--model protein2genome '
        '--score 400 '
        '{input.bro_peptides} '
        '{input.mh_genome} '
        '> {output.exonerate_res} '
        '2> {log} '

rule filter_bro_peptides:
    input:
        peptide_db = 'data/peptide_dbs/Mhyp.faa',
        peptide_hit_ids = 'output/viral_nr_blastp_r/Mh/mh_bro_peptide_ids.txt'
    output:
        bro_peptides = 'output/mh_exonerate/bro_peptides.faa'
    threads:
        50
    singularity:
        bbduk_container
    log:
        'output/logs/filter_bro_peptides.log'
    shell:
        'filterbyname.sh '
        'in={input.peptide_db} '
        'include=t '
        'names={input.peptide_hit_ids} '
        'substring=name '
        'ignorejunk=t '
        'out={output.pot_viral_peptides} '
        '&> {log}'

rule blastp_nr:
    input:
        pot_viral_peptides = 'output/viral_blastp/{sample}_potential_viral_peptides.faa'
    output:
        blastp_res = 'output/nr_blastp/{sample}_blastp.outfmt3'
    params:
        blast_db = 'bin/db/blastdb/nr/nr'
    threads:
        50
    log:
        'output/logs/nr_blastp_{sample}.log'
    shell:
        'blastp '
        '-query {input.pot_viral_peptides} '
        '-db {params.blast_db} '
        '-num_threads {threads} '
        '-evalue 1e-05 '
        '-outfmt "6 std salltitles" > {output.blastp_res} '
        '2> {log}'

##technically meant for nt sequence not aa sequence - should probably use something else??
rule filter_peptides:
    input:
        peptide_db = 'data/peptide_dbs/{sample}.faa',
        peptide_hit_ids = 'output/viral_blastp/{sample}_peptide_hit_ids.txt'
    output:
        pot_viral_peptides = 'output/viral_blastp/{sample}_potential_viral_peptides.faa'
    threads:
        50
    singularity:
        bbduk_container
    log:
        'output/logs/{sample}_filter_pot_viral_peptides.log'
    shell:
        'filterbyname.sh '
        'in={input.peptide_db} '
        'include=t '
        'names={input.peptide_hit_ids} '
        'substring=name '
        'ignorejunk=t '
        'out={output.pot_viral_peptides} '
        '&> {log}'

rule filter_peptide_ids:
    input:
        blastp_res = 'output/viral_blastp/{sample}_blastp.outfmt3'
    output:
        peptide_hit_ids = 'output/viral_blastp/{sample}_peptide_hit_ids.txt'
    singularity:
        tidyverse_container
    log:
        'output/logs/peptide_hit_ids_{sample}.log'
    script:
        'src/peptide_hit_id_lists.R'

rule blastp_viral:
    input:
        query = 'data/peptide_dbs/{sample}.faa',
        gi_list = 'data/gi_lists/virus.gi.txt'
    output:
        blastp_res = 'output/viral_blastp/{sample}_blastp.outfmt3'
    params:
        blast_db = 'bin/db/blastdb/nr/nr'
    threads:
        50
    log:
        'output/logs/viral_blastp_{sample}.log'
    shell:
        'blastp '
        '-query {input.query} '
        '-db {params.blast_db} '
        '-gilist {input.gi_list} '
        '-num_threads {threads} '
        '-evalue 1e-05 '
        '-outfmt "6 std salltitles" > {output.blastp_res} '
        '2> {log}'








