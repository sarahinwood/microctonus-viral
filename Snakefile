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
        expand('output/nr_blastp/{sample}_blastp.outfmt3', sample=all_samples)

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
        gi_list = 'data/bracovirus_ichnovirus_sequence.gi.txt'
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








