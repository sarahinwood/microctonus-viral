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
        ##blast results
        expand('output/nr_blastp/{sample}_blastp.outfmt3', sample=all_samples),
        ##exonerate to map genome bro genes identified with blast onto M.hyp genome
        #'output/mh_exonerate/genome_bro/exonerate_bro_scaffolds.txt',
        ##exonerate to map transcriptome bro genes onto M.hyp genome
        'output/mh_exonerate/transcriptome_bro/mh_trinity_broN_exonerate.out',
        ##exonerate to map transcriptome baculoviridae genes onto M.hyp genome
        'output/mh_exonerate/transcriptome_baculoviridae/vul_baculo_exonerate.out',
        ##blast other peptides on the same contigs as genome bro genes
        'output/mh_exonerate/genome_bro_scaffolds/bro_scaffold_peptides_blastp.outfmt3',
        ###viral peptide exonerate
       	'output/mh_exonerate/genome_viral_scaffolds/viral_scaffold_peptides_blastp.outfmt3',
       	##interproscan for all peptides on viral scaffolds
       	'output/mh_exonerate/genome_viral_scaffolds/interproscan/all_viral_scaffold_peptides.fasta.tsv'


##for exonerate
##could use this to decrease output size while keeping vulgar output --showalignment no
##--bestn report best n matches for each query - kind of handy to be abke to search for bro scaffold no.s in full res currently


rule mh_trans_baculoviridae_grep_res:
    input:
        baculoviridae_exonerate = 'output/mh_exonerate/transcriptome_baculoviridae/mh_baculoviridae_exonerate.out'
    output:
        vulgar_baculoviridae_exonerate = 'output/mh_exonerate/transcriptome_baculoviridae/vul_baculo_exonerate.out'
    shell:
        'egrep -i "vulgar:" {input.baculoviridae_exonerate} > {output.vulgar_baculoviridae_exonerate}'

#where do other transriptome baculoviridae genes align to in M.hyp genome?
rule mh_transcriptome_baculoviridae_exonerate:
	input:
		baculoviridae_genes = 'output/mh_exonerate/transcriptome_baculoviridae/baculoviridae_genes.fasta',
		mh_genome = 'data/Mh_assembly.fa'
	output:
		exonerate_res = 'output/mh_exonerate/transcriptome_baculoviridae/mh_baculoviridae_exonerate.out'
	threads:
		20
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
		baculoviridae_gene_ids = 'output/mh_exonerate/transcriptome_baculoviridae/baculoviridae_gene_ids.txt',
		mh_transcriptome = 'data/mh_transcriptome.fasta'
	output:
		baculoviridae_genes = 'output/mh_exonerate/transcriptome_baculoviridae/baculoviridae_genes.fasta'
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

#where to M.hyp transcriptome bro genes align to in M.hyp genome?
rule mh_transcriptome_broN_exonerate:
	input:
		broN_genes = 'output/mh_exonerate/transcriptome_bro/bro-n_domain_trinity_genes.fasta',
		mh_genome = 'data/Mh_assembly.fa'
	output:
		exonerate_res = 'output/mh_exonerate/transcriptome_bro/mh_trinity_broN_exonerate.out'
	threads:
		20
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

##what blastp hits do the other peptides on same scaffolds as bro genes get?
rule blast_bro_scaffold_peptides:
    input:
        bro_scaffold_peptides = 'output/mh_exonerate/genome_bro_scaffolds/bro_scaffold_peptides.fasta'
    output:
        blastp_res = 'output/mh_exonerate/genome_bro_scaffolds/bro_scaffold_peptides_blastp.outfmt3'
    params:
        blast_db = 'bin/db/blastdb/nr/nr'
    threads:
        50
    log:
        'output/logs/nr_blastp_bro_scaffold_peptides.log'
    shell:
        'blastp '
        '-query {input.bro_scaffold_peptides} '
        '-db {params.blast_db} '
        '-num_threads {threads} '
        '-evalue 1e-05 '
        '-outfmt "6 std salltitles" > {output.blastp_res} '
        '2> {log}'

##filter out peptides that map onto bro scaffolds but are NOT bro peptides
rule filter_bro_scaffold_peptides:
    input:
        peptide_list = 'output/mh_exonerate/genome_bro_scaffolds/bro_scaffold_peptide_ids.txt',
        peptide_db = 'data/peptide_dbs/Mhyp.faa'
    output:
        bro_scaffold_peptides = 'output/mh_exonerate/genome_bro_scaffolds/bro_scaffold_peptides.fasta'
    threads:
        50
    singularity:
        bbduk_container
    log:
        'output/logs/filter_bro_scaffold_peptides.log'
    shell:
        'filterbyname.sh '
        'in={input.peptide_db} '
        'include=t '
        'names={input.peptide_list} '
        'ignorejunk=t '
        'out={output.bro_scaffold_peptides} '
        '&> {log}'

#generate list of peptides that map onto bro scaffolds and are NOT the bro peptides
rule list_peptides_on_bro_scaffolds:
    input:
        exonerate_bro_scaffolds_res = 'output/mh_exonerate/genome_bro_scaffolds/bro_scaffold_exonerate.out',
        bro_final_table = 'output/mh_exonerate/genome_bro/bro_genome_exonerate_full_res.csv'
    output:
        bspe_final_table = 'output/mh_exonerate/genome_bro_scaffolds/exonerate_table.csv',
        peptide_list = 'output/mh_exonerate/genome_bro_scaffolds/bro_scaffold_peptide_ids.txt'
    singularity:
        tidyverse_container
    log:
        'output/logs/list_peptides_on_bro_scaffolds.log'
    script:
        'src/list_peptides_on_bro_scaffolds.R'

##exonerate on bro genome scaffolds to see what other peptides are on the same scaffolds as the bro genes?
rule exonerate_bro_scaffolds:
    input:
        bro_scaffolds = 'output/mh_exonerate/genome_bro_scaffolds/bro_scaffolds.fasta',
        mh_peptides = 'data/peptide_dbs/Mhyp.faa'
    output:
        exonerate_res = 'output/mh_exonerate/genome_bro_scaffolds/bro_scaffold_exonerate.out'
    threads:
        20
    log:
        'output/logs/bro_scaffold_exonerate.log'
    shell:
        'bin/exonerate-2.2.0-x86_64/bin/exonerate '
        '--model protein2genome '
        '--score 400 '
        '--showalignment no '
        '{input.mh_peptides} '
        '{input.bro_scaffolds} '
        '> {output.exonerate_res} '
        '2> {log} '

##Pull out genome scaffolds that bro genes map onto
rule filter_bro_scaffolds:
    input:
        bro_scaffold_list = 'output/mh_exonerate/genome_bro/bro_scaffold_ids.txt',
        mh_genome = 'data/Mh_assembly.fa'
    output:
        bro_scaffolds = 'output/mh_exonerate/genome_bro_scaffolds/bro_scaffolds.fasta'
    threads:
        50
    singularity:
        bbduk_container
    log:
        'output/logs/filter_bro_scaffolds.log'
    shell:
        'filterbyname.sh '
        'in={input.mh_genome} '
        'include=t '
        'names={input.bro_scaffold_list} '
        'out={output.bro_scaffolds} '
        '&> {log}'

##generate list of genome scaffolds that bro genes map onto
rule list_bro_scaffold_ids:
    input:
        exonerate_res = 'output/mh_exonerate/genome_bro/mh_genome_bro_exonerate_vulgar.out'
    output:
        bro_scaffold_ids = 'output/mh_exonerate/genome_bro/bro_scaffold_ids.txt',
        bro_final_table = 'output/mh_exonerate/genome_bro/bro_genome_exonerate_full_res.csv'
    singularity:
        tidyverse_container
    log:   
        'output/logs/list_bro_scaffold_ids.log'
    script:
        'src/list_bro_scaffold_ids.R'

##get exonerate res in format easier to manipulate in R
rule mh_genome_bro_grep_res:
    input:
        genome_bro_exonerate = 'output/mh_exonerate/genome_bro/mh_genome_bro_exonerate.out'
    output:
        vulgar_bro_exonerate = 'output/mh_exonerate/genome_bro/mh_genome_bro_exonerate_vulgar.out'
    shell:
        'egrep -i "vulgar:" {input.genome_bro_exonerate} > {output.vulgar_bro_exonerate}'

##where do bro peptides identified with reciprocal blast map to in M.hyp genome?
rule mh_genome_bro_exonerate:
    input:
        bro_peptides = 'output/mh_exonerate/genome_bro/bro_peptides.faa',
        mh_genome = 'data/Mh_assembly.fa'
    output:
        exonerate_res = 'output/mh_exonerate/genome_bro/mh_genome_bro_exonerate.out'
    threads:
        20
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
        ##made list of bro peptides myself as annotation name not all in same format in R
        peptide_hit_ids = 'output/viral_nr_blastp_r/Mh/mh_bro_peptide_ids.txt'
    output:
        bro_peptides = 'output/mh_exonerate/genome_bro/bro_peptides.faa'
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
        'ignorejunk=t '
        'out={output.bro_peptides} '
        '&> {log}'


###########################
##interproscan for all peptides on viral scaffolds
rule interproscan_viral_scaffold_peptides:
	input:
		viral_scaffold_peptides = 'output/mh_exonerate/genome_viral_scaffolds/all_viral_scaffold_peptides.fasta'
	output:
		interpro_tsv = 'output/mh_exonerate/genome_viral_scaffolds/interproscan/all_viral_scaffold_peptides.fasta.tsv'
	params:
		outdir = 'output/mh_exonerate/genome_viral_scaffolds/interproscan'
	threads:
		50
	log:
		'output/logs/interproscan_viral_scaffold_peptides.log'
	shell:
		'bin/interproscan-5.31-70.0/interproscan.sh '
		'--input {input.viral_scaffold_peptides} '
		'--seqtype p '
		'--output-dir {params.outdir} '
		'--cpu {threads} '
		'--goterms '
		'2> {log} '

##filter out all peptides on viral scaffolds
rule filter_all_viral_scaffold_peptides:
    input:
        peptide_list = 'output/mh_exonerate/genome_viral_scaffolds/viral_scaffold_all_peptide_ids.txt',
        peptide_db = 'data/peptide_dbs/Mhyp.faa'
    output:
        viral_scaffold_peptides = 'output/mh_exonerate/genome_viral_scaffolds/all_viral_scaffold_peptides.fasta'
    threads:
        50
    singularity:
        bbduk_container
    log:
        'output/logs/filter_all_viral_scaffold_peptides.log'
    shell:
        'filterbyname.sh '
        'in={input.peptide_db} '
        'include=t '
        'names={input.peptide_list} '
        'ignorejunk=t '
        'out={output.viral_scaffold_peptides} '
        '&> {log}'

##what blastp hits do the other peptides on same scaffolds as bro genes get?
rule blast_viral_scaffold_peptides:
    input:
        viral_scaffold_peptides = 'output/mh_exonerate/genome_viral_scaffolds/viral_scaffold_nv_peptides.fasta'
    output:
        blastp_res = 'output/mh_exonerate/genome_viral_scaffolds/viral_scaffold_peptides_blastp.outfmt3'
    params:
        blast_db = 'bin/db/blastdb/nr/nr'
    threads:
        50
    log:
        'output/logs/nr_blastp_viral_scaffold_peptides.log'
    shell:
        'blastp '
        '-query {input.viral_scaffold_peptides} '
        '-db {params.blast_db} '
        '-num_threads {threads} '
        '-evalue 1e-05 '
        '-outfmt "6 std salltitles" > {output.blastp_res} '
        '2> {log}'

##filter out peptides that map onto viral scaffolds but are NOT the blastp viral peptides
rule filter_viral_scaffold_peptides:
    input:
        peptide_list = 'output/mh_exonerate/genome_viral_scaffolds/viral_scaffold_non_viral_peptide_ids.txt',
        peptide_db = 'data/peptide_dbs/Mhyp.faa'
    output:
        viral_scaffold_peptides = 'output/mh_exonerate/genome_viral_scaffolds/viral_scaffold_nv_peptides.fasta'
    threads:
        50
    singularity:
        bbduk_container
    log:
        'output/logs/filter_viral_scaffold_peptides.log'
    shell:
        'filterbyname.sh '
        'in={input.peptide_db} '
        'include=t '
        'names={input.peptide_list} '
        'ignorejunk=t '
        'out={output.viral_scaffold_peptides} '
        '&> {log}'

#generate list of peptides that map onto viral scaffolds but are NOT the viral peptides from the recip-blastp
rule list_peptides_on_viral_scaffolds:
    input:
        exonerate_viral_scaffolds_res = 'output/mh_exonerate/genome_viral_scaffolds/viral_scaffold_exonerate.out',
        viral_exonerate_table = 'output/mh_exonerate/genome_viral/viral_genome_exonerate_full_res.csv'
    output:
        vspe_final_table = 'output/mh_exonerate/genome_viral_scaffolds/exonerate_table.csv',
        non_viral_peptide_list = 'output/mh_exonerate/genome_viral_scaffolds/viral_scaffold_non_viral_peptide_ids.txt',
        all_peptides_list = 'output/mh_exonerate/genome_viral_scaffolds/viral_scaffold_all_peptide_ids.txt'
    singularity:
        tidyverse_container
    log:
        'output/logs/list_peptides_on_viral_scaffolds.log'
    script:
        'src/list_peptides_on_viral_scaffolds.R'

#exonerate on viral genome scaffolds to see what other peptides are on the same scaffolds as the viral genes?
rule exonerate_viral_scaffolds:
    input:
        viral_scaffolds = 'output/mh_exonerate/genome_viral_scaffolds/viral_scaffolds.fasta',
        mh_peptides = 'data/peptide_dbs/Mhyp.faa'
    output:
        exonerate_res = 'output/mh_exonerate/genome_viral_scaffolds/viral_scaffold_exonerate.out'
    threads:
        20
    log:
        'output/logs/viral_scaffold_exonerate.log'
    shell:
        'bin/exonerate-2.2.0-x86_64/bin/exonerate '
        '--model protein2genome '
        '--score 400 '
        '--showalignment no '
        '{input.mh_peptides} '
        '{input.viral_scaffolds} '
        '> {output.exonerate_res} '
        '2> {log} '

#pull out genome scaffolds that viral genes map onto
rule filter_viral_scaffolds:
    input:
        viral_scaffold_list = 'output/mh_exonerate/genome_viral/viral_scaffold_ids.txt',
        mh_genome = 'data/Mh_assembly.fa'
    output:
        viral_scaffolds = 'output/mh_exonerate/genome_viral_scaffolds/viral_scaffolds.fasta'
    threads:
        50
    singularity:
        bbduk_container
    log:
        'output/logs/filter_viral_scaffolds.log'
    shell:
        'filterbyname.sh '
        'in={input.mh_genome} '
        'include=t '
        'names={input.viral_scaffold_list} '
        'out={output.viral_scaffolds} '
        '&> {log}'

#generate list of genome scaffolds that viral genes map onto
rule list_viral_scaffold_ids:
    input:
        viral_exonerate_res = 'output/mh_exonerate/genome_viral/viral_peptide_exonerate.out'
    output:
        viral_scaffold_ids = 'output/mh_exonerate/genome_viral/viral_scaffold_ids.txt',
        viral_exonerate_table = 'output/mh_exonerate/genome_viral/viral_genome_exonerate_full_res.csv'
    singularity:
        tidyverse_container
    log:   
        'output/logs/list_viral_scaffold_ids.log'
    script:
        'src/list_viral_scaffold_ids.R'

#get exonerate res in format easier to manipulate in R
rule viral_peptide_exonerate_grep_res:
    input:
        viral_exonerate = 'output/mh_exonerate/genome_viral/viral_peptide_exonerate.out'
    output:
        vulgar_viral_exonerate = 'output/mh_exonerate/genome_viral/viral_peptide_exonerate_vulgar.out'
    shell:
        'egrep -i "vulgar:" {input.viral_exonerate} > {output.vulgar_viral_exonerate}'

#find which scaffolds all of the viral peptides map to - LOOK HERE AND SEE IF GENES APPEAR TO HAVE INTRONS!!!!
rule viral_peptide_exonerate:
	input:
		viral_peptides = 'output/mh_exonerate/genome_viral/viral_peptides.faa',
		mh_genome = 'data/Mh_assembly.fa'
	output:
		exonerate_res = 'output/mh_exonerate/genome_viral/viral_peptide_exonerate.out'
	threads:
		20
	log:
		'output/logs/viral_peptide_exonerate.log'
	shell:
		'bin/exonerate-2.2.0-x86_64/bin/exonerate '
        '--model protein2genome '
        '--score 400 '
        '{input.viral_peptides} '
        '{input.mh_genome} '
        '> {output.exonerate_res} '
        '2> {log} '	

#filter out all 40 Mh peptides that had viral hits in recip-blastp
rule filter_viral_peptides:
	input:
		peptide_db = 'data/peptide_dbs/Mhyp.faa',
		viral_peptide_ids = 'output/viral_nr_blastp_r/Mh/Mh_viral_peptide_list.txt'
	output:
		viral_peptides = 'output/mh_exonerate/genome_viral/viral_peptides.faa'
	threads:
		50
	singularity:
		bbduk_container
	log:
		'output/logs/filter_viral_peptides.log'
	shell:
	 	'filterbyname.sh '
        'in={input.peptide_db} '
        'include=t '
        'names={input.viral_peptide_ids} '
        'ignorejunk=t '
        'out={output.viral_peptides} '
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

#technically meant for nt sequence not aa sequence - should probably use something else??
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








