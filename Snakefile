#!/usr/bin/env python3

import pathlib2
import pandas
import os

#############
# FUNCTIONS #
#############

def resolve_path(x):
    return(str(pathlib2.Path(x).resolve(strict=False)))

####for peptide db files####

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

###########
# GLOBALS #
###########

####for peptide dbs####
peptide_sample_key_file = 'data/peptide_sample_key.csv'
peptide_sample_key = pandas.read_csv(peptide_sample_key_file)
peptide_dir = 'data/peptide_dbs'
all_dbs = sorted(set(peptide_sample_key['Name']))

##############
# CONTAINERS #
##############

tidyverse_container = 'shub://TomHarrop/singularity-containers:r_3.5.0'

#########
# SETUP #
#########

####for peptide db files####
all_faa = find_db_files(peptide_dir)
all_faa_files = sorted(set(peptide_sample_key['Name']))

#########
# RULES #
#########

rule target:
    input:
        ##blast results
        expand('output/nr_blastp/{sample}_blastp.outfmt3', sample=all_faa_files),
        ##exonerate to map transcriptome bro genes onto M.hyp genome
        'output/mh_exonerate/transcriptome_bro/mh_trinity_broN_exonerate.out',
        ##exonerate to map transcriptome baculoviridae genes onto M.hyp genome
        'output/mh_exonerate/transcriptome_baculoviridae/vul_baculo_exonerate.out',
        ###viral peptide exonerate
       	'output/mh_exonerate/genome_viral_scaffolds/viral_scaffold_peptides_blastp.outfmt3',
       	##interproscan for all peptides on viral scaffolds
       	'output/mh_exonerate/genome_viral_scaffolds/interproscan/all_viral_scaffold_peptides.fasta.tsv',
        ##for trial of exonerate with gff3
        'output/mh_exonerate/genome_viral/exonerate_gff_trial/exonerate.output.gff',
        ##blastx of full viral scaffolds
        'output/mh_exonerate/genome_viral_scaffolds/blastx/blastx_viral_scaffolds.outfmt3'

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

#blastx viral scaffolds to see what hits they get and if peptides predicted by blastx differ in location from the peptide dbs
rule blastx_viral_scaffolds:
    input:
        viral_scaffolds = 'output/mh_exonerate/genome_viral_scaffolds/viral_scaffolds.fasta'
    output:
        blastp_res = 'output/mh_exonerate/genome_viral_scaffolds/blastx/blastx_viral_scaffolds.outfmt3'
    params:
        blast_db = 'bin/db/blastdb/nr/nr'
    threads:
        50
    log:
        'output/logs/blastx_viral_scaffolds.log'
    shell:
        'blastx '
        '-query {input.viral_scaffolds} '
        '-db {params.blast_db} '
        '-num_threads {threads} '
        '-evalue 1e-05 '
        '-outfmt "6 std salltitles" > {output.blastp_res} '
        '2> {log}'

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

##reformat exonerate gff output
rule process_exonerate_gff3:
    input:
        exonerate_gff_output = 'output/mh_exonerate/genome_viral/exonerate_gff_trial/viral_peptide_exonerate.out'
    output:
        gff3 = 'output/mh_exonerate/genome_viral/exonerate_gff_trial/exonerate.output.gff'
    threads:
        20
    shell:
        'src/process_exonerate_gff3.pl {input.exonerate_gff_output} > {output.gff3}'

##same but with gff to see what that output looks like
rule viral_peptide_exonerate_gff:
    input:
        viral_peptides = 'output/mh_exonerate/genome_viral/viral_peptides.faa',
        mh_genome = 'data/Mh_assembly.fa'
    output:
        exonerate_res = 'output/mh_exonerate/genome_viral/exonerate_gff_trial/viral_peptide_exonerate.out'
    threads:
        20
    log:
        'output/logs/viral_peptide_exonerate_gff3.log'
    shell:
        'bin/exonerate-2.2.0-x86_64/bin/exonerate '
        '--model protein2genome '
        '--score 400 '
        '{input.viral_peptides} '
        '{input.mh_genome} '
        '> {output.exonerate_res} '
        '--showtargetgff yes '
        '--showalignment no '
        '--ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n" '
        '--bestn 1 '
        '2> {log} ' 

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








