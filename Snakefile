#!/usr/bin/env python3

import peppy

#############
# FUNCTIONS #
#############

def get_peptide_dbs(wildcards):
    input_keys = ['peptide_db']
    my_pep = pep.get_sample(wildcards.sample).to_dict()
    return {k: my_pep[k] for k in input_keys}

###########
# GLOBALS #
###########

##this parses the config & sample key files into an object named pep
pepfile: 'data/config.yaml'
##can now use this to generate list of all samples
all_samples = pep.sample_table['sample_name']

##############
# CONTAINERS #
##############

tom_tidyverse_container = 'shub://TomHarrop/singularity-containers:r_3.5.0'
tidyverse_container = 'docker://rocker/tidyverse'
bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'

#########
# RULES #
#########

rule target:
    input:
        ##recip blast -> all peptides on viral scaffolds blast results
        expand('output/nr_analysis_{species}/{species}_blastp_peptides_viral_scaffolds.outfmt6', species=all_samples),
        ##interproscan results
        expand('output/nr_analysis_{species}/{species}_interproscan.tsv', species=all_samples)
        ##exonerate to map transcriptome bro genes onto M.hyp genome
        #'output/mh_exonerate/transcriptome_bro/mh_trinity_broN_exonerate.out',
        ##exonerate to map transcriptome baculoviridae genes onto M.hyp genome
        #'output/mh_exonerate/transcriptome_baculoviridae/vul_baculo_exonerate.out',
        ###viral peptide exonerate
        #'output/mh_exonerate/genome_viral_scaffolds/viral_scaffold_peptides_blastp.outfmt6',
        ##interproscan for all peptides on viral scaffolds
        #'output/mh_exonerate/genome_viral_scaffolds/interproscan/all_viral_scaffold_peptides.fasta.tsv',
        ##for trial of exonerate with gff3 - errors in output - must fix manually!
        #'output/mh_exonerate/genome_viral/exonerate_gff_trial/exonerate.output.gff',
        ##blastx of full viral scaffolds
        #'output/mh_exonerate/genome_viral_scaffolds/blastx/blastx_viral_scaffolds.outfmt6'

rule mh_trans_baculoviridae_grep_res:
    input:
        baculoviridae_exonerate = 'output/mh_exonerate/transcriptome_baculoviridae/mh_baculoviridae_exonerate.out'
    output:
        vulgar_baculoviridae_exonerate = 'output/mh_exonerate/transcriptome_baculoviridae/vul_baculo_exonerate.out'
    shell:
        'egrep -i "vulgar:" {input.baculoviridae_exonerate} > {output.vulgar_baculoviridae_exonerate}'

#where do other transcriptome baculoviridae genes align to in M.hyp genome?
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

##should also look at genes and whether they have introns - genes predicted with eukaryotic translation code so need to re-predict with bacterial

#############################################
## what else on scaffolds with viral hits? ##
#############################################

rule interproscan_peptides_viral_scaffolds:
    input:
        peptides_viral_scaffolds = 'output/nr_analysis_{species}/{species}_peptides_viral_scaffolds.faa'
    output:
        interpro_tsv = 'output/nr_analysis_{species}/{species}_interproscan.tsv'
    threads:
        20
    log:
        'output/logs/{species}/interproscan_peptides_viral_scaffolds.log'
    shell:
        'bin/interproscan-5.51-85.0/interproscan.sh '
        '--input {input.peptides_viral_scaffolds} '
        '--formats TSV '
        '--outfile {output.interpro_tsv} '
        '--goterms '
        '--cpu {threads} '
        '2> {log}'

rule blastp_peptides_viral_scaffolds:
    input:
        peptides_viral_scaffolds = 'output/nr_analysis_{species}/{species}_peptides_viral_scaffolds.faa'
    output:
        blastp_res = 'output/nr_analysis_{species}/{species}_blastp_peptides_viral_scaffolds.outfmt6'
    params:
        blast_db = 'bin/db/blastdb/nr/nr'
    threads:
        20
    log:
        'output/logs/{species}/blastp_peptides_viral_scaffolds.log'
    shell:
        'blastp '
        '-query {input.peptides_viral_scaffolds} '
        '-db {params.blast_db} '
        '-num_threads {threads} '
        '-evalue 1e-05 '
        '-max_target_seqs 1 '
        '-outfmt "6 std staxids salltitles" > {output.blastp_res} '
        '2> {log}'

rule filter_peptides_viral_scaffolds:
    input:
        peptide_db = 'data/final_microctonus_assemblies_annotations/{species}.proteins.fa',
        peptides_viral_scaffolds = 'output/nr_analysis_{species}/{species}_peptides_viral_scaffolds.txt'
    output:
        peptides_viral_scaffolds = 'output/nr_analysis_{species}/{species}_peptides_viral_scaffolds.faa'
    threads:
        20
    singularity:
        bbduk_container
    log:
        'output/logs/{species}/filter_viral_peptides.log'
    shell:
        'filterbyname.sh '
        'in={input.peptide_db} '
        'include=t '
        'names={input.peptides_viral_scaffolds} '
        'ignorejunk=t '
        'out={output.peptides_viral_scaffolds} '
        '&> {log}'  

rule link_viral_peptides_to_scaffolds:
    input:
        viral_peptide_list = 'output/nr_analysis_{species}/{species}_viral_peptides.txt',
        gff = 'data/final_microctonus_assemblies_annotations/{species}.gff3'
    output:
        peptides_viral_scaffolds = 'output/nr_analysis_{species}/{species}_peptides_viral_scaffolds.txt'
    singularity:
        tidyverse_container
    log:
        'output/logs/{species}/link_viral_peptides_to_scaffolds.log'
    script:
        'src/list_peptides_on_viral_scaffolds.R'

#############################
## reciprocal viral blastp ##
#############################

rule blastp_nr_viral_res:
    input:
        nr_blastp_res = 'output/nr_blastp/{species}_blastp.outfmt6'
    output:
        blastp_viral_res = 'output/nr_analysis_{species}/{species}_nr_blast_viral.csv',
        viral_peptide_list = 'output/nr_analysis_{species}/{species}_viral_peptides.txt'
    singularity:
        tidyverse_container
    log:
        'output/logs/{species}/blastp_nr_viral_res.log'
    script:
        'src/nr_blastp_viral_analysis.R'

##should test with '-max_target_seqs 1 ' and compare res
rule blastp_nr:
    input:
        pot_viral_peptides = 'output/viral_blastp/{species}_potential_viral_peptides.faa'
    output:
        blastp_res = 'output/nr_blastp/{species}_blastp.outfmt6'
    params:
        blast_db = 'bin/db/blastdb/nr/nr'
    threads:
        50
    log:
        'output/logs/{species}/nr_blastp.log'
    shell:
        'blastp '
        '-query {input.pot_viral_peptides} '
        '-db {params.blast_db} '
        '-num_threads {threads} '
        '-evalue 1e-05 '
        '-outfmt "6 std staxids salltitles" > {output.blastp_res} '
        '2> {log}'

rule filter_pot_viral_peptides:
    input:
        peptide_db = 'data/final_microctonus_assemblies_annotations/{species}.proteins.fa',
        peptide_hit_ids = 'output/viral_blastp/{species}_peptide_hit_ids.txt'
    output:
        pot_viral_peptides = 'output/viral_blastp/{species}_potential_viral_peptides.faa'
    threads:
        50
    singularity:
        bbduk_container
    log:
        'output/logs/{species}/filter_pot_viral_peptides.log'
    shell:
        'filterbyname.sh '
        'in={input.peptide_db} '
        'include=t '
        'names={input.peptide_hit_ids} '
        'substring=name '
        'ignorejunk=t '
        'out={output.pot_viral_peptides} '
        '&> {log}'

rule filter_pot_viral_peptide_ids:
    input:
        blastp_res = 'output/viral_blastp/{species}_blastp.outfmt6'
    output:
        peptide_hit_ids = 'output/viral_blastp/{species}_peptide_hit_ids.txt'
    singularity:
        tidyverse_container
    log:
        'output/logs/{species}/peptide_hit_ids.log'
    script:
        'src/peptide_hit_id_lists.R'

##should test with '-max_target_seqs 1 ' and compare res
rule blastp_viral:
    input:
        query = 'data/final_microctonus_assemblies_annotations/{species}.proteins.fa',
        gi_list = 'data/gi_lists/virus.gi.txt'
    output:
        blastp_res = 'output/viral_blastp/{species}_blastp.outfmt6'
    params:
        blast_db = 'bin/db/blastdb/nr/nr'
    threads:
        50
    log:
        'output/logs/{species}/viral_blastp.log'
    shell:
        'blastp '
        '-query {input.query} '
        '-db {params.blast_db} '
        '-gilist {input.gi_list} '
        '-num_threads {threads} '
        '-evalue 1e-05 '
        '-outfmt "6 std staxids salltitles" > {output.blastp_res} '
        '2> {log}'
