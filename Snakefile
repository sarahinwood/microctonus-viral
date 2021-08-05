#!/usr/bin/env python3

import peppy
import pathlib2

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
all_species = pep.sample_table['sample_name']

##note if using genome files - use {species}.scaffolds.fa files (not {species}.assembly.fa)
##they have PGA scaffodls renamed to match GFF3

##############
# CONTAINERS #
##############

tom_tidyverse_container = 'shub://TomHarrop/singularity-containers:r_3.5.0'
tidyverse_container = 'docker://rocker/tidyverse'
bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
orffinder_container = 'docker://unlhcc/orffinder:0.4.3'
blast_container = 'docker://ncbi/blast:2.12.0'

##need to run actually using container
prodigal_container = 'docker://biocontainers/prodigal:v1-2.6.3-4-deb_cv1'

#########
# RULES #
#########

rule target:
    input:
        ##recip blast -> all peptides on viral scaffolds blast results
        expand('output/viral_contigs_blastp/{species}/{species}_blastp_peptides_viral_contigs.outfmt6', species=all_species),
        expand('output/final_blast_tables/csv/{species}.csv', species=['Mh', 'MO', 'FR']),
        ##interproscan results - IR results all on Hi-C contigs so not running
        expand('output/interpro/{species}_interpro_table.csv', species=['Mh', 'MO', 'FR']),
        ##prodigal predictions
        expand('output/prodigal/{species}/blastp_gff.csv', species=['Mh', 'MO', 'FR']),
        ##DNA virus contig stats
        expand('output/bbstats/{species}_bb_stats.out', species=['Mh', 'MO', 'FR'])

##############
## interpro ##
##############

rule interpro_analysis:
    input:
        interpro_res = 'output/interpro/{species}_interproscan.tsv'
    output:
        interpro_table = 'output/interpro/{species}_interpro_table.csv'
    threads:
        20
    script:
        'src/prodigal_interpro.R'

##prodigal adds * as stop codon - need to remove first
rule interproscan_peptides_viral_contigs:
    input:
        peptides_viral_contigs = 'output/interpro/{species}_protein_translations.faa'
    output:
        interpro_tsv = 'output/interpro/{species}_interproscan.tsv'
    threads:
        20
    log:
        'output/logs/{species}/interproscan_peptides_viral_contigs.log'
    shell:
        'bin/interproscan-5.51-85.0/interproscan.sh '
        '--input {input.peptides_viral_contigs} '
        '--formats TSV '
        '--outfile {output.interpro_tsv} '
        '--goterms '
        '--cpu {threads} '
        '2> {log}'

###################################
## re-predict viral contig genes ##
###################################

rule prodigal_blastp_analysis:
    input:
        prodigal_blastp_res = 'output/prodigal/{species}/prodigal_blastp.outfmt6',
        prodigal_gff_file = 'output/prodigal/{species}/gene_predictions.gff'
    output:
        prodigal_blastp_res_table = 'output/prodigal/{species}/prodigal_blastp_res.csv',
        blastp_gff = 'output/prodigal/{species}/blastp_gff.csv'
    threads:
        20
    log:
        'output/logs/{species}/prodigal_blastp_analysis.log'
    script:
        'src/prodigal_blastp_analysis.R'

##Do I run this step again with lower evalue threshold for all 3 species?
rule blast_prodigal_predictions:
    input:
        prodigal = 'output/prodigal/{species}/protein_translations.faa'
    output:
        blastp_res = 'output/prodigal/{species}/prodigal_blastp.outfmt6'
    params:
        blast_db = 'bin/db/blastdb/nr/nr'
    threads:
        20
    singularity:
        blast_container
    log:
        'output/logs/{species}/prodigal_blastp.log'
    shell:
        'blastp '
        '-query {input.prodigal} '
        '-db {params.blast_db} '
        '-num_threads {threads} '
        '-evalue 1e-05 '
        '-outfmt "6 std staxids salltitles" > {output.blastp_res} '
        '2> {log}'

##re-predict viral genes using bacterial translation code
rule prodigal:
    input:
        viral_contigs = 'output/viral_contigs_blastp/{species}/{species}_DNA_virus_contigs.faa'
    output:
        protein_translations = 'output/prodigal/{species}/protein_translations.faa',
        nucleotide_seq = 'output/prodigal/{species}/nucleotide_seq.fasta',
        gene_predictions = 'output/prodigal/{species}/gene_predictions.gff'
    log:
        'output/logs/{species}/prodigal.log'
    threads:
        1
    singularity:
        prodigal_container
    shell:
        'prodigal '
        '-i {input.viral_contigs} '
        '-a {output.protein_translations} '
        '-d {output.nucleotide_seq} '
        '-f gff '
        '-p meta '
        '-o {output.gene_predictions} '
        '2> {log} ' 

########################
## viral contig stats ##
########################

rule bb_stats:
    input:
        DNA_virus_contigs = 'output/viral_contigs_blastp/{species}_DNA_virus_contigs.faa'
    output:
        stats = 'output/bbstats/{species}_bb_stats.out'
    log:
        'output/logs/bbstats/{species}_viral_contigs.log'
    singularity:
        bbduk_container
    shell:
        'stats.sh '
        'in={input.DNA_virus_contigs} '
        'out={output.stats} '
        '2> {log}'

rule filter_DNA_viral_contigs:
    input:
        peptide_db = 'data/final_genomes/{species}.fa',
        DNA_virus_contigs = 'output/viral_contigs_blastp/{species}/{species}_DNA_virus_contigs.txt'
    output:
        DNA_virus_contigs = 'output/viral_contigs_blastp/{species}/{species}_DNA_virus_contigs.faa'
    threads:
        20
    singularity:
        bbduk_container
    log:
        'output/logs/{species}/filter_viral_contigs.log'
    shell:
        'filterbyname.sh '
        'in={input.peptide_db} '
        'include=t '
        'names={input.DNA_virus_contigs} '
        'ignorejunk=t '
        'out={output.DNA_virus_contigs} '
        '&> {log}'  

#############################################
## what else on scaffolds with viral hits? ##
#############################################

rule blastp_viral_contigs_analysis:
    input:
        viral_contig_blast_res = 'output/viral_contigs_blastp/{species}/{species}_blastp_peptides_viral_contigs.outfmt6',
        nr_viral_blast_res = 'output/nr_analysis/{species}/{species}_nr_blast_viral.csv',
        gff_file = 'data/final_microctonus_assemblies_annotations/{species}.gff3'
    output:
        full_blast_table = 'output/final_blast_tables/csv/{species}.csv'
    log:
        'output/logs/{species}/viral_contig_blast_analysis.log'
    threads:
        20
    singularity:
        tidyverse_container
    script:
        'src/blastp_viral_contigs_analysis.R'

rule blastp_peptides_viral_contigs:
    input:
        peptides_viral_contigs = 'output/viral_contigs_blastp/{species}/{species}_{peptides_viral_contigs.faa}'
    output:
        blastp_res = 'output/viral_contigs_blastp/{species}/{species}_blastp_peptides_viral_contigs.outfmt6'
    params:
        blast_db = 'bin/db/blastdb/nr/nr'
    threads:
        20
    singularity:
        blast_container
    log:
        'output/logs/{species}/blastp_peptides_viral_contigs.log'
    shell:
        'blastp '
        '-query {input.peptides_viral_contigs} '
        '-db {params.blast_db} '
        '-num_threads {threads} '
        '-evalue 1e-05 '
        '-max_target_seqs 1 '
        '-outfmt "6 std staxids salltitles" > {output.blastp_res} '
        '2> {log}'

rule filter_peptides_viral_contigs:
    input:
        peptide_db = 'data/final_microctonus_assemblies_annotations/{species}.proteins.fa',
        peptides_viral_contigs = 'output/viral_contigs_blastp/{species}/{species}_peptides_viral_contigs.txt'
    output:
        peptides_viral_contigs = 'output/viral_contigs_blastp/{species}/{species}_peptides_viral_contigs.faa'
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
        'names={input.peptides_viral_contigs} '
        'ignorejunk=t '
        'out={output.peptides_viral_contigs} '
        '&> {log}'  

##removes scaffolds that are main-genome from hic to avoid blasting every gene on them
##and peptides already ID'd in viral BlastP - no need to Blast again
rule list_peptides_on_viral_contigs:
    input:
        viral_peptide_list = 'output/nr_analysis/{species}/{species}_viral_peptides.txt',
        gff_file = 'data/final_microctonus_assemblies_annotations/{species}.gff3',
        hic_scaffold_list = 'data/hi-c_genomes/{species}_hic_scaffold_ids.txt',
        virus_info_table = 'output/nr_analysis/{species}/{species}_nr_blast_viral_plot.csv'
    output:
        peptides_viral_contigs = 'output/viral_contigs_blastp/{species}/{species}_peptides_viral_contigs.txt',
        contig_counts = 'output/viral_contigs_blastp/{species}/{species}_viral_contig_counts.csv',
        contig_to_viral_genome = 'output/viral_contigs_blastp/{species}/{species}_contig_viral_genome.csv',
        DNA_virus_contigs = 'output/viral_contigs_blastp/{species}/{species}_DNA_virus_contigs.txt'
    singularity:
        tidyverse_container
    log:
        'output/logs/{species}/link_viral_peptides_to_scaffolds.log'
    script:
        'src/list_peptides_on_viral_contigs.R'

#############################
## reciprocal viral blastp ##
#############################

rule nr_blastp_viral_analysis:
    input:
        nr_blastp_res = 'output/nr_blastp/{species}_blastp.outfmt6',
        taxid_list = 'output/taxids/species_virus_taxids.txt'
    output:
        blastp_viral_res = 'output/nr_analysis/{species}/{species}_nr_blast_viral.csv',
        viral_peptide_list = 'output/nr_analysis/{species}/{species}_viral_peptides.txt'
    singularity:
        tidyverse_container
    log:
        'output/logs/{species}/blastp_nr_viral_res.log'
    script:
        'src/nr_blastp_viral_analysis.R'

rule blastp_nr:
    input:
        pot_viral_peptides = 'output/viral_blastp/{species}_potential_viral_peptides.faa'
    output:
        blastp_res = 'output/nr_blastp/{species}_blastp.outfmt6'
    params:
        blast_db = 'bin/db/blastdb/nr/nr'
    threads:
        50
    singularity:
        blast_container
    log:
        'output/logs/{species}/nr_blastp.log'
    shell:
        'blastp '
        '-query {input.pot_viral_peptides} '
        '-db {params.blast_db} '
        '-num_threads {threads} '
        '-evalue 1e-05 '
        '-max_target_seqs 1 '
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

##rename script to match rule
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

rule blastp_viral:
    input:
        query = 'data/final_microctonus_assemblies_annotations/{species}.proteins.fa',
        taxid_list = 'output/taxids/species_virus_taxids.txt'
    output:
        blastp_res = 'output/viral_blastp/{species}_blastp.outfmt6'
    params:
        blast_db = 'bin/db/blastdb/nr/nr'
    threads:
        50
    singularity:
        blast_container
    log:
        'output/logs/{species}/viral_blastp.log'
    shell:
        'blastp '
        '-query {input.query} '
        '-db {params.blast_db} '
        '-taxidlist {input.taxid_list} '
        '-num_threads {threads} '
        '-evalue 1e-05 '
        '-max_target_seqs 1 '
        '-outfmt "6 std staxids salltitles" > {output.blastp_res} '
        '2> {log}'

##filter to below genus level as blast only uses species level or below
rule filter_virus_taxids:
    input:
        'output/taxids/virus_taxids.txt'
    output:
        'output/taxids/species_virus_taxids.txt'
    shell:
        'cat {input} | '
        "bin/taxonkit-0.8.0 filter "
        "-L genus "
        '> {output}'

##generate list of viral taxids
rule list_virus_taxids:
    output:
        virus_taxid_list = 'output/taxids/virus_taxids.txt'
    shell:
        'bin/taxonkit-0.8.0 list '
        '--ids 10239 '
        '--indent "" '
        '> {output}'

##to download nr db
##nohup update_blastdb.pl --num_threads 50 --decompress nr &> nohup.out &