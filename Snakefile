#!/usr/bin/env python3

import peppy

#############
# FUNCTIONS #
#############

def get_peptide_dbs(wildcards):
    input_keys = ['peptide_db']
    my_pep = pep.get_sample(wildcards.sample).to_dict()
    return {k: my_pep[k] for k in input_keys}

##needed to get BUSCO running in new folder
def resolve_path(x):
    return(str(pathlib2.Path(x).resolve(strict=False)))

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
busco_container = 'docker://ezlabgva/busco:v4.0.2_cv1'
tidyverse_container = 'docker://rocker/tidyverse'
bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'

#########
# RULES #
#########

rule target:
    input:
        ##recip blast -> all peptides on viral scaffolds blast results
        expand('output/nr_analysis_{species}/{species}_blastp_peptides_viral_scaffolds.outfmt6', species=all_species),
        ##interproscan results
        expand('output/nr_analysis_{species}/{species}_interproscan.tsv', species=all_species),
        ##busco on non-hic genomes
        expand('output/busco/{not_hic_species}/short_summary.specific.endopterygota_odb10.{not_hic_species}.txt', filter=['MO', 'IR']),
        ##GC content plots
        expand('output/bbstats/{hic_species}/all_gc_boxplot.pdf', hic_species=['Mh', 'FR']),
        expand('output/bbstats/{not_hic_species}/all_gc_boxplot.pdf', not_hic_species=['MO', 'IR'])

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

###########################
## GC content of contigs ##
###########################

rule plot_GC_content_not_hic:
    input:
        gc = 'output/bb_stats/{not_hic_species}/gc.txt',
        gc_hist = 'output/bb_stats/{not_hic_species}/gc_hist.out',
        viral_contig_list = 'output/nr_analysis_{not_hic_species}/{not_hic_species}_viral_scaffold_counts.csv',
        busco = 'output/busco/{not_hic_species}/full_table_{not_hic_species}.tsv'
    output:
        hic_only = 'output/bbstats/{not_hic_species}/hic_gc_boxplot.pdf',
        all_contigs = 'output/bbstats/{not_hic_species}/all_gc_boxplot.pdf',
        gc_histogram = 'output/bbstats/{not_hic_species}/gc_histogram.pdf'
    log:
        'output/logs/bbstats/{not_hic_species}_plot_gc.log'
    script:
        'src/plot_GC_content_not_hic.R'

rule plot_GC_content_hic:
    input:
        gc = 'output/bb_stats/{hic_species}/gc.txt',
        gc_hist = 'output/bb_stats/{hic_species}/gc_hist.out',
        viral_contig_list = 'output/nr_analysis_{hic_species}/{hic_species}_viral_scaffold_counts.csv',
        hic_scaffold_list = 'data/hi-c_genomes/{species}_hic_scaffold_ids.txt'
    output:
        hic_only = 'output/bbstats/{hic_species}/hic_gc_boxplot.pdf',
        all_contigs = 'output/bbstats/{hic_species}/all_gc_boxplot.pdf',
        gc_histogram = 'output/bbstats/{hic_species}/gc_histogram.pdf'
    log:
        'output/logs/bbstats/{hic_species}_plot_gc.log'
    script:
        'src/plot_GC_content_hic.R'

##for two assemblies without hi-c need busco to compare busco containing to viral
rule busco_non_hic_genomes:
    input:
        genome = 'data/final_genomes/{not_hic_species}.fa'
    output:
        'output/busco/{not_hic_species}/full_table_{not_hic_species}.tsv'
    log:
        str(pathlib2.Path(resolve_path('output/logs/'),
                            'busco_{not_hic_species}.log'))
    params:
        wd = 'output/busco',
        genome = lambda wildcards, input: resolve_path(input.genome)
    singularity:
        busco_container
    threads:
        20
    shell:
        'cd {params.wd} || exit 1 ; '
        'busco '
        '--force '
        '--in {params.genome} '
        '--out {wildcards.not_hic_species} '
        '--lineage endopterygota_odb10 '
        '--cpu {threads} '
        '--augustus_species tribolium '
        '--mode transcriptome '
        '-f '
        '&> {log} '

rule bb_stats:
    input:
        genome = 'data/final_genomes/{species}.fa'
    output:
        gc = 'output/bb_stats/{species}/gc.txt',
        stats = 'output/bb_stats/{species}/bb_stats.out',
        gc_hist = 'output/bb_stats/{species}/gc_hist.out'
    log:
        'output/logs/bbstats/{species}.log'
    singularity:
        bbduk_container
    shell:
        'stats.sh '
        'in={input.genome} '
        'out={output.stats} '
        'gc={output.gc} '
        'gcformat=4 '
        'gchist={output.gc_hist} '

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

##removes scaffolds that are main-genome from hic to avoid blasting every gene on them
rule list_peptides_on_viral_scaffolds:
    input:
        viral_peptide_list = 'output/nr_analysis_{species}/{species}_viral_peptides.txt',
        gff = 'data/final_microctonus_assemblies_annotations/{species}.gff3',
        hic_scaffold_list = 'data/hi-c_genomes/{species}_hic_scaffold_ids.txt'
    output:
        peptides_viral_scaffolds = 'output/nr_analysis_{species}/{species}_peptides_viral_scaffolds.txt',
        scaffold_counts = 'output/nr_analysis_{species}/{species}_viral_scaffold_counts.csv'
    singularity:
        tidyverse_container
    log:
        'output/logs/{species}/link_viral_peptides_to_scaffolds.log'
    script:
        'src/list_peptides_on_viral_scaffolds.R'

#############################
## reciprocal viral blastp ##
#############################

rule nr_blastp_viral_analysis:
    input:
        nr_blastp_res = 'output/nr_blastp/{species}_blastp.outfmt6',
        taxid_list = 'output/taxids/species_virus_taxids.txt'
    output:
        blastp_viral_res = 'output/nr_analysis_{species}/{species}_nr_blast_viral.csv',
        viral_peptide_list = 'output/nr_analysis_{species}/{species}_viral_peptides.txt'
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