#!/usr/bin/env Rscript

#######
# LOG #
#######

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

#############
# LIBRARIES #
#############

library(data.table)
library(dplyr)

###########
# GLOBALS #
###########

nr_blastp_res <- snakemake@input[["nr_blastp_res"]]
taxid_list <- snakemake@input[['taxid_list']]

########
# MAIN #
########

nr_blastp <- fread(nr_blastp_res)
viral_taxids <- fread(taxid_list, header=FALSE)

setnames(nr_blastp, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14"),
         new=c("peptide_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "taxid", "annotation"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(nr_blastp, peptide_id, evalue, -bit_score)
##extract result with lowest evalue for each peptide, sorted by bit-value in case of e-value tie
min_evalues <- nr_blastp[,.SD[which.min(evalue)], by=peptide_id]
##extract hits to viruses
virus_annots <- dplyr::filter(nr_blastp, taxid %in% viral_taxids$V1)
fwrite(virus_annots, snakemake@output[["blastp_viral_res"]])

##list of peptides with viral hits
virus_peptide_list <- (virus_annots$peptide_id)
fwrite(list(virus_peptide_list), snakemake@output[["viral_peptide_list"]])

#write log
sessionInfo()
