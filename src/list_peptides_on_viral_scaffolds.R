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

###########
# GLOBALS #
###########

viral_peptide_list <- snakemake@input[["viral_peptide_list"]]
gff_file <- snakemake@input[["gff"]]

########
# MAIN #
########

viral_peptides <- fread(viral_peptide_list, header=FALSE)
gff <- fread(gff_file)
##set column names
setnames(gff, old=c("##gff-version 3", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9"), new=c("scaffold_id", "source", "type", "start", "end", "score", "strand", "phase", "attributes"))

##filter down to mRNA rows only
gff_mrna <- subset(gff, gff$type=="mRNA")
##split attributes column
gff_mrna$peptide_id <- tstrsplit(gff_mrna$attributes, "=", keep=c(2))
gff_mrna$peptide_id <- tstrsplit(gff_mrna$peptide_id, ";", keep=c(1))

##gff table with viral peptide IDs
viral_peptide_gff <- subset(gff_mrna, peptide_id %in% viral_peptides$V1)
##list scaffolds with viral peptides
viral_scaffolds <- data.table(unique(viral_peptide_gff$scaffold_id))

##subset mRNA gff to get all peptides on scaffolds with viral peptides
viral_scaffold_peptides <- subset(gff_mrna, scaffold_id %in% viral_scaffolds$V1)
fwrite(list(viral_scaffold_peptides$peptide_id), snakemake@output[["peptides_viral_scaffolds"]])

