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

viral_peptide_list <- snakemake@input[["viral_peptide_list"]]
gff_file <- snakemake@input[["gff_file"]]
hic_scaffold_list <- snakemake@input[["hic_scaffold_list"]]

########
# MAIN #
########

viral_peptides <- fread(viral_peptide_list, header=FALSE)
hic_scaffold_ids <- fread(hic_scaffold_list, header=FALSE)
##note that contig names are changed from final assembly file but are same contigs
gff <- fread(gff_file, skip=1, header=FALSE)
##set column names
setnames(gff, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9"),
	new=c("contig_id", "source", "type", "start", "end", "score", "strand", "phase", "attributes"))

##filter down to mRNA rows only
gff_mrna <- subset(gff, gff$type=="mRNA")
##split attributes column
gff_mrna$peptide_id <- tstrsplit(gff_mrna$attributes, "=", keep=c(2))
gff_mrna$peptide_id <- tstrsplit(gff_mrna$peptide_id, ";", keep=c(1))

##gff table with viral peptide IDs
viral_peptide_gff <- subset(gff_mrna, peptide_id %in% viral_peptides$V1)
##list contigs with viral peptides
viral_contigs <- data.table(unique(viral_peptide_gff$contig_id))

##subset mRNA gff to get all peptides on contigs with viral peptides
viral_contig_peptides <- subset(gff_mrna, contig_id %in% viral_contigs$V1)
##count how many peptides from each contig
contig_counts <- viral_contig_peptides %>% count(contig_id)

##filter out main genome hic contigs - too many peptides to blast
viral_contigs_not_hic <- subset(viral_contig_peptides, !(contig_id %in% hic_scaffold_ids$V1))
##remove peptides with viral BlastP hits - have already searched them
viral_contigs_peptides <- subset(viral_contigs_not_hic, !(peptide_id %in% viral_peptides$V1))

##output
fwrite(contig_counts, snakemake@output[["contig_counts"]])
fwrite(list(unique(viral_contigs_not_hic$contig_id)), snakemake@output[["viral_contig_ids"]])
##write list of peptides
fwrite(list(viral_contigs_peptides$peptide_id), snakemake@output[["peptides_viral_contigs"]])

# write log
sessionInfo()