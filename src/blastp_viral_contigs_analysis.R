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

viral_contig_blast_res <- snakemake@input[["viral_contig_blast_res"]]
nr_viral_blast_res <- snakemake@input[["nr_viral_blast_res"]]
gff_file <- snakemake@input[["gff_file"]]

########
# MAIN #
########

viral_contig_blast <- fread(viral_contig_blast_res)
nr_viral_blast <- fread(nr_viral_blast_res)
gff <- fread(gff_file, skip=1, header=FALSE)

##viral_contig_blast res
setnames(viral_contig_blast, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14"),
         new=c("peptide_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "taxid", "annotation"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(viral_contig_blast, peptide_id, evalue, -bit_score)
##extract result with lowest evalue for each peptide, sorted by bit-value in case of e-value tie
min_evalues <- viral_contig_blast[,.SD[which.min(evalue)], by=peptide_id]

##join with viral blast
blast_full_table <- full_join(min_evalues, nr_viral_blast)

##set GFF column names
setnames(gff, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9"),
         new=c("contig_id", "source", "type", "start", "end", "score", "strand", "phase", "attributes"))

##filter down to mRNA rows only
gff_mrna <- subset(gff, gff$type=="mRNA")
##split attributes column
gff_mrna$peptide_id <- tstrsplit(gff_mrna$attributes, "=", keep=c(2))
gff_mrna$peptide_id <- tstrsplit(gff_mrna$peptide_id, ";", keep=c(1))
##peptide to contig
peptide_to_contig <- gff_mrna[,c(1,10)]

##blast res with scaffold labels
final_table <- merge(blast_full_table, peptide_to_contig, all.x=TRUE, by="peptide_id")
fwrite(final_table, snakemake@output[["full_blast_table"]])

# write log
sessionInfo()
