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

exonerate_bro_scaffolds_res <- snakemake@input[["exonerate_bro_scaffolds_res"]]
bro_final_table <- snakemake@input[["bro_final_table"]]

########
# MAIN #
########

##exonerate to see what peptides map to bro scaffolds
##to then blast peptides and see what other genes are around bro genes in these scaffolds
bro_scaffold_peptides_exonerate <- fread(exonerate_bro_scaffolds_res, header=FALSE, sep="")
bro_scaffold_peptides_exonerate <- bro_scaffold_peptides_exonerate[-c(1,2,67),]
bspe1 <- bro_scaffold_peptides_exonerate[,tstrsplit(V1, "M ", fixed=TRUE, keep=1)]
bspe2 <- bspe1[,tstrsplit(V1, "vulgar:", fixed=TRUE, keep=2)]
bspe3 <- bspe2[,tstrsplit(V1, " ", fixed=TRUE, keep=c(1,2,3,4,5,6,7,8,9,10))]
bspe_final <- select(bspe3, -V1, -V5, -V9)
setnames(bspe_final, old=c("V2", "V3", "V4", "V6", "V7", "V8", "V10"), new=c("peptide_id (query)", "query_start", "query_end", "target", "target_start", "target_end", "raw_score"))
fwrite(bspe_final, snakemake@output[["bspe_final_table"]])

list_peptides_on_bro_scaffolds <- data.table(unique(bspe_final$`peptide_id (query)`))
bro_final_table <- fread(bro_final_table)
list_bro_peptides <- data.table(unique(bro_final_table$`peptide_id (query)`))

##filter out peptides on bro scaffolds that are NOT bro peptides
peptides_not_bro <- list_peptides_on_bro_scaffolds$V1[!(list_peptides_on_bro_scaffolds$V1 %in% list_bro_peptides$V1)]
fwrite(list(peptides_not_bro), snakemake@output[["peptide_list"]])

