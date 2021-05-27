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

exonerate_viral_scaffolds_res <- snakemake@input[["exonerate_viral_scaffolds_res"]]
viral_exonerate_table <- snakemake@input[["viral_exonerate_table"]]

########
# MAIN #
########

##exonerate to see what peptides map to viral scaffolds
##to then blast peptides and see what other genes are around viral genes in these scaffolds
viral_scaffold_peptides_exonerate <- fread(exonerate_viral_scaffolds_res, header=FALSE, sep="")
viral_scaffold_peptides_exonerate <- viral_scaffold_peptides_exonerate[-c(1,2,123),]
vspe1 <- viral_scaffold_peptides_exonerate[,tstrsplit(V1, "M ", fixed=TRUE, keep=1)]
vspe2 <- vspe1[,tstrsplit(V1, "vulgar:", fixed=TRUE, keep=2)]
vspe3 <- vspe2[,tstrsplit(V1, " ", fixed=TRUE, keep=c(1,2,3,4,5,6,7,8,9,10))]
vspe_final <- select(vspe3, -V1, -V5, -V9)
setnames(vspe_final, old=c("V2", "V3", "V4", "V6", "V7", "V8", "V10"), new=c("peptide_id (query)", "query_start", "query_end", "target", "target_start", "target_end", "raw_score"))
fwrite(vspe_final, snakemake@output[["vspe_final_table"]])

##both viral and other peptides on viral-containing scaffolds
list_peptides_on_viral_scaffolds <- data.table(unique(vspe_final$`peptide_id (query)`))
fwrite(list(list_peptides_on_viral_scaffolds$V1), snakemake @output[["all_peptides_list"]])
viral_exonerate_table <- fread(viral_exonerate_table)
list_viral_peptides <- data.table(unique(viral_exonerate_table$`peptide_id (query)`))

##filter out peptides on viral scaffolds that are NOT viral peptides
peptides_not_viral <- list_peptides_on_viral_scaffolds$V1[!(list_peptides_on_viral_scaffolds$V1 %in% list_viral_peptides$V1)]
fwrite(list(peptides_not_viral), snakemake@output[["non_viral_peptide_list"]])

