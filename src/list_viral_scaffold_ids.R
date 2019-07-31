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

viral_exonerate_res <- snakemake@input[["viral_exonerate_res"]]

########
# MAIN #
########

##make tidy table of exonerate res
viral_exonerate_res <- fread(viral_exonerate_res, header=FALSE, sep="")
viral_table1 <- viral_exonerate_res[,tstrsplit(V1, "M ", fixed=TRUE, keep=1)]
viral_table2 <- viral_table1[,tstrsplit(V1, "vulgar:", fixed=TRUE, keep=2)]
viral_table3 <- viral_table2[,tstrsplit(V1, " ", fixed=TRUE, keep=c(1,2,3,4,5,6,7,8,9,10))]
viral_table4 <- select(viral_table3, -V1, -V5, -V9)
setnames(viral_table4, old=c("V2", "V3", "V4", "V6", "V7", "V8", "V10"), new=c("peptide_id (query)", "query_start", "query_end", "target", "target_start", "target_end", "raw_score"))
viral_exonerate_table <- viral_table4[,.SD[which.max(raw_score)], by=`peptide_id (query)`]
fwrite(viral_exonerate_table, snakemake@output[["viral_exonerate_table"]])
##make list of scaffolds viral genes map onto
viral_scaffold_ids <- viral_exonerate_table$target
fwrite(list(viral_scaffold_ids), snakemake@output[["viral_scaffold_ids"]])

#write log
sessionInfo()