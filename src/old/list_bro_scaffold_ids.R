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

exonerate_res <- snakemake@input[["exonerate_res"]]

########
# MAIN #
########

genome_bro_res <- fread(exonerate_res, header=FALSE, sep="")
bro_table1 <- genome_bro_res[,tstrsplit(V1, "M ", fixed=TRUE, keep=1)]
bro_table2 <- bro_table1[,tstrsplit(V1, "vulgar:", fixed=TRUE, keep=2)]
bro_table3 <- bro_table2[,tstrsplit(V1, " ", fixed=TRUE, keep=c(1,2,3,4,5,6,7,8,9,10))]
bro_final_table <- select(bro_table3, -V1, -V5, -V9)
setnames(bro_final_table, old=c("V2", "V3", "V4", "V6", "V7", "V8", "V10"), new=c("peptide_id (query)", "query_start", "query_end", "target", "target_start", "target_end", "raw_score"))
fwrite(bro_final_table, snakemake@output[["bro_final_table"]])

##write list of scaffolds bro mapped to
bro_scaffolds <- (unique(bro_final_table$target))
fwrite((list(bro_scaffolds)), snakemake@output[["bro_scaffold_ids"]])

#write log
sessionInfo()