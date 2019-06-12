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

blastp_res <- snakemake@input[["blastp_res"]]

########
# MAIN #
########

blastp_res_table <- fread(blastp_res)
fwrite(unique(blastp_res_table[,list(V1)]), snakemake@output[["peptide_hit_ids"]], col.names = FALSE)

#write log
sessionInfo()