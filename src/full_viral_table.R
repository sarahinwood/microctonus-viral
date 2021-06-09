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

interpro_res <- snakemake@input[["interpro_res"]]
blastp_gff_file <- snakemake@input[["blastp_gff_file"]]

########
# MAIN #
########

interpro <- fread(interpro_res)
blastp_gff <- fread(blastp_gff_file)




#write log
sessionInfo()