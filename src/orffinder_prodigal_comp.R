library(rtracklayer)
library(data.table)
library(dplyr)

##orffinder predictions into table
mh_orffinder <- fread("output/orffinder/Mh/Mh_orffinder_coords.out", header=FALSE, sep = ":")
mh_orffinder$V1 <- tstrsplit(mh_orffinder$V1, ">lcl\\|", keep=c(2))
mh_orffinder$ID <- tstrsplit(mh_orffinder$V1, "_", keep=c(1))
mh_orffinder$scaffold_no <- tstrsplit(mh_orffinder$V1, "_", keep=c(3))
mh_orffinder$seqid <- paste("scaffold_", mh_orffinder$scaffold_no, sep="")
mh_orffinder$end <- tstrsplit(mh_orffinder$V3, " ", keep=c(1))
mh_orffinder$partial <- tstrsplit(mh_orffinder$V3, " ", keep=c(5))

##organise table
mh_orffinder_table <- mh_orffinder[,c(6,2,7,8,4)]
setnames(mh_orffinder_table, old=c("V2"), new=c("start"))
mh_orffinder_table$source <- paste("ORFfinder_v0.4.3")
mh_orffinder_table$type <- paste("CDS")
mh_orffinder_table$end <- as.numeric(mh_orffinder_table$end)
##if start less than end, +, otherwise -
mh_orffinder_table$strand <- ifelse((mh_orffinder_table$start<mh_orffinder_table$end), "+", "-")
final_mh_orffinder_table <- mh_orffinder_table[,c(1,6,7,8,2,3,5,4)]

##compare with prodigal
prodigal_gff <- data.table(data.frame(readGFF("output/prodigal/Mh/gene_predictions.gff")))
viral_predictions <- full_join(final_mh_orffinder_table, prodigal_gff)


##getorf
mh_getorf <- fread("output/getorf/Mh/Mh_getorf_coords.out", header=FALSE, sep = " ", fill=TRUE)
mh_getorf$ID <- tstrsplit(mh_getorf$V1, ">", keep=c(2))
mh_getorf$seqid <- tstrsplit(mh_getorf$ID, "_", keep=c(2))
mh_getorf$seqid <- paste("scaffold_", mh_getorf$seqid, sep="")
mh_getorf$start <- tstrsplit(mh_getorf$V2, "\\[", keep=c(2))
mh_getorf$end <- tstrsplit(mh_getorf$V4, "\\]", keep=c(1))
mh_getorf$end <- as.numeric(unlist(mh_getorf$end))
mh_getorf$start <- lapply(mh_getorf$start, as.numeric)
mh_getorf$start <- as.numeric(unlist(mh_getorf$start))
mh_getorf$strand <- tstrsplit(mh_getorf$V5, "\\(", keep=c(2))
mh_getorf$source <- "getorf_v6.6.0.0"
mh_getorf_table <- mh_getorf[,c(7,8,9,10,11,12)]


all_viral_predictions <- full_join(viral_predictions, mh_getorf_table)
fwrite(all_viral_predictions, "output/getorf/Mh/orffinder_getorf_prodigal.csv")
