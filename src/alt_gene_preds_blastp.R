library(data.table)
library(dplyr)

############
## getORF ## 23 LbFV, 6 DaFV, 48 total blast hits, 131 total predictions
############

getorf_blastp <- fread("output/getorf/Mh/getorf_blastp.outfmt6")
setnames(getorf_blastp, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14"),
         new=c("ID", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "taxid", "annotation"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(getorf_blastp, ID, evalue, -bit_score)
##extract result with lowest evalue for each peptide, sorted by bit-value in case of e-value tie
getorf_min_evalues <- getorf_blastp[,.SD[which.min(evalue)], by=ID]
full_preds_blast <- merge(getorf_min_evalues, mh_getorf_table, by="ID", all.y=TRUE)
full_preds_blast$orf_no <- tstrsplit(full_preds_blast$ID, "_", keep=c(3))
fwrite(full_preds_blast, "output/getorf/Mh/getorf_table_blastp.csv")

###############
## ORFfinder ## 21 LbFV, 6 DaFV, 41 total blast hits, 118 total predictions
###############

orffinder_blastp <- fread("output/orffinder/Mh/orffinder_blastp.outfmt6")
setnames(orffinder_blastp, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14"),
         new=c("id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "taxid", "annotation"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(orffinder_blastp, id, evalue, -bit_score)
##extract result with lowest evalue for each peptide, sorted by bit-value in case of e-value tie
orffinder_min_evalues <- orffinder_blastp[,.SD[which.min(evalue)], by=id]
fwrite(orffinder_min_evalues, "output/orffinder/Mh/orffinder_blastp.csv")

###############
## GeneMarkS ## 23 LbFV, 6 DaFV, 46 total blast hits, 73 total predictions
###############

genemarks_blastp <- fread("output/genemarks/genemarks_blastp.outfmt6")
setnames(genemarks_blastp, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14"),
         new=c("ID", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "taxid", "annotation"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(genemarks_blastp, ID, evalue, -bit_score)
##extract result with lowest evalue for each peptide, sorted by bit-value in case of e-value tie
genemarks_min_evalues <- genemarks_blastp[,.SD[which.min(evalue)], by=ID]
