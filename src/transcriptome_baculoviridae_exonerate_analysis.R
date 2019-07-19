library('data.table')
library('dplyr')

baculoviridae_exonerate <- fread('output/mh_exonerate/vul_baculo_exonerate.out',header=FALSE, sep="")
bac_exonerate_res_table1 <- baculoviridae_exonerate[,tstrsplit(V1, " M ", fixed=TRUE, keep=1)]
bac_exonerate_res_table2 <- bac_exonerate_res_table1[,tstrsplit(V1, "vulgar: ", fixed=TRUE, keep=2)]
bac_exonerate_res_table3 <- bac_exonerate_res_table2[,tstrsplit(V1, " ", fixed=TRUE, keep=c(1,2,3,4,5,6,7,8,9))]

bac_exonerate_res_final_res_table <- select(bac_exonerate_res_table3, -V4, -V8)
setnames(bac_exonerate_res_final_res_table, old=c("V1", "V2", "V3", "V5", "V6", "V7", "V9"), new=c("trinity_id", "query_start", "query_end", "target", "target_start", "target_end", "raw_score"))
bac_exon_max_scores <- bac_exonerate_res_final_res_table[,.SD[which.max(raw_score)], by="trinity_id"]
write.csv(bac_exon_max_scores, "output/mh_exonerate/baculoviridae_exonerate_max_scores.csv")


genome_bro_res <- fread('output/mh_exonerate/mh_genome_bro_exonerate.out',header=FALSE, sep="")
bro_table1 <- genome_bro_res[,tstrsplit(V1, "M ", fixed=TRUE, keep=1)]
bro_table2 <- bro_table1[,tstrsplit(V1, "vulgar:", fixed=TRUE, keep=2)]
bro_table3 <- bro_table2[,tstrsplit(V1, " ", fixed=TRUE, keep=c(1,2,3,4,5,6,7,8,9,10))]
bro_final_table <- select(bro_table3, -V1, -V5, -V9)
setnames(bro_final_table, old=c("V2", "V3", "V4", "V6", "V7", "V8", "V10"), new=c("peptide_id (query)", "query_start", "query_end", "target", "target_start", "target_end", "raw_score"))
write.csv(bro_final_table, "output/mh_exonerate/bro_genome_exonerate_full_res.csv")

bro_scaffolds <- (unique(bro_final_table$target))
fwrite((list(bro_scaffolds)), "output/mh_exonerate/bro_scaffold_ids.txt")
