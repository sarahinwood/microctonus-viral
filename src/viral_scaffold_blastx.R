library("data.table")
library("dplyr")

viral_scaffold_blastx <- fread("output/mh_exonerate/genome_viral_scaffolds/blastx/blastx_viral_scaffolds.outfmt3")
setnames(viral_scaffold_blastx, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"), new=c("scaffold_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))
setorder(viral_scaffold_blastx, scaffold_id, query_start)
fwrite(viral_scaffold_blastx, "output/mh_exonerate/genome_viral_scaffolds/blastx/blastx_viral_scaffolds.csv")

##filter out any LbFV hits
scaff_bx_leptopilina <- dplyr::filter(viral_scaffold_blastx, grepl('Leptopilina boulardi filamentous virus', annotation))
fwrite(scaff_bx_leptopilina, "output/mh_exonerate/genome_viral_scaffolds/blastx/LbFV_blastx_annots.csv")

##scaff.8624 & scaff80 - are there any viral hits for these scaffolds?
scaff_8624_80 <- dplyr::filter(viral_scaffold_blastx, grepl('Scaffold8624|Scaffold80', scaffold_id))
