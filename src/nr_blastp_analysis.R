library("data.table")

##M.hyp results
mhyp_nr_blastp <- fread("output/nr_blastp/Mhyp_blastp.outfmt3")
setnames(mhyp_nr_blastp, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"), new=c("peptide_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))
##extract result with lowest evalue for each peptide
mh_min_evalues <- mhyp_nr_blastp[,.SD[which.min(evalue)], by=peptide_id]

##M.aeth FR results
Ma_FR_nr_blastp <- fread("output/nr_blastp/Maeth_FR_blastp.outfmt3")
setnames(Ma_FR_nr_blastp, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"), new=c("peptide_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))
##extract result with lowest evalue for each peptide
Ma_FR_min_evalues <- Ma_FR_nr_blastp[,.SD[which.min(evalue)], by=peptide_id]

##M.aeth IE results
Ma_IE_nr_blastp <- fread("output/nr_blastp/Maeth_IE_blastp.outfmt3")
setnames(Ma_IE_nr_blastp, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"), new=c("peptide_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))
##extract result with lowest evalue for each peptide
Ma_IE_min_evalues <- Ma_IE_nr_blastp[,.SD[which.min(evalue)], by=peptide_id]

##M.aeth MA results
Ma_MA_nr_blastp <- fread("output/nr_blastp/Maeth_MA_blastp.outfmt3")
setnames(Ma_MA_nr_blastp, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"), new=c("peptide_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))
##extract result with lowest evalue for each peptide
Ma_MA_min_evalues <- Ma_MA_nr_blastp[,.SD[which.min(evalue)], by=peptide_id]
