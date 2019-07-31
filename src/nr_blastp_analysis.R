library("data.table")
library("dplyr")
library("ggplot2")

##M.hyp results
mhyp_nr_blastp <- fread("output/nr_blastp/Mhyp_blastp.outfmt3")
setnames(mhyp_nr_blastp, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"), new=c("peptide_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(mhyp_nr_blastp, peptide_id, evalue, -bit_score)
##extract result with lowest evalue for each peptide - what if multiple rows with lowest min?
Mh_min_evalues <- mhyp_nr_blastp[,.SD[which.min(evalue)], by=peptide_id]

Mh_virus <- dplyr::filter(Mh_min_evalues, grepl('virus', annotation))
fwrite(Mh_virus, "output/viral_nr_blastp_r/Mh/Mh_viral.csv")
Mh_virus_peptide_list <- (Mh_virus$peptide_id)
fwrite(list(Mh_virus_peptide_list), "output/viral_nr_blastp_r/Mh/Mh_viral_peptide_list.txt")
Mh_leptopilina <- dplyr::filter(Mh_min_evalues, grepl('Leptopilina boulardi filamentous virus', annotation))
fwrite(Mh_leptopilina, "output/viral_nr_blastp_r/Mh/Mh_leptopilina.csv")

##plot viral species
viral_species <- fread("output/viral_nr_blastp_r/Mh/Mh_viral_species.csv")
annots_per_species <- viral_species[,length(unique(peptide_id)), by=species]
ggplot(annots_per_species, aes(x=reorder(species, -V1), y=V1))+
  theme(axis.text.x = element_text(angle = 65, hjust = 1, face = "italic")) +
  geom_col()+xlab("Species")+ylab("Number of Annotations")

##filter out baculovirus hits from inital viral blastp and compare to better hits from nrblastp
mhyp_viral_blast <- fread('output/viral_blastp/Mhyp_blastp.outfmt3')
setnames(mhyp_viral_blast, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"), new=c("peptide_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(mhyp_viral_blast, peptide_id, evalue, -bit_score)
##extract result with lowest evalue for each peptide - what if multiple rows with lowest min?
Mh_viral_min_evalues <- mhyp_viral_blast[,.SD[which.min(evalue)], by=peptide_id]
##filter baculo hits from first viral blastp
Mh_viral_blast_baculo <- dplyr::filter(Mh_viral_min_evalues, grepl('baculovirus | granulovirus | nucleopolyhedrovirus', annotation))
##merge with nr blastp hits for baculo hit peptides to compare annots
merged_blastp_baculo_hits <- merge(Mh_viral_blast_baculo, Mh_min_evalues, by.x='peptide_id', by.y='peptide_id', all.x=TRUE, all.y=FALSE)


##M.aeth FR results
Ma_FR_nr_blastp <- fread("output/nr_blastp/Maeth_FR_blastp.outfmt3")
setnames(Ma_FR_nr_blastp, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"), new=c("peptide_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(Ma_FR_nr_blastp, peptide_id, evalue, -bit_score)
##extract result with lowest evalue for each peptide
Ma_FR_min_evalues <- Ma_FR_nr_blastp[,.SD[which.min(evalue)], by=peptide_id]
Ma_FR_virus <- dplyr::filter(Ma_FR_min_evalues, grepl('virus', annotation))
fwrite(Ma_FR_virus, "output/viral_nr_blastp_r/Ma_FR/Ma_FR_viral.csv")

##M.aeth IE results
Ma_IE_nr_blastp <- fread("output/nr_blastp/Maeth_IE_blastp.outfmt3")
setnames(Ma_IE_nr_blastp, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"), new=c("peptide_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(Ma_IE_nr_blastp, peptide_id, evalue, -bit_score)
##extract result with lowest evalue for each peptide
Ma_IE_min_evalues <- Ma_IE_nr_blastp[,.SD[which.min(evalue)], by=peptide_id]
Ma_IE_virus <- dplyr::filter(Ma_IE_min_evalues, grepl('virus', annotation))
fwrite(Ma_IE_virus, "output/viral_nr_blastp_r/Ma_IE/Ma_IE_viral.csv")

##M.aeth MA results
Ma_MA_nr_blastp <- fread("output/nr_blastp/Maeth_MA_blastp.outfmt3")
setnames(Ma_MA_nr_blastp, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"), new=c("peptide_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(Ma_MA_nr_blastp, peptide_id, evalue, -bit_score)
##extract result with lowest evalue for each peptide
Ma_MA_min_evalues <- Ma_MA_nr_blastp[,.SD[which.min(evalue)], by=peptide_id]
Ma_MA_virus <- dplyr::filter(Ma_MA_min_evalues, grepl('virus', annotation))
fwrite(Ma_MA_virus, "output/viral_nr_blastp_r/Ma_MA/Ma_MA_viral.csv")
