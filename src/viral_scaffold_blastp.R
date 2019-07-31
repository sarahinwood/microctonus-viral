library("data.table")
library("dplyr")
library("ggplot2")

viral_scaffold_nr_blastp <- fread("output/mh_exonerate/genome_viral_scaffolds/viral_scaffold_peptides_blastp.outfmt3")
setnames(viral_scaffold_nr_blastp, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"), new=c("peptide_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(viral_scaffold_nr_blastp, peptide_id, evalue, -bit_score)
fwrite(viral_scaffold_nr_blastp, "output/mh_exonerate/genome_viral_scaffolds/blast_all_res.csv")
##extract result with lowest evalue for each peptide - what if multiple rows with lowest min?
blast_min_evalues <- viral_scaffold_nr_blastp[,.SD[which.min(evalue)], by=peptide_id]
fwrite(blast_min_evalues, "output/mh_exonerate/genome_viral_scaffolds/blast_min_evalues.csv")

##read in exonerate res so we can see scaffold no. with annot
exonerate <- fread("output/mh_exonerate/genome_viral_scaffolds/exonerate_table.csv")
exonerate_max_scores <- exonerate[,.SD[which.max(raw_score)], by=`peptide_id (query)`]

##read in viral res from recip-blast to merge with this also
viral_recip <- fread("output/viral_nr_blastp_r/Mh/Mh_viral.csv")

##merge exonerate with blast results
exonerate_viral_annots <- merge(exonerate_max_scores, viral_recip, by.x="peptide_id (query)", by.y="peptide_id", all=TRUE)
exonerate_all_annots <- merge (exonerate_viral_annots, blast_min_evalues, by.x="peptide_id (query)", by.y="peptide_id", all=TRUE)
fwrite(exonerate_all_annots, "output/mh_exonerate/genome_viral_scaffolds/exonerate_w_annots.csv")

peptide_to_scaffold <- select(exonerate_all_annots, "peptide_id (query)", "target")

exonerate_taxa1 <- exonerate_all_annots[,tstrsplit(annotation.x, "<>", fixed=TRUE, keep=1), by="peptide_id (query)"]
exonerate_taxa2 <- exonerate_taxa1[,tstrsplit(V1, "]", fixed=TRUE), by="peptide_id (query)"]

exonerate_taxa1_2 <- exonerate_all_annots[,tstrsplit(annotation.y, "<>", fixed=TRUE, keep=1), by="peptide_id (query)"]
exonerate_taxa2_2 <- exonerate_taxa1_2[,tstrsplit(V1, "]", fixed=TRUE), by="peptide_id (query)"]

interpro_res <- fread('output/mh_exonerate/genome_viral_scaffolds/interproscan/all_viral_scaffold_peptides.fasta.tsv', header=FALSE, fill=TRUE)
setnames(interpro_res, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14"), new=c("peptide_id", "Sequence_MD5_digest", "Sequence_Length", "Analysis", "Signature_Accession", "Signature_Description", "Start_Location", "Stop_Location", "Score", "Status", "Date", "Interpro_Annotation_Accession", "Interpro_Annotation_Description", "GO_Annotation"))
interpro_scaffold <- merge (interpro_res, peptide_to_scaffold, by.x="peptide_id", by.y="peptide_id (query)", all.x=TRUE)
fwrite(interpro_scaffold, "output/mh_exonerate/genome_viral_scaffolds/interproscan/interpro_res_scaffolds.csv")

taxa <- fread("output/mh_exonerate/genome_viral_scaffolds/taxa.csv")
annots_per_taxa <- taxa[,length(unique(`peptide_id (query)`)), by=V1.x]
ggplot(annots_per_taxa, aes(x=reorder(V1.x, -V1), y=V1))+
  theme(axis.text.x = element_text(angle = 65, hjust = 1, face = "italic")) +
  geom_col()+xlab("Genera")+ylab("Number of BlastX Annotations")

##see if some of the peptides on viral scaffolds had viral hit in first blastp
viral_blast <- fread("output/viral_blastp/Mhyp_blastp.outfmt3")
v_min_e <- viral_blast[,.SD[which.min(V11)], by=V1]
v <- merge(blast_min_evalues, v_min_e, by.x="peptide_id", by.y="V1", all.x=TRUE, all.y=FALSE)
