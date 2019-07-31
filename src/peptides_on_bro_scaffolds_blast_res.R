library("data.table")
library("VennDiagram")

##blast hits for peptides on same contigs as bro genes
bscaffold_blastp <- fread("output/mh_exonerate/genome_bro_scaffolds/bro_scaffold_peptides_blastp.outfmt3")
setnames(bscaffold_blastp, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"), new=c("peptide_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(bscaffold_blastp, peptide_id, evalue, -bit_score)
##extract result with lowest evalue for each peptide - what if multiple rows with lowest min?
bscaffold_min_evalues <- bscaffold_blastp[,.SD[which.min(evalue)], by=peptide_id]
bscaffold_peptide_ids <- bscaffold_min_evalues$peptide_id

bro_peptides <- fread('output/viral_nr_blastp_r/Mh/mh_bro_peptide_ids.txt', header=FALSE)
bro_peptide_ids <- bro_peptides$V1
viral_recip_blast_peptides <- fread('output/viral_nr_blastp_r/Mh/Mh_viral.csv')
viral_recip_ids <- viral_recip_blast_peptides$peptide_id
leptopilina <- fread("output/viral_nr_blastp_r/Mh/Mh_leptopilina.csv")
leptopilina_ids <- leptopilina$peptide_id

#Draw Venn Diagram
Set1 <- RColorBrewer::brewer.pal(3, "Set1")
vd <- venn.diagram(x = list("Bro peptides"=bro_peptide_ids, "same scaffold as Bro"=bscaffold_peptide_ids, "Viral hit in recip-blastp"=viral_recip_ids), filename=NULL, fill=Set1, alpha=0.5, cex = 1, cat.cex=1, lwd=1)
grid.newpage()
grid.draw(vd)

##leptopilina virus venn
vd2 <- venn.diagram(x = list("Leptopilina boulardi filamentous virus annot"=leptopilina_ids, "same scaffold as Bro"=bscaffold_peptide_ids, "Viral hit in recip-blastp"=viral_recip_ids), filename=NULL, fill=Set1, alpha=0.5, cex = 1, cat.cex=1, lwd=1)
grid.newpage()
grid.draw(vd2)

tidy_exonerate_res <- fread("output/mh_exonerate/genome_bro_scaffolds/exonerate_table.csv")
max_score_exonerate_res <- tidy_exonerate_res[,.SD[which.max(raw_score)], by="peptide_id (query)"]
exonerate_bscaffold_blast_res <- merge(max_score_exonerate_res, bscaffold_min_evalues, by.x="peptide_id (query)", by.y="peptide_id", all.x=TRUE, all.y=TRUE)
recip_blast_res <- fread("output/viral_nr_blastp_r/Mh/Mh_viral.csv")
blast_scaffold_annots <- merge(exonerate_bscaffold_blast_res, recip_blast_res, by.x="peptide_id (query)", by.y="peptide_id", all.x=TRUE, all.y=FALSE)
fwrite(blast_scaffold_annots, "output/mh_exonerate/genome_bro_scaffolds/bro_scaffold_peptide_annots.csv")
