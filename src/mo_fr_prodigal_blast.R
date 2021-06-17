library(data.table)
library(dplyr)

mo <- fread("output/prodigal/MO/blastp_gff.csv", na.string="")
fr <- fread("output/prodigal/FR/blastp_gff.csv", na.string="")
mo_annots <- subset(mo, !(is.na(annotation)))
fr_annots <- subset(fr, !(is.na(annotation)))

annot_overlap <- intersect(mo_annots$annotation, fr_annots$annotation)
