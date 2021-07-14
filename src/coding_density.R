library(rtracklayer)
library(data.table)

#Mh viral contigs bbstats length
mh_viral_contig_length <- 92650

prodigal_gff <- data.frame(readGFF("output/prodigal/Mh/gene_predictions.gff"))
prodigal_gff$length <- (prodigal_gff$end)-(prodigal_gff$start)+1
sum(prodigal_gff$length)/mh_viral_contig_length

sum(!(prodigal_gff$partial=="00"))

#n genes - 98 (22 incomplete at one or both edges though)
#length ranges 108-3867
#start MET and end in STOP - 5 genes are complete but don't start with ATG - viruses only ATG
#how many of these genes don't even have homology to LbFV?
#91.7% coding density in 92.6kb total (LbFV 80.0% in 111kb)
#33.8% GC content - 66.2% AT (LbFV 78.7% AT)

#repeat content at all? particularly on edges of contigs?

##spread of genes on either strand


##orffinder
##118 genes
