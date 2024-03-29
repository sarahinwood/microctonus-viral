library(data.table)
library(dplyr)
library(ggplot2)
library(viridis)

blast_res_files <- list.files(path="output", pattern="*blastp_gff_plot.csv",
                              full.names=TRUE, recursive=TRUE)
##read into single table for plotting
tables <- lapply(blast_res_files, fread)
full_table_no <- do.call(rbind, tables)
#outliers <- list("scaffold_90", "scaffold_995", "scf1605")
#full_table_no <- subset(full_table, !(scaffold_id %in% outliers))
##count no. of each viral family
families <- count(full_table_no, family, species, euk_or_pro, sort=TRUE)

##reorder
families$species <- factor(families$species, levels=c("M. aethiopoides French", "M. aethiopoides Moroccan", "M. hyperodae"))

eukaryotic <- subset(families, euk_or_pro=="Eukaryotic")
prokaryotic <- subset(families, euk_or_pro=="Prokaryotic")

##plot
plot_euk <- ggplot(eukaryotic, aes(x=species, y=n, fill=reorder(family, n)))+
  scale_fill_viridis(discrete=TRUE, direction=-1)+
  geom_col()+
  theme_bw()+
  xlab("Species")+
  ylab("Number of peptides with BlastP hits")+
  labs(fill="Eukaryotic families")+
  guides(fill = guide_legend(reverse = TRUE))+
  coord_flip()+
  theme(axis.text.y = element_text(hjust = 1, face = "italic"))+ylim(0,50)

##need to save as svg and edit in inkscape to un-italicise French, Moroccan
ggsave(file="output/prodigal/eukaryotic_families.svg", plot=plot_euk, width=8, height=4)

##plot
plot_pro <- ggplot(prokaryotic, aes(x=species, y=n, fill=reorder(family, n)))+
  scale_fill_viridis(discrete=TRUE, direction=-1)+
  geom_col()+
  theme_bw()+
  xlab("Species")+
  ylab("Number of peptides with BlastP hits")+
  labs(fill="Prokaryotic families")+
  guides(fill = guide_legend(reverse = TRUE))+
  coord_flip()+
  theme(axis.text.y = element_text(hjust = 1, face = "italic"))+ylim(0,50)

##need to save as svg and edit in inkscape to un-italicise French, Moroccan
ggsave(file="output/prodigal/prokaryotic_families.svg", plot=plot_pro, width=8, height=4)
