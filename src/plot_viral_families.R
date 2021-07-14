library(data.table)
library(dplyr)
library(ggplot2)
library(viridis)

blast_res_files <- list.files(path="output", pattern="*_nr_blast_viral_plot.csv",
                              full.names=TRUE, recursive=TRUE)
##read into single table for plotting
tables <- lapply(blast_res_files, fread)
full_table <- do.call(rbind, tables)

## DNA ##
dna_full_table <- subset(full_table, viral_genome == "DNA")
##count no. of each viral family
dna_viral_families <- count(dna_full_table, virus_family, species, sort=TRUE)
##reorder species
dna_viral_families$species <- factor(dna_viral_families$species, levels=c("M. aethiopoides Irish" ,"M. aethiopoides French", "M. aethiopoides Moroccan", "M. hyperodae"))

##plot
dna_plot <- ggplot(dna_viral_families, aes(x=species, y=n, fill=reorder(virus_family, n)))+
  scale_fill_viridis(discrete=TRUE, direction=-1)+
  geom_col()+
  theme_bw()+
  xlab("Species")+
  ylab("Number of peptides with viral hits")+
  labs(fill="DNA virus family")+
  guides(fill = guide_legend(reverse = TRUE))+
  coord_flip()+
  theme(axis.text.y = element_text(hjust = 1, face = "italic"))+
  ylim(0,40)

##need to save as svg and edit in inkscape to un-italicise Irish, French, Moroccan
ggsave(file="output/dna_virus_families.svg", plot=dna_plot, width=8, height=3)

## RNA ##
rna_full_table <- subset(full_table, viral_genome == "RNA")
##count no. of each viral family
rna_viral_families <- count(rna_full_table, virus_family, species, sort=TRUE)
##reorder species
rna_viral_families$species <- factor(rna_viral_families$species, levels=c("M. aethiopoides Irish" ,"M. aethiopoides French", "M. aethiopoides Moroccan", "M. hyperodae"))


##plot
rna_plot <- ggplot(rna_viral_families, aes(x=species, y=n, fill=reorder(virus_family, n)))+
  scale_fill_viridis(discrete=TRUE, direction=-1)+
  geom_col()+
  theme_bw()+
  xlab("Species")+
  ylab("Number of peptides with viral hits")+
  labs(fill="RNA virus family")+
  guides(fill = guide_legend(reverse = TRUE))+
  coord_flip()+
  theme(axis.text.y = element_text(hjust = 1, face = "italic"))+
  ylim(0,40)

##need to save as svg and edit in inkscape to un-italicise Irish, French, Moroccan
ggsave(file="output/rna_virus_families.svg", plot=rna_plot, width=8, height=3)
