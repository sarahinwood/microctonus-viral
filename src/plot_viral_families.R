library(data.table)
library(dplyr)
library(ggplot2)
library(viridis)

blast_res_files <- list.files(path="output", pattern="*_nr_blast_viral_plot.csv",
                              full.names=TRUE, recursive=TRUE)
##read into single table for plotting
tables <- lapply(blast_res_files, fread)
full_table <- do.call(rbind, tables)
##count no. of each viral family
viral_families <- count(full_table, virus_family, species, sort=TRUE)

##reorder species
viral_families$species <- factor(viral_families$species, levels=c("M. aethiopoides Irish" ,"M. aethiopoides French", "M. aethiopoides Moroccan", "M. hyperodae"))

##plot
plot <- ggplot(viral_families, aes(x=species, y=n, fill=reorder(virus_family, n)))+
  scale_fill_viridis(discrete=TRUE, direction=-1)+
  geom_col()+
  theme_bw()+
  xlab("Species")+
  ylab("Number of peptides with viral hits")+
  labs(fill="Virus family")+
  guides(fill = guide_legend(reverse = TRUE))+
  coord_flip()+
  theme(axis.text.y = element_text(hjust = 1, face = "italic"))

##need to save as svg and edit in inkscape to un-italicise Irish, French, Moroccan
ggsave(file="output/virus_families.svg", plot=plot, width=8, height=3)
