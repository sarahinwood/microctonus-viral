library(data.table)
library(dplyr)
library(ggplot2)
library(viridis)

blast_res_files <- list.files(path="output/final_blast_tables/eukaryotic", pattern="*.csv",
                              full.names=TRUE, recursive=TRUE)
##read into single table for plotting
tables <- lapply(blast_res_files, fread, strip.white=FALSE)
full_table <- do.call(rbind, tables)
##count no. of each viral family
families <- count(full_table, family, species, sort=TRUE)

##reorder
families$species <- factor(families$species, levels=c("M. aethiopoides French", "M. aethiopoides Moroccan", "M. hyperodae"))

##plot
plot_euk <- ggplot(families, aes(x=species, y=n, fill=reorder(family, n)))+
  scale_fill_viridis(discrete=TRUE, direction=-1)+
  geom_col()+
  theme_bw()+
  xlab("Species")+
  ylab("Number of peptides with BlastP hits")+
  labs(fill="Eukaryotic families")+
  guides(fill = guide_legend(reverse = TRUE))+
  coord_flip()+
  theme(axis.text.y = element_text(hjust = 1, face = "italic"))+
  scale_y_continuous(breaks=c(0,2,4,6,8,10,12))

##need to save as svg and edit in inkscape to un-italicise French, Moroccan
ggsave(file="output/final_blast_tables/eukaryotic_families.svg", plot=plot_euk, width=8, height=4)
