# Packages
library(cowplot)
library(ggplot2)
library(vegan)
library(tidyverse)
library(ape)
library(reshape)
library(RColorBrewer)
library(phyloseq)
library(microbiome)
library(magrittr)

source("0_Data_preparation/0_Util_fonctions.R")

# Colors 
endoS=c("Spiroplasma","Wolbachia")
cols=c('Endosymbiont'="#FFD92F","Putative gut symbiont"="#6ab41f")

phyloseq_obj = readRDS("Phyloseq_object.RDS")

sample_data(phyloseq_obj) <- sample_data(phyloseq_obj) %>% 
  as("data.frame")  %>%
  separate(col = host_scientific_name,sep = "_", into = c("Host_genus",NA),remove = F)  %>%
  mutate(Host_species=gsub("_"," ",host_scientific_name))

  
dna <- Biostrings::DNAStringSet(taxa_names(phyloseq_obj))
names(dna) <- taxa_names(phyloseq_obj)
phyloseq_obj <- merge_phyloseq(phyloseq_obj, dna)
taxa_names(phyloseq_obj) <- paste0("ASV", seq(ntaxa(phyloseq_obj)),"_",tax_tibble(phyloseq_obj)$Genus)
tax_table(phyloseq_obj) <- 
phyloseq_obj

taxo <- tax_tibble(phyloseq_obj)

# Create different Phyloseq objects 
#endoS_PS <- phyloseq_obj %>% 
#  subset_taxa(Genus%in%endoS)
#endoS_PS

# alternative for endoSymbionts: 
endoS_PS <- phyloseq_obj %>% 
  subset_taxa(Genus=="Spiroplasma"|Order=="Rickettsiales"|Order=="Chlamydiales") #|Order=="Diplorickettsiales")
endoS_PS

Spiro_PS <- phyloseq_obj %>% 
  subset_taxa(Genus=="Spiroplasma")
Spiro_PS

Wolba_PS <- phyloseq_obj %>% 
  subset_taxa(Genus=="Wolbachia")
Wolba_PS

# gutB_PS  <- phyloseq_obj %>% 
#  subset_taxa(!Genus%in%endoS)

gut_subset <- subset(t(otu_table(phyloseq_obj)), 
                     !colnames(otu_table(phyloseq_obj)) %in% taxa_names(endoS_PS))
gutB_PS  <- merge_phyloseq(t(gut_subset), 
                           tax_table(phyloseq_obj), 
                           sample_data(phyloseq_obj),
                           refseq(phyloseq_obj),
                           phy_tree(phyloseq_obj))

gutB_PS 
taxoGut <- tax_tibble(gutB_PS)
taxoGutTable <- taxoGut %>% 
  mutate(readCounts = taxa_sums(gutB_PS)) %>% 
  group_by(Order) %>% 
  summarise(n=n(),readCount=sum(readCounts))


# General Depth plots 
Info_depth = as(sample_data(endoS_PS), "data.frame") %>% 
  left_join(tibble(depth_endo=sample_sums(endoS_PS),sample_alias=names(sample_sums(endoS_PS)))) %>% 
  left_join(tibble(depth_gut=sample_sums(gutB_PS),sample_alias=names(sample_sums(gutB_PS))))  %>% 
  left_join(tibble(depth_Wolba=sample_sums(Wolba_PS),sample_alias=names(sample_sums(Wolba_PS)))) %>% 
  left_join(tibble(depth_Spiro=sample_sums(Spiro_PS),sample_alias=names(sample_sums(Spiro_PS))))

