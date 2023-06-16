# Packages
library(cowplot)
library(ggplot2)
library(vegan)
library(tidyverse)
library(ape)
library(reshape)
library(RColorBrewer)
library(phyloseq)
source("Scripts/V_submitted/Public/0_Data_preparation/0_Util_fonctions.R")

# Colors 
endoS=c("Spiroplasma","Wolbachia")
cols=c('Endosymbiont'="#FFD92F",'Gut symbiont'="#6ab41f")
hostsepciers_display = expand_grid(color=brewer.pal(n = 8, name="Set1"),shape=c(15,17,18))[1:22,] 

phyloseq_obj = readRDS("MyData/1_Important_Not_submitted_data/Phyloseq_object.RDS")

# Create different Phyloseq objects 
endoS_PS <- phyloseq_obj %>% 
  subset_taxa(Genus%in%endoS)
endoS_PS

Spiro_PS <- phyloseq_obj %>% 
  subset_taxa(Genus=="Spiroplasma")
Spiro_PS

Wolba_PS <- phyloseq_obj %>% 
  subset_taxa(Genus=="Wolbachia")
Wolba_PS

gutB_PS  <- phyloseq_obj %>% 
  subset_taxa(!Genus%in%endoS)

gutB_PS 

# General Depth plots 
Info_depth = as(sample_data(endoS_PS), "data.frame") %>% 
  left_join(tibble(depth_endo=sample_sums(endoS_PS),sample_alias=names(sample_sums(endoS_PS)))) %>% 
  left_join(tibble(depth_gut=sample_sums(gutB_PS),sample_alias=names(sample_sums(gutB_PS))))  %>% 
  left_join(tibble(depth_Wolba=sample_sums(Wolba_PS),sample_alias=names(sample_sums(Wolba_PS)))) %>% 
  left_join(tibble(depth_Spiro=sample_sums(Spiro_PS),sample_alias=names(sample_sums(Spiro_PS))))

