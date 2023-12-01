################################################################################
################################################################################
#   Supp figure TAXONOMIC composition 
################################################################################
################################################################################

###################
## DATA AND CODE ##
###################


rm(list=ls())

# Load packages and functions
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(fantaxtic)

source("Scripts/V_submitted/Revision/0_Data_preparation/5_Load_Data.R")

# Get the most abundant phyla and the most abundant families within those phyla

top_nestedC <- nested_top_taxa(subset_samples(phyloseq_obj,suborder_host=="Caelifera"),
                              top_tax_level = "Order",
                              nested_tax_level = "Family",
                              n_top_taxa = 4,  n_nested_taxa = 3)

top_nestedC$top_taxa

# Plot the relative abundances at two levels.
Family_plotC = plot_nested_bar(ps_obj = top_nestedC$ps_obj, 
                              top_level = "Order", nested_level = "Family") + 
  MyTheme + ylab("Relative read counts") +
  facet_grid(~host_scientific_name, scales = "free_x", space= "free_x",drop=T) + 
  theme(legend.position="bottom",legend.key.size=unit(.45, "cm"),
        strip.text  = element_text(size=5,angle=90),
        strip.background=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ggtitle("Caelifera") 

Family_plotC

#Family_plot
ggsave(filename =  "Supp_Data_in_SM_doc/Supp_Fig_3_Taxonomy_GenusC.pdf",
       plot = Family_plotC, 
       width = 350, height = 100, device = 'pdf',units="mm")


top_nestedE <- nested_top_taxa(subset_samples(phyloseq_obj,suborder_host=="Ensifera"),
                               top_tax_level = "Order",
                               nested_tax_level = "Family",
                               n_top_taxa = 8,  n_nested_taxa = 6)

top_nestedE$top_taxa

# Plot the relative abundances at two levels.
Family_plotE = plot_nested_bar(ps_obj = top_nestedE$ps_obj, 
                               top_level = "Order", nested_level = "Family") + 
  MyTheme + ylab("Relative read counts") +
  facet_grid(~host_scientific_name, scales = "free_x", space= "free_x",drop=T) + 
  theme(legend.position="bottom",legend.key.size=unit(.45, "cm"),
        strip.text  = element_text(size=5,angle=90),
        strip.background=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(legend.position="bottom",legend.key.size=unit(.45, "cm")) +
  ggtitle("Ensifera")

#Family_plot
ggsave(filename =  "Supp_Data_in_SM_doc/Supp_Fig_3_Taxonomy_GenusE.pdf",
       plot = Family_plotE, 
       width = 350, height = 100, device = 'pdf',units="mm")

