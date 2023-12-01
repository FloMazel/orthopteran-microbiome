#######################################
#  #
#######################################

# Packages and code 
library(tidyverse)
library(phyloseq)
library(RColorBrewer)
library(ape)

source("0_Data_preparation/0_Util_fonctions.R")

path_dada2 = "" # path for dada2 output or wherever you stored the output of the previous scripts
setwd(path_dada2)

# Load data 
OTU_table <- read.table("seqtab.nochim_selectedASVs_1000reads.txt")
dim(OTU_table) # 336 1957

Taxonomy <- read.table("ASV_taxonomy_filtered.txt")

metadata = read_csv("Submitted_metadata.csv") %>% 
  subset(sample_alias %in% rownames(OTU_table)) # 337 samples - remove samples with <1000 reads from metadata

phylogeny <- read.tree("ASV_filtered_aligned.tree")


# Format Phyloseq object
MT <- as.data.frame(metadata); rownames(MT)=MT$sample_alias
PT <- as.matrix(Taxonomy); rownames(PT)=PT[,'seq']

phyloseq_obj=phyloseq(tax_table(PT),otu_table(OTU_table,taxa_are_rows=F),sample_data(MT),phy_tree(phylogeny))
phyloseq_obj

saveRDS(phyloseq_obj,"Phyloseq_object.RDS")
