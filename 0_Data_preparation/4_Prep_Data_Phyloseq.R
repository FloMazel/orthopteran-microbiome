#######################################
#  #
#######################################

# Packages and code 
library(tidyverse)
library(phyloseq)
library(RColorBrewer)
library(ape)

source("Scripts/V3/0_Data_preparation/0_Util_fonctions.R")

# Load data 
OTU_table <- read.table("MyData/1_Important_Not_submitted_data/seqtab.nochim_selectedASVs_1000reads.txt")
dim(OTU_table) # 336 1957

Taxonomy <- read.table("MyData/1_Important_Not_submitted_data/ASV_taxonomy_filtered.txt")

metadata = read_csv("MyData/2_Submission_data/Submitted_metadata.csv") %>% 
  subset(sample_alias %in% rownames(OTU_table)) # 337 samples - remove samples with <1000 reads from metadata

phylogeny <- read.tree("MyData/1_Important_Not_submitted_data/ASV_filtered_aligned.tree")


# Format Phyloseq object
MT <- as.data.frame(metadata); rownames(MT)=MT$sample_alias
PT <- as.matrix(Taxonomy); rownames(PT)=PT[,'seq']

phyloseq_obj=phyloseq(tax_table(PT),otu_table(OTU_table,taxa_are_rows=F),sample_data(MT),phy_tree(phylogeny))
phyloseq_obj

saveRDS(phyloseq_obj,"MyData/1_Important_Not_submitted_data/Phyloseq_object.RDS")
