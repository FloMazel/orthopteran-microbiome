library(ape)
library(Biostrings)
library(seqinr)
library(tidyverse)
library(genbankr)
library(ggtree)

path_dada2 = "" # path for dada2 output

# This script is used to build a phylogeny for 1/all ASVs using fastTree and 2/Wolbachia ASVs suign a mroe refined algorithm 

##################################
##################################
#     1/ ALL ASVs PHYLOGENY        #
##################################
##################################

#align sequences
# to run on the terminal 
"cd /Users/fmazel/Data/Criquets/Processed_data/dada2_data_submission/"
"mafft ASV_filtered.fasta > ASV_filtered_aligned.fasta # v7.490"

aligned = read.fasta( paste0(path_dada2,"ASV_filtered_aligned.fasta")
write.fasta(sequences = as.list(aligned),
            names = export_taxonomy$seq, #eport taxonomy files from "2_ASVs_selection.R" script
            file.out = paste0(path_dada2,"ASV_filtered_aligned.fasta")

# Reconstruct Fast Tree V2.1.11 
# to run on the terminal 
"FastTree -gtr -nt ASV_filtered_aligned.fasta > ASV_filtered_aligned.tree" # to run on the terminal 


##################################
##################################
#     2/ WOLBACHIA PHYLOGENY        #
##################################
##################################
meta_backBone = read.table(file = "Centennial_review_Metadata.txt",sep = "\t",header = T) %>% 
  mutate(Seq_ID=NCBIAccession) 
# to download from Kaur R, Shropshire JD, Cross KL, Leigh B, Mansueto AJ, Stewart V, et al. Living in the endosymbiotic world of Wolbachia: A centennial review. Cell Host Microbe 2021; 29: 879–893. 


##########################
#     GENBANK  IDs       #
##########################

# Non circular Genome metadata (link seq ID to assembly ID)
path_to_master_assembly_project = "" # tath were genoées from above ref were downloaded
master_assembly_list = tibble(assembly_name = list.files(path_to_master_assembly_project)) %>% 
  separate(assembly_name,sep = "_",into=c(NA,"GCA","ASN",NA),remove = F) %>% 
  mutate(GCA_name = paste0("GCA_",GCA),
         type="non_circular") 

master_assembly_list$Seq_ID <- sapply(X = list.files(path_to_master_assembly_project,full.names = T),function(x){readGenBank(x)@version["accession.version"]})
master_assembly_list = master_assembly_list %>% select(GCA_name,Seq_ID,type)

# Circular genomes 
report_path="MyData/WolbachiaGenome/Centenial_reviewList/genome_assemblies_asm_struct/ncbi-genomes-2022-09-29/"
NCBI_metadata = tibble(report_name = list.files(report_path)) %>% 
  mutate(GCA_name = paste0("GCA_",unlist(lapply(report_name,function(x){strsplit(x,"_")[[1]][2]}))),
         type="circular") %>% 
  subset(!GCA_name%in%master_assembly_list$GCA_name)
NCBI_metadata$Seq_ID = sapply(NCBI_metadata$report_name,function(x){read.table(paste0(report_path,x))[[1]]})
NCBI_metadata = NCBI_metadata %>% select(GCA_name,Seq_ID,type)

# All downloaded genomes 
path_genome <- "MyData/WolbachiaGenome/Centenial_reviewList/ncbi-genomes-2022-09-23/Genomes"

download_list = tibble(fasta_file_name = list.files(path_genome,pattern = "gz.fa"),
                       fasta_file_path = list.files(path_genome,full.names = T,pattern = "gz.fa")) %>% 
  mutate(GCA_name = paste0("GCA_",unlist(lapply(fasta_file_name,function(x){strsplit(x,"_")[[1]][2]})))) %>% 
  left_join(rbind(master_assembly_list,NCBI_metadata))


##########################
# CONSTRUCT BACkBONE TREE #
##########################

# Put GCA name in Assembly fasta header
for (i in 1:dim(download_list)[1]){
  fa = read.fasta(download_list$fasta_file_path[i])
  write.fasta(sequences = fa ,
              names = paste(download_list$GCA_name[i],names(fa),sep="_"),
              file.out = paste0(download_list$fasta_file_path[i],".fa")) }

# Concatenate all genomes 
"cat Genomes/*.fna.gz.fa > all_Wolbachia.fna" # to run on the terminal 

# Extract rRNA from representative set of Wolbachia
"barrnap --outseq all_Wolbachia_rRNA.fasta all_Wolbachia.fna" # to run on the terminal 

# Retrieve 16S only 
rRNA <- "MyData/WolbachiaGenome/Centenial_reviewList/ncbi-genomes-2022-09-23/all_Wolbachia_rRNA.fasta"
p16S =  "MyData/WolbachiaGenome/Centenial_reviewList/all_Wolbachia_16S.fasta"

rRNA = read.fasta(rRNA)
names_extracted_rRNA <- names(rRNA)
to_keep <- grep("16S",names_extracted_rRNA)
seqs16S <- rRNA [to_keep]
seqs16S <- seqs16S [sapply(seqs16S,length)>1200] # keep only long sequences, <1200 pb

write.fasta(sequences = seqs16S,names = names(seqs16S),file.out = p16S)


# Align
# to run on the terminal 
"mafft all_Wolbachia_16S.fasta > all_Wolbachia_16S_aligned.fasta"

# Build tree (IQtree)
# to run on the terminal 
"iqtree -s all_Wolbachia_16S_aligned.fasta -m TEST -bb 1000 -alrt 1000 -nt AUTO"

# Check bakcbone trees
BackBone = read.tree("MyData/WolbachiaGenome/Centenial_reviewList/all_Wolbachia_16S_aligned.fasta.treefile")
plot(BackBone)

tips = tibble(raw_tips=BackBone$tip.label,
              GCA_name=paste0("GCA_",unlist(lapply(BackBone$tip.label,function(x){strsplit(x,"_")[[1]][5]}))),
              multiple16S=duplicated(GCA_name)) %>% 
  left_join(download_list) %>% 
  left_join(meta_backBone) %>% 
  subset(!is.na(Supergroup)) %>% 
  mutate(name_group=paste(Supergroup,GCA_name,sep="_")) 

# 16S copy number 
copies_16s <- tips %>% group_by(GCA_name) %>% summarise(n16s=n(),host=unique(Host))
summary(copies_16s$n16s)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.000   1.000   1.000   1.183   1.000   6.000 

tips = tips %>% 
  select(GCA_name,raw_tips,Supergroup) %>% 
  group_by(GCA_name) %>% sample_n(1) # keep only one 16S copy per genome 

BackBone = keep.tip(BackBone,tips$raw_tips)
BackBone$tip.label = tips$GCA_name[match(BackBone$tip.label,tips$raw_tips)]
plot(BackBone)

p = ggtree(BackBone)

p = p %<+% tips + 
  geom_tiplab(aes(fill = factor(Supergroup)),
              color = "black", # color for label font
              geom = "label",  # labels not text
              label.padding = unit(0.15, "lines"), # amount of padding around the labels
              label.size = 0)  # size of label border

##########################
# ADD V4 AMPLICON DATA   #
##########################

# add fragments to the alignment
Fasta_Wolbachia = "MyData/WolbachiaGenome/Centenial_reviewList/Wolbachia_sequences.fa"

Taxonomy <- read.table("ASV_taxonomy.txt")
Wolba =  Taxonomy %>%
  subset(Genus=="Wolbachia")

revComp = reverseComplement(DNAStringSet(Wolba$seq))
write.fasta(as.list(revComp),names=paste0("Wolba_ASV_",1:length(Wolba$seq)),file.out = Fasta_Wolbachia)
Wolba = Wolba %>% 
  mutate(namesPhylo=paste0("Wolba_ASV_",1:length(Wolba$seq)))



# align Amplicon to template 
# to run on the terminal 
"mafft --keeplength --addfragments Wolbachia_sequences.fa all_Wolbachia_16S_aligned.fasta > Wolbachia_sequences_aligned_CentenialR.fa"

# Build tree (IQtree) with constrain 
# to run on the terminal 
"iqtree -s Wolbachia_sequences_aligned_CentenialR.fa -m TEST -bb 1000 -alrt 1000 -nt AUTO -g all_Wolbachia_16S_aligned.fasta.treefile"

# collapse unsupported nodes
# to run on the terminal 
"iqtree -t Wolbachia_sequences_aligned_CentenialR.fa.treefile -minsupnew 90/95"

# Check Full Tree
FullTree = read.tree("MyData/WolbachiaGenome/Centenial_reviewList/Wolbachia_sequences_aligned_CentenialR.fa.treefile.collapsed")
FullTree = read.tree("MyData/WolbachiaGenome/Centenial_reviewList/Wolbachia_sequences_aligned_CentenialR.fa.treefile")

plot(FullTree)

tips_full = tibble(raw_tips=FullTree$tip.label) %>% 
  left_join(tips) %>% 
  mutate(Group=if_else(is.na(Supergroup),"ASVs",Supergroup))

(BackBone,tips$raw_tips)
= tips$GCA_name[match(BackBone$tip.label,tips$raw_tips)]

##############################
#      Supp. Figure 2
##############################
p = ggtree(FullTree)
p = p %<+% tips_full + 
  geom_tippoint(aes(color = factor(Group)))  # size of label border
p
ggsave(filename = "Supp_figures/Wolbchia_Full_tree.pdf",plot = p, device = "pdf",height = 15,width = 10,limitsize = FALSE)

# get only ASV tree
FullTreeNoC = read.tree("MyData/WolbachiaGenome/Centenial_reviewList/Wolbachia_sequences_aligned_CentenialR.fa.treefile")
ASV_tree = keep.tip(FullTreeNoC,Wolba$namesPhylo)
ASV_tree
write.tree(ASV_tree,file = "MyData/Microbiome_error1/Wolbachia_ASV.tre")
