rm(list=ls())

# Load data
library(ggplot2)
library(ggpubr)
library(MicEco)
library(phyloseq)
library(microbiome)
library(magrittr)
library(ape)

source("Scripts/V3/0_Data_preparation/0_Util_fonctions.R")
source("Scripts/V3/1_Main_figures/Functions_figures.R")
source("Scripts/V3/0_Data_preparation/5_Load_Data.R")
Host_Tree = read.tree("MyData/1_Important_Not_submitted_data/RAxML_bipartitions.result.31_tricked.newick")

# Rarefy data
set.seed(19)
N_raref = 1000

endoS_PS_rarefied <- endoS_PS %>% 
  rarefy_even_depth(N_raref)
endoS_PS_rarefied #259 samples

gutB_PS_rarefied <-  gutB_PS %>% 
  rarefy_even_depth(N_raref)
gutB_PS_rarefied # 167 samples

PhyloSeq_endo=endoS_PS_rarefied
PhyloSeq_gut=gutB_PS_rarefied

# beta 
betaM="bray"
exit_file="Redaction/Submission/ISMEcom_revision/Figure_3.pdf"
ylabT="Bray-curtis dissimilarity between host species"

set.seed(19)
N_raref = 1000
betaM="unifrac"
exit_file="Redaction/Submission/ISMEcom_revision/Supp_Data/SUpp_Fig_Mantel_Unifrac.pdf"
ylabT="Unifrac dissimilarity between host species"

set.seed(19)
N_raref = 1000
betaM="jaccard"
exit_file="Redaction/Submission/ISMEcom_revision/Supp_Data/SUpp_Fig_Mantel_Jaccard.pdf"
ylabT="Jaccard dissimilarity between host species"


############################################ 
###########    Endosymbiont.     ########### 
############################################

# Prune samples  
sample_endo_unique=c("MEP2-E7","MEP3-D6","MEP4-G5","MEP1-D7") # to remove
bray_PhyloSeq_endo  <- prune_samples(!PhyloSeq_endo@sam_data$sample_alias %in% sample_endo_unique,PhyloSeq_endo) %>% 
  phyloseq::distance(method = betaM)#, binary = T)

MDS_bray_endoS_PS_rarefied <- bray_PhyloSeq_endo  %>% metaMDS() 
MDS_data_bray_endoS_PS_rarefied <-  MDS_bray_endoS_PS_rarefied$points %>% as_tibble(rownames = "sample_alias") %>% 
  left_join(Info_depth)

# Get phylogenetic distance 
Pruned_Tree = keep.tip(phy = Host_Tree,tip = intersect(unique(MDS_data_bray_endoS_PS_rarefied$host_scientific_name),Host_Tree$tip.label))
Pruned_Tree # 21 host 
PhyDist=cophenetic(Pruned_Tree)
reshaped_D_phy=as.data.frame(melt(as.matrix(PhyDist),varnames = c("host_scientific_name.x","host_scientific_name.y"))) %>% 
  mutate(Phy_dist=value) %>% 
  select(-value)

# Reshape distance between pair of samples / species 
reshaped_Beta_endo = as.data.frame(melt(as.matrix(bray_PhyloSeq_endo),varnames = c("sample_alias_1","sample_alias_2"))) %>% 
  subset(!sample_alias_1==sample_alias_2) %>% 
  mutate(Beta=value) %>% select(-value) %>% 
  left_join(select(MDS_data_bray_endoS_PS_rarefied,sample_alias,host_scientific_name,suborder_host), by=c("sample_alias_1"="sample_alias")) %>% 
  left_join(select(MDS_data_bray_endoS_PS_rarefied,sample_alias,host_scientific_name,suborder_host), by=c("sample_alias_2"="sample_alias")) %>% 
  left_join(reshaped_D_phy) %>% 
  mutate(Comparisons = ifelse(host_scientific_name.x==host_scientific_name.y,"Same species",
                              ifelse(suborder_host.x==suborder_host.y,"Same Order different species","Different order")))

reshaped_Beta_endo_speciesMean = reshaped_Beta_endo %>% 
  group_by(host_scientific_name.x,host_scientific_name.y) %>% 
  summarise(Phy_dist=mean(Phy_dist),
            Beta=median(Beta),
          Comparisons=unique(Comparisons)) %>% 
  mutate(host_scientific_name.x=as.character(host_scientific_name.x),
         host_scientific_name.y=as.character(host_scientific_name.y))


MeanBetaEndo = select(reshaped_Beta_endo_speciesMean,host_scientific_name.x,host_scientific_name.y,Beta) %>% 
  spread(host_scientific_name.x, Beta, fill=NA) %>% 
  column_to_rownames(var="host_scientific_name.y") %>% 
  as.matrix

# Mantel test 
sp = colnames(MeanBetaEndo)
mantel_test_endo = mantel(MeanBetaEndo[sp,sp],PhyDist[sp,sp])
mantel_test_endo$statistic
mantel_test_endo$signif

Mantel_plot_endo = reshaped_Beta_endo_speciesMean  %>% 
  subset(!Comparisons=="Same species") %>% 
  ggplot(aes(y=Beta,x=Phy_dist)) + geom_point(alpha=.5,color="grey") +  MyTheme + ylim(.3,1)+
  geom_smooth(method="lm",color="black") + ggtitle("A. Endosymbiont community") +
  xlab("Phylogenetic distance between host species") + ylab(ylabT) + 
  annotate("text", x=.9, y=.6, size =P2/(72.27 / 25.4), label= paste0("Pearson r=",round(mantel_test_endo$statistic,3))) + 
  annotate("text", x=.9, y=.55, size =P2/(72.27 / 25.4),label= paste0("Mantel p=",mantel_test_endo$signif)) 

Mantel_plot_endo


############################################ 
###########    Gut symbionts    ########### 
############################################


# Prune samples  
sample_gut_unique=c("MEP2-B7","MEP1-G7","MEP2-A4","MEP4-D12")# to remove

bray_PhyloSeq_gut  <- prune_samples(!PhyloSeq_gut@sam_data$sample_alias %in% sample_endo_unique,PhyloSeq_gut) %>% 
  phyloseq::distance(method = betaM)#, binary=T)
MDS_bray_gut_PS_rarefied <- bray_PhyloSeq_gut  %>% metaMDS() 
MDS_data_bray_gut_PS_rarefied <-  MDS_bray_gut_PS_rarefied$points %>% as_tibble(rownames = "sample_alias") %>% 
  left_join(Info_depth)

# Get phylogenetic distance 
Pruned_Tree = keep.tip(phy = Host_Tree,tip = intersect(unique(MDS_data_bray_gut_PS_rarefied$host_scientific_name),Host_Tree$tip.label))
Pruned_Tree # 22 host 
PhyDist=cophenetic(Pruned_Tree)
reshaped_D_phy=as.data.frame(melt(as.matrix(PhyDist),varnames = c("host_scientific_name.x","host_scientific_name.y"))) %>% 
  mutate(Phy_dist=value) %>% 
  select(-value)

# Reshape distance between pair of samples / species 
reshaped_Beta_gut = as.data.frame(melt(as.matrix(bray_PhyloSeq_gut),varnames = c("sample_alias_1","sample_alias_2"))) %>% 
  subset(!sample_alias_1==sample_alias_2) %>% 
  mutate(Beta=value) %>% select(-value) %>% 
  left_join(select(MDS_data_bray_gut_PS_rarefied,sample_alias,host_scientific_name,suborder_host), by=c("sample_alias_1"="sample_alias")) %>% 
  left_join(select(MDS_data_bray_gut_PS_rarefied,sample_alias,host_scientific_name,suborder_host), by=c("sample_alias_2"="sample_alias")) %>% 
  left_join(reshaped_D_phy) %>% 
  mutate(Comparisons = ifelse(host_scientific_name.x==host_scientific_name.y,"Same species",
                              ifelse(suborder_host.x==suborder_host.y,"Same Order different species","Different order")))

reshaped_Beta_gut_speciesMean = reshaped_Beta_gut %>% 
  group_by(host_scientific_name.x,host_scientific_name.y) %>% 
  summarise(Phy_dist=mean(Phy_dist),
            Beta=median(Beta),
            Comparisons=unique(Comparisons)) %>% 
  mutate(host_scientific_name.x=as.character(host_scientific_name.x),
         host_scientific_name.y=as.character(host_scientific_name.y))


MeanBetagut = select(reshaped_Beta_gut_speciesMean,host_scientific_name.x,host_scientific_name.y,Beta) %>% 
  spread(host_scientific_name.x, Beta, fill=NA) %>% 
  column_to_rownames(var="host_scientific_name.y") %>% 
  as.matrix

# Mantel test 
sp = colnames(MeanBetagut)
mantel_test_gut = mantel(MeanBetagut[sp,sp],PhyDist[sp,sp])
mantel_test_gut$statistic
mantel_test_gut$signif

Mantel_plot_gut = reshaped_Beta_gut_speciesMean  %>% 
  subset(!Comparisons=="Same species") %>% 
  ggplot(aes(y=Beta,x=Phy_dist)) + geom_point(alpha=.5,color="grey") +  MyTheme + ylim(.3,1)+
  geom_smooth(method="lm",color="black") + ggtitle("B. Putative gut symbiont community") +
  xlab("Phylogenetic distance between host species") + ylab(ylabT) + 
  annotate("text", x=.9, y=.6, size =P2/(72.27 / 25.4), label= paste0("Pearson r=",round(mantel_test_gut$statistic,3))) + 
  annotate("text", x=.9, y=.55, size =P2/(72.27 / 25.4),label= paste0("Mantel p=",mantel_test_gut$signif)) 

Mantel_plot_gut

########################©########################©
########################©########################©
Mantel_plot_gut

model_plot=plot_grid(Mantel_plot_endo,Mantel_plot_gut,
                     nrow = 1)
model_plot
ggsave(plot = model_plot ,filename = exit_file,device="pdf",height = 80,width = 180,units="mm")

