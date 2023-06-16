library(ape)

source("Scripts/3_Load_Data.R")

# Host Tree
Host_Tree = read.tree("MyData/Host/RAxML_bipartitions.result.31_tricked.newick")
Pruned_Tree = keep.tip(phy = Host_Tree,tip = intersect(unique(Info_depth$Species),Host_Tree$tip.label))
Pruned_Tree # 21 host 
PhyDist=cophenetic(Pruned_Tree)
reshaped_D_phy=as.data.frame(melt(as.matrix(PhyDist),varnames = c("Species.x","Species.y"))) %>% 
  mutate(Phy_dist=value) %>% 
  select(-value)

# Wolba phylogeny 
Taxonomy <- read.table("MyData/Microbiome_error1/ASV_taxonomy.txt")
Wolba =  Taxonomy %>%
  subset(Genus=="Wolbachia")
Wolba = Wolba %>% 
  mutate(namesPhylo=paste0("Wolba_ASV_",1:length(Wolba$seq)))

wolba_phylo = read.tree("MyData/Microbiome_error1/Wolbachia_ASV.tre")
wolba_phylo

wolba_phylo$tip.label = Wolba$seq[match(wolba_phylo$tip.label,Wolba$namesPhylo)]

# Distributiob data 
N_min = 499
Wolba_samples = Info_depth %>% subset(depth_Wolba>N_min) %>% pull(SampleID)
Wolba_PS <- Wolba_PS %>% subset_samples(SampleID %in% Wolba_samples)
Wolba_PS2 = merge_phyloseq(Wolba_PS, wolba_phylo)
Wolba_PS2 
betaM="wunifrac"
#betaM="bray"
bray_Wolba_PS  <- Wolba_PS2 %>% 
  phyloseq::distance(betaM)

reshaped_Beta_endo = as.data.frame(melt(as.matrix(bray_Wolba_PS),varnames = c("SampleID_1","SampleID_2"))) %>% 
  subset(!SampleID_1==SampleID_2) %>% 
  mutate(Beta=value) %>% select(-value) %>% 
  left_join(select(Info_depth,SampleID,Species,Suborder), by=c("SampleID_1"="SampleID")) %>% 
  left_join(select(Info_depth,SampleID,Species,Suborder), by=c("SampleID_2"="SampleID")) %>% 
  left_join(reshaped_D_phy) %>% 
  mutate(Comparisons = ifelse(Species.x==Species.y,"Same species",
                              ifelse(Suborder.x==Suborder.y,"Same Order different species","Different order")))


reshaped_Beta_endo_speciesMean = reshaped_Beta_endo %>% 
  group_by(Species.x,Species.y) %>% 
  summarise(Phy_dist=mean(Phy_dist),
            Beta=median(Beta),
            Comparisons=unique(Comparisons)) %>% 
  mutate(Species.x=as.character(Species.x),
         Species.y=as.character(Species.y))


MeanBetaEndo = select(reshaped_Beta_endo_speciesMean,Species.x,Species.y,Beta) %>% 
  spread(Species.x, Beta, fill=NA) %>% 
  column_to_rownames(var="Species.y") %>% 
  as.matrix

# Mantel test 
sp = colnames(MeanBetaEndo)
mantel_test_endo = mantel(MeanBetaEndo[sp,sp],PhyDist[sp,sp])
mantel_test_endo$statistic
mantel_test_endo$signif

Mantel_plot_endo = reshaped_Beta_endo_speciesMean  %>% 
  subset(!Comparisons=="Same species") %>% 
  ggplot(aes(y=Beta,x=Phy_dist)) + geom_point(alpha=.5,color="grey") +  MyTheme +
  geom_smooth(method="lm",color="black") + ggtitle("Wolbachia") +
  xlab("Phylogenetic distance") + ylab("Weighted Unifrac Dissimilarity") + 
  annotate("text", x=.9, y=.6, size =P2/(72.27 / 25.4), label= paste0("Pearson r=",round(mantel_test_endo$statistic,3))) + 
  annotate("text", x=.9, y=.55, size =P2/(72.27 / 25.4),label= paste0("Mantel p=",mantel_test_endo$signif)) 

Mantel_plot_endo

ggsave(plot = Mantel_plot_endo ,filename = "Redaction/V1/Supp_figures/Mantel_wUnifravc_Wolbachia.pdf",device="pdf",height = 80,width = 80,units="mm")


