# Contraint PERMANOVA per site 

rm(list=ls())

# Load data
library(ggplot2)
library(ggpubr)
library(MicEco)
library(phyloseq)
library(microbiome)
library(magrittr)

source("0_Data_preparation/0_Util_fonctions.R")
source("1_Main_figures/0_Functions_figures.R")
source("0_Data_preparation/5_Load_Data.R")

# Rarefy data
set.seed(19)
N_raref = 1000

endoS_PS
endoS_PS_rarefied <- endoS_PS %>% 
  rarefy_even_depth(N_raref)
endoS_PS_rarefied #259 samples

gutB_PS_rarefied <-  gutB_PS %>% 
  rarefy_even_depth(N_raref)
gutB_PS_rarefied # 167 samples

PhyloSeq_endo=endoS_PS_rarefied
PhyloSeq_gut=gutB_PS_rarefied

permanova_plot = list()

# Main figure 
FigFileName = "Supp_Data_in_SM_doc/Constranied_permanova_plot.pdf"

for (betaM in c("bray","jaccard","Aitchison","unifrac")) {
################################
#### Fig 1A. Endosymbiont   ####
################################
  
  # Prune samples  
  sample_endo_unique=c("MEP2-E7","MEP3-D6","MEP4-G5","MEP1-D7") # to remove
  bray_PhyloSeq_endo  <- prune_samples(!PhyloSeq_endo@sam_data$sample_alias %in% sample_endo_unique,PhyloSeq_endo)
  
  # Beta 
  
  if (betaM=="jaccard"){
    beta_bray_PhyloSeq_endo  <- bray_PhyloSeq_endo %>% 
      phyloseq::distance(betaM, binary = TRUE)
  } else if (betaM=="Aitchison"){
    beta_bray_PhyloSeq_endo  <- bray_PhyloSeq_endo %>% 
      microbiome::transform("clr") %>% 
      phyloseq::distance("euclidean")
  } else {  beta_bray_PhyloSeq_endo  <- bray_PhyloSeq_endo %>% phyloseq::distance(betaM)}
  
  hist(c(beta_bray_PhyloSeq_endo ))
  print(betaM)
  
  # Ordination 
  if (betaM=="Aitchison") { 
    MDS_bray_PhyloSeq_endo <- beta_bray_PhyloSeq_endo %>% pcoa() %>% 
      extract2("vectors") %>% magrittr::extract(,1:2) %>% as_tibble(rownames = "sample_alias") %>% 
      dplyr::rename(MDS1=Axis.1,MDS2=Axis.2)
  } else { 
    MDS_bray_PhyloSeq_endo <- beta_bray_PhyloSeq_endo %>% metaMDS() %>% extract2("points") %>% 
      as_tibble(rownames = "sample_alias")}
  
  MDS_data_bray_PhyloSeq_endo <-  MDS_bray_PhyloSeq_endo  %>% 
    left_join(Info_depth)
  
  head(MDS_data_bray_PhyloSeq_endo)
  
  
  ################################
  #### Fig 1B. Gut symbiont   ####
  ################################
  
  # Prune samples
  sample_gut_unique=c("MEP2-B7","MEP1-G7","MEP2-A4","MEP4-D12")# to remove
  bray_gut_PS  <- prune_samples(!PhyloSeq_gut@sam_data$sample_alias %in% sample_gut_unique,PhyloSeq_gut)
  
  
  if (betaM=="Aitchison"){
    beta_bray_gut_PS  <- bray_gut_PS  %>% 
      microbiome::transform("clr") %>% 
      phyloseq::distance("euclidean")
  } else if (betaM=="jaccard"){
    beta_bray_gut_PS  <- bray_gut_PS %>% 
      phyloseq::distance(betaM, binary = TRUE)
  } else {beta_bray_gut_PS = bray_gut_PS %>% distance(betaM)}
  
  
  # Ordination 
  if (betaM=="Aitchison") { 
    MDS_bray_gut_PS <- beta_bray_gut_PS %>% pcoa() %>% 
      extract2("vectors") %>% magrittr::extract(,1:2) %>% as_tibble(rownames = "sample_alias") %>% 
      dplyr::rename(MDS1=Axis.1,MDS2=Axis.2)
  } else { 
    MDS_bray_gut_PS <- beta_bray_gut_PS %>% metaMDS() %>% extract2("points") %>% 
      as_tibble(rownames = "sample_alias")}
  
  MDS_data_bray_gut <-  MDS_bray_gut_PS %>% 
    left_join(Info_depth)
  
  
  ################################
  #### Fig 1C. Model Result   ####
  ################################
  
  # model endosymbionts CONSTRAINED 
  model_endoCONSTRAINED = with(MDS_data_bray_PhyloSeq_endo,adonis2(beta_bray_PhyloSeq_endo ~ host_scientific_name + host_sex,
                       by="margin", strata = sample_site_ID ))
  model_endoCONSTRAINED
  model_endoCONSTRAINED2 = adonis_OmegaSq(model_endoCONSTRAINED,partial = T)
  model_endoCONSTRAINED2
  
  # model gut symbionts CONSTRAINED 
  model_gutCONSTRAINED  = with(MDS_data_bray_gut, adonis2(beta_bray_gut_PS ~ host_scientific_name + host_sex ,
                      by="margin", strata = sample_site_ID ))
  model_gutCONSTRAINED  
  model_gutCONSTRAINED2 = adonis_OmegaSq(model_gutCONSTRAINED )
  model_gutCONSTRAINED2
  
  # grouping results 
  
  Permanova_res = unnest(model_gutCONSTRAINED2) %>% 
    subset(!is.na(F)) %>% 
    mutate(Predictor=c("Host species","Sex"),
           Symbiont = rep("Putative gut symbiont",2)) %>% 
    mutate(Signif = ifelse(`Pr(>F)`>.05,"",ifelse(`Pr(>F)`>.01,"*",'**')))
  
  Permanova_res2 = unnest(model_endoCONSTRAINED2) %>% 
    subset(!is.na(F)) %>% 
    mutate(Predictor=c("Host species","Sex"),
           Symbiont = rep("Endosymbiont",2)) %>% 
    mutate(Signif = ifelse(`Pr(>F)`>.05,"",ifelse(`Pr(>F)`>.01,"*",'**')))
  
  Permanova_res = rbind(Permanova_res,Permanova_res2)
  
  Permanova_res
  
  permanova_plot[[betaM]] = Permanova_res %>% 
    subset(!is.na(F)) %>% 
    ggplot(aes(y=parOmegaSq,x=reorder(Predictor,dplyr::desc(parOmegaSq)),fill=Symbiont))+
    geom_bar(stat='identity',position="dodge")+
    geom_text(aes(label=Signif), position=position_dodge(width=0.9), vjust=.5,size=6)+ ylim(c(0,0.7))+
    MyTheme+scale_fill_manual(values = cols)+#scale_color_manual(values = cols)+
    xlab('Predictor of microbiota beta-diversity')+ylab(expression(Omega^2*' (effect size) of perMANOVA'))+
    ggtitle(betaM)+
    theme(legend.position="bottom",legend.title = element_blank())
  
  permanova_plot[[betaM]]
  
  }
  ################################
  ####     FINAL FIGURE       ####
  ################################
  

  Fig_constrained_permanova = plot_grid(permanova_plot[["bray"]],
                                        permanova_plot[["jaccard"]],
                                        permanova_plot[["Aitchison"]],
                                        permanova_plot[["unifrac"]],
                                        nrow = 2)
  #model_plot
  ggsave(plot = Fig_constrained_permanova ,filename = FigFileName,device="pdf",height = 195,width = 195,units="mm")
  
  
