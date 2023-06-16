############################################
############################################
# Figure 2 on host species specificity  
############################################
############################################

# Summary 
# Main figure
# Supp figures 


###################
## DATA AND CODE ##
###################


rm(list=ls())

# Load data
library(ggplot2)
library(ggpubr)
library(MicEco)
library(phyloseq)
library(microbiome)
library(magrittr)

source("Scripts/V_submitted/Public/0_Data_preparation/0_Util_fonctions.R")
source("Scripts/V_submitted/Public/1_Main_figures/0_Functions_figures.R")
source("Scripts/V_submitted/Public/0_Data_preparation/5_Load_Data.R")

# Rarefy data
set.seed(19)
N_raref = 1000

endoS_PS_rarefied <- endoS_PS %>% 
  rarefy_even_depth(N_raref)
endoS_PS_rarefied #259 samples

gutB_PS_rarefied <-  gutB_PS %>% 
  rarefy_even_depth(N_raref)
gutB_PS_rarefied # 167 samples


gutSamples = sample_data(gutB_PS_rarefied)@row.names
endoSamples = sample_data(endoS_PS_rarefied)@row.names

length(unique(c(gutSamples,endoSamples))) # 334 
length(intersect(gutSamples,endoSamples)) #92

92/259
92/167


PhyloSeq_endo=endoS_PS_rarefied
PhyloSeq_gut=gutB_PS_rarefied


# Main figure (rarefied data + bray crutis)
betaM="bray"
FigFileName = "Redaction/Submission/ISMEcom_revision/Main_Figure/Figure_2.pdf"
Specificity_figure_stats(gutB_PS_rarefied,endoS_PS_rarefied,betaM,FigFileName)

# Supp Mat Jaccard (rarefied data + jaccard)
betaM="jaccard"
FigFileName = "Redaction/Submission/ISMEcom_revision/Supp_Data/Supp_Fig_10_JaccardRaref.pdf"
Specificity_figure_stats(gutB_PS_rarefied,endoS_PS_rarefied,betaM,FigFileName)

# Supp Mat Aitchison
betaM="Aitchison"
FigFileName = "Redaction/Submission/ISMEcom_revision/Supp_Data/Supp_Fig_XX_AitchisonRaref.pdf"
Specificity_figure_stats(gutB_PS_rarefied,endoS_PS_rarefied,betaM,FigFileName)

# Supp Mat Unifrac
#betaM="wunifrac"
#FigFileName = "Redaction/Submission/ISMEcom_revision/Supp_Data/Supp_Fig_XX_W_UnifracRaref.pdf"
#Specificity_figure_stats(gutB_PS_rarefied,endoS_PS_rarefied,betaM,FigFileName)

# Supp Mat Unifrac
betaM="unifrac"
FigFileName = "Redaction/Submission/ISMEcom_revision/Supp_Data/Supp_Fig_XX_UnifracRaref.pdf"
Specificity_figure_stats(gutB_PS_rarefied,endoS_PS_rarefied,betaM,FigFileName)

# Supp Mat Non Raref
betaM="bray"
FigFileName = "Redaction/Submission/ISMEcom_revision/Supp_Data/Supp_Fig_9_BrayNon_raref.pdf"

endo_samples = Info_depth %>% subset(depth_endo>999) %>% pull(sample_alias)
endoS_PS_rarefied <- endoS_PS %>% subset_samples(sample_alias %in% endo_samples)

gut_samples = Info_depth %>% subset(depth_gut>999) %>% pull(sample_alias)
gutB_PS_rarefied <-  gutB_PS %>% subset_samples(sample_alias %in% gut_samples)

Specificity_figure_stats(gutB_PS_rarefied,endoS_PS_rarefied,betaM,FigFileName)
