library(cowplot)
library(vegan)
library(tidyverse)
library(phyloseq)

endoS=c("Spiroplasma","Wolbachia")
cols=c('Endosymbiont'="#FFD92F",'Gut symbiont'="#6ab41f")

source("Scripts/V_submitted/Public/0_Data_preparation/5_Load_Data.R")


# Richness estimates
Info_depth = Info_depth %>% 
  mutate(sample_alias2 = (str_replace(sample_alias,"-","."))) %>% 
  left_join(tibble(sample_alias2 = rownames(estimate_richness(Wolba_PS,measures = c("Chao1"))),Chao1_Wolba = estimate_richness(Wolba_PS,measures = c("Chao1"))$Chao1)) %>% 
  left_join(tibble(sample_alias2 = rownames(estimate_richness(Wolba_PS,measures = c("Observed"))),Observed_Wolba = estimate_richness(Wolba_PS,measures = c("Observed"))$Observed)) %>%
  
  left_join(tibble(sample_alias2 = rownames(estimate_richness(Spiro_PS,measures = c("Chao1"))),Chao1_Spiro = estimate_richness(Spiro_PS,measures = c("Chao1"))$Chao1)) %>% 
  left_join(tibble(sample_alias2 = rownames(estimate_richness(Spiro_PS,measures = c("Observed"))),Observed_Spiro = estimate_richness(Spiro_PS,measures = c("Observed"))$Observed)) %>%
  
  left_join(tibble(sample_alias2 = rownames(estimate_richness(endoS_PS,measures = c("Chao1"))),Chao1_endoS = estimate_richness(endoS_PS,measures = c("Chao1"))$Chao1)) %>% 
  left_join(tibble(sample_alias2 = rownames(estimate_richness(endoS_PS,measures = c("Observed"))),Observed_endoS = estimate_richness(endoS_PS,measures = c("Observed"))$Observed)) %>%
  
  left_join(tibble(sample_alias2 = rownames(estimate_richness(gutB_PS,measures = c("Chao1"))),Chao1_gutB = estimate_richness(gutB_PS,measures = c("Chao1"))$Chao1)) %>% 
  left_join(tibble(sample_alias2 = rownames(estimate_richness(gutB_PS,measures = c("Observed"))),Observed_gutB = estimate_richness(gutB_PS,measures = c("Observed"))$Observed))
  
# Plot Richness estimates

Gut_Alpha_Plot = Info_depth %>% 
  ggplot(aes(y = Chao1_gutB, x = host_scientific_name,color=suborder_host)) + geom_boxplot() + 
  MyTheme + facet_wrap(vars(suborder_host),scales="free_x") +ylab("Chao1 estimate of richness") + ggtitle("Putative gut Symbiont")

EndoS_Alpha_Plot = Info_depth %>% 
  ggplot(aes(y = Chao1_endoS, x = host_scientific_name,color=suborder_host)) + geom_boxplot() + 
  MyTheme + facet_wrap(vars(suborder_host),scales="free_x") +ylab("Chao1 estimate of richness") + ggtitle("Endosymbiont")

Wolba_Alpha_Plot = Info_depth %>% 
  ggplot(aes(y = Chao1_Wolba, x = host_scientific_name,color=suborder_host)) + geom_boxplot() + 
  MyTheme + facet_wrap(vars(suborder_host),scales="free_x") +ylab("Chao1 estimate of richness") + ggtitle("Wolbachia")

Spiro_Alpha_Plot = Info_depth %>% 
  ggplot(aes(y = Chao1_Spiro, x = host_scientific_name,color=suborder_host)) + geom_boxplot() + 
  MyTheme + facet_wrap(vars(suborder_host),scales="free_x") +ylab("Chao1 estimate of richness") + ggtitle("Spiroplasma")


alpha_plot = plot_grid(Gut_Alpha_Plot,
          Wolba_Alpha_Plot,
          Spiro_Alpha_Plot,
          align = "h",ncol=1)

ggsave("Redaction/Submission/ISMEcom_revision/Supp_Data_in_SM_doc/Suop_fig_5_alpha_plot.pdf",alpha_plot,device = 'pdf',
       width = 175,height =300  ,units="mm")

# Get some estimates of mean richness 

Info_depth %>% 
  subset(Chao1_endoS > 0) %>% 
  summarise(mean(Chao1_endoS),sd(Chao1_endoS))

Info_depth %>% 
  subset(Chao1_Spiro > 0) %>% 
  summarise(mean(Chao1_Spiro),sd(Chao1_Spiro))


Info_depth %>% 
  subset(Chao1_Wolba > 0) %>% 
  summarise(mean(Chao1_Wolba),sd(Chao1_Wolba))

Info_depth %>% 
  subset(Chao1_gutB > 0) %>% 
  summarise(mean(Chao1_gutB),sd(Chao1_gutB))

# Get some info on relative abundances 

Info_depth = Info_depth %>%
  left_join(tibble(maxAb_Wolba = apply(as.matrix(Wolba_PS@otu_table),1,max)/apply(as.matrix(Wolba_PS@otu_table),1,sum),sample_alias = rownames(Wolba_PS@otu_table))) %>% 
  left_join(tibble(maxAb_Spiro = apply(as.matrix(Spiro_PS@otu_table),1,max)/apply(as.matrix(Spiro_PS@otu_table),1,sum),sample_alias = rownames(Spiro_PS@otu_table)))

Info_depth %>% 
 subset(Chao1_Wolba > 0) %>% 
  summarise(mean(maxAb_Wolba),sd(maxAb_Wolba)) # 78% relative read count for most abundant Wolba ASVs when Wolba present 

Info_depth %>% 
  subset(Chao1_Spiro > 0) %>% 
  summarise(mean(maxAb_Spiro),sd(maxAb_Spiro)) # 99% relative read count for most abundant Wolba ASVs when Wolba present 


Wolba = Info_depth %>% 
  subset(Chao1_Wolba > 0) %>% 
  ggplot(aes(y=maxAb_Wolba,x="Wolbachia")) + geom_boxplot() + ylim(c(0,1)) +
  MyTheme +ylab("Relative read count of the most abundant endosymbiont ASV")

Spiro = Info_depth %>% 
  subset(Chao1_Spiro > 0) %>% 
  ggplot(aes(y=maxAb_Spiro,x="Spiroplasma")) + geom_boxplot() + ylim(c(0,1)) +
  MyTheme +ylab("Relative read count of the most abundant endosymbiont ASV")


domainant_plot = plot_grid(Wolba,
                           Spiro,
                       align = "h",ncol=2)

domainant_plot

ggsave("Redaction/Submission/ISMEcom_revision/Supp_Data_in_SM_doc/Supp_Fig6_dominantEndoS.pdf",domainant_plot,device = 'pdf',
       width = 175,height =80  ,units="mm")


