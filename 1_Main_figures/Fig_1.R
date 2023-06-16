################################################################################
################################################################################
# Figure 1 on relative proportion of endosymbionts VS putative gut symbionts 
################################################################################
################################################################################

# Summary 
# Produce fig.1 
# Provide corresponding statistical tests 


###################
## DATA AND CODE ##
###################

rm(list=ls())
library(cowplot)
library(gridExtra)
library(tidyverse)
library(phyloseq)

source("Scripts/V_submitted/Public/0_Data_preparation/5_Load_Data.R")

ensifera_image <- svgparser::read_svg("Redaction/PhyloPic.8f55e2ec.Birgit-Lang.Ensifera.svg")
caelifera_image <- svgparser::read_svg("Redaction/PhyloPic.7c142ec5.Birgit-Lang.Caelifera.svg")
# Copyright 
# Caelifera http://phylopic.org/image/7c142ec5-aebb-495d-80fa-1b575090d5db/
# Ensifera http://phylopic.org/image/8f55e2ec-f2ea-407a-bfbf-077aa10b5d36/


###################
##    Fig 1A     ##
###################

# Define endosymbionts VS putative gut symbionts
endoS=c("Spiroplasma","Wolbachia")
cols=c('Endosymbiont'="#FFD92F",'Gut symbiont'="#6ab41f")

# Summarize data at the genus level 
Taxo_counts_Genus <- phyloseq_obj %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt()

Summary_Genus = Taxo_counts_Genus %>% 
  mutate(Endosymbiont = if_else(Genus%in%endoS,"Endosymbiont","Gut symbiont")) %>% 
  group_by(Endosymbiont) %>% 
  summarise(Sum_RelAb_bySample=sum(Abundance)/336)

Summary_Fam = Taxo_counts_Genus %>% 
  group_by(Family) %>% 
  summarise(Sum_RelAb_bySample=sum(Abundance)/336)

#


# Plot relative proportions 
TitleAdjust = -5
VadjustX =3 

Fig1A <- Summary_Genus %>%  
  ggplot(aes(y=Sum_RelAb_bySample,x=1, fill=Endosymbiont)) +  
  geom_bar(stat="identity",position = "stack") +
  scale_fill_manual(values=cols) +
  xlab("Total") + ylab("Relative read counts") + ggtitle("A. Total") +
  MyTheme + guides(fill=guide_legend(ncol=2)) +
  theme(axis.text.x.bottom =element_blank(), 
        axis.ticks.x=element_blank(),
        axis.line = element_blank(),
        axis.title.x = element_text(vjust =VadjustX ),
        plot.title = element_text(vjust = TitleAdjust)) +
  
  theme(legend.position="bottom",
        legend.title = element_blank(),
        legend.box = "horizontal",
        axis.text.x =element_text(size=P2,angle=0,hjust=0.5)) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1,size=3)))

legend <- cowplot::get_legend(Fig1A)

Fig1A  <- Fig1A  + theme(legend.position="none")
Fig1A

#ggsave(Fig1A,filename = "Redaction/V1/Main_Figure/Fig1A.pdf", ,device = 'pdf', width = 60,height =60  ,units="mm")

###################
##    Fig 1BC    ##
###################

Fig1Theme = theme(axis.title.x = element_text(vjust=VadjustX),
                  axis.title.y = element_blank(),
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank(),
                  axis.ticks = element_blank(),
                  axis.line = element_blank())

Summary_Endo_perSample = Taxo_counts_Genus %>% 
  mutate(Endosymbiont = if_else(Genus%in%endoS,"Endosymbiont","Gut symbiont")) %>% 
  group_by(Endosymbiont,Sample) %>% 
  summarise(Abundance=sum(Abundance),
  Suborder=unique(suborder_host)) 

Summary_Endo_perSample 

Fig1B <- Summary_Endo_perSample  %>% 
  subset(Suborder=="Ensifera") %>% 
  mutate(Abu_Endosymbiont=ifelse(Endosymbiont=='Endosymbiont', Abundance,0)) %>% 
  ggplot(aes(x = fct_reorder(Sample, Abu_Endosymbiont) , y = Abundance, fill=Endosymbiont)) + 
  geom_bar(position="fill",stat = "identity") +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  xlab("Samples") + 
  scale_fill_manual(values=cols) + ggtitle("B. Ensifera") +
  MyTheme + theme(legend.position="none",plot.title = element_text(vjust = TitleAdjust)) + 
  Fig1Theme +   annotation_custom(ensifera_image,xmin = 60,xmax = 2200,ymin=.8) 

Fig1B

Fig1C <- Summary_Endo_perSample  %>%
  subset(Suborder=="Caelifera") %>% 
  mutate(Abu_Endosymbiont=ifelse(Endosymbiont=='Endosymbiont', Abundance,0)) %>% 
  ggplot(aes(x = fct_reorder(Sample, Abu_Endosymbiont) , y = Abundance, fill=Endosymbiont)) + 
  geom_bar(position="fill",stat = "identity") +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative read counts") + xlab("Samples") + 
  scale_fill_manual(values=cols) + ggtitle("C. Caelifera") +
  MyTheme + theme(legend.position="none",plot.title = element_text(vjust = TitleAdjust))+
  Fig1Theme +   annotation_custom(caelifera_image,xmin = 200,xmax = 10000,ymin=.8,ymax=1) 

Fig1C

###################
##     FINAL     ##
###################

fig1top <- plot_grid(Fig1A,Fig1B,Fig1C,align = "h",
                  ncol=3,rel_widths=c(1,4,4))

fig1 = plot_grid(fig1top,legend,
                 ncol=1,rel_heights = c(7,1))

fig1

#ggsave("Redaction/Submission/ISMEcom_revision/Fig_1.pdf",fig1,device = 'pdf',
#       width = 175,height =80  ,units="mm")


###################
##    TEST       ##
###################

Summary_Categ_Sample = Taxo_counts_Genus %>% 
  mutate(Endosymbiont = if_else(Genus%in%endoS,"Endosymbiont","Gut symbiont")) %>% 
  group_by(Endosymbiont,Sample) %>% 
  summarise(RRC=sum(Abundance),Suborder=unique(suborder_host),Species=unique(host_scientific_name),Sex=unique(host_sex),Elevation=mean(Elevation)) %>% 
  pivot_wider(names_from = Endosymbiont,values_from=RRC) %>% 
  mutate(Ratio_EG = Endosymbiont/`Gut symbiont`,
         DominanceEndo=Endosymbiont>.5)


# Some basics stats on mean rel read counts 
Summary_Categ_Sample %>% 
  group_by(Suborder) %>% 
  summarise(mean(Endosymbiont),sd(Endosymbiont))

#  Suborder  `mean(Endosymbiont)` `sd(Endosymbiont)`
#1 Caelifera                0.814              0.271
#2 Ensifera                 0.155              0.277

# Some basics stats on prevalence of endoS


Summary_EndoS_Sample = Taxo_counts_Genus %>% 
  group_by(Genus,Sample) %>% 
  summarise(RRC=sum(Abundance),Suborder=unique(suborder_host),Species=unique(host_scientific_name),Sex=unique(host_sex),Elevation=mean(Elevation)) %>% 
  subset( Genus %in% endoS)

Summary_EndoS_Sample %>% 
  subset(RRC>.1) %>% 
  group_by(Genus) %>% 
  summarise(n_distinct(Species))
# Spiroplasma                    18 species with at least one indiviudal with >10# real. read counts 
# Wolbachia                      19 species with at least one indiviudal with >10# real. read counts 

# Some stats 
kruskal.test(Summary_Categ_Sample$Endosymbiont,Summary_Categ_Sample$Suborder)
#Kruskal-Wallis chi-squared = 136.41, df = 1, p-value < 2.2e-16

