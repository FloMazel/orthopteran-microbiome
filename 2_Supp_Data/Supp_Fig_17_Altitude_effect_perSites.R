rm(list=ls())

# Load data
library(ggplot2)
library(ggpubr)
library(MicEco)
library(phyloseq)
library(microbiome)
library(magrittr)

source("Scripts/V_submitted/Public/0_Data_preparation/0_Util_fonctions.R")
source("Scripts/V_submitted/Public/0_Data_preparation/5_Load_Data.R")
#library(tidyterra)

#tally table elevation 
elevation_sites = Info_depth %>% 
  group_by(sample_site_ID) %>% sample_n(1) %>% 
  select(sample_site_ID, Elevation)


# Define key parameters 
focal_sp = c("Chorthippus_parallelus","Euthystira_brachyptera")
set.seed(19)
N_raref = 1000
#grad <- hypso.colors(10, "dem_poster")

# Rarefy
endoS_PS_rarefied <- endoS_PS %>% 
  rarefy_even_depth(N_raref) 
gutB_PS_rarefied <-  gutB_PS %>% 
  rarefy_even_depth(N_raref) 

# Define metric and corresponding output 
betaM="bray"
FigFileName = "Redaction/Submission/ISMEcom_revision/Supp_Data_in_SM_doc/Altitude_supp_fig_Bray.pdf"

betaM="Aitchison"
FigFileName = "Redaction/Submission/ISMEcom_revision/Supp_Data_in_SM_doc/Altitude_supp_fig_Aitchison.pdf"

export_plot = list()
# Species focused on 


#species_generalist = c("Euthystira_brachyptera")
#species_generalist = c("Chorthippus_parallelus")


for (species_generalist in focal_sp) {

### Subset sample list 

endoS_PS_rarefied_oneSpecies = endoS_PS_rarefied %>% 
  subset_samples(host_scientific_name %in% species_generalist)
meta_endo_rarefied = as_tibble(phyloseq::sample_data(endoS_PS_rarefied_oneSpecies )) 

gutB_PS_rarefied_oneSpecies = gutB_PS_rarefied %>% 
  subset_samples(host_scientific_name %in% species_generalist)
meta_gut_rarefied = as_tibble(phyloseq::sample_data(gutB_PS_rarefied_oneSpecies)) 


################################
### PLOT SAMPLING DESIGN / BALANCE 
################################

gut_tally = meta_gut_rarefied %>% 
  group_by(host_scientific_name, Elevation) %>% 
  summarise(n=n()) %>% 
  mutate(data="Gut symbiont")

endo_tally = meta_endo_rarefied %>% 
  group_by(host_scientific_name, Elevation) %>% 
  summarise(n=n()) %>% 
  mutate(data="Endosymbiont")

tally = rbind(endo_tally,gut_tally)

Fig1 = tally %>% ggplot(aes(y=n,x=factor(Elevation),fill=data)) + geom_bar(stat="identity",position="dodge") +
  MyTheme+scale_fill_manual(values = cols) + 
  xlab('Elevation')+ylab("# of samples") + 
  ggtitle(paste0("Sample distribution, ", str_replace(species_generalist, "_", " ") )) +
  theme(legend.position="bottom",legend.title = element_blank())


################################
#### Fig 1A. Endosymbiont   ####
################################

if (betaM=="Aitchison"){
  beta_bray_PhyloSeq_endo  <- endoS_PS_rarefied_oneSpecies %>% 
    microbiome::transform("clr") %>% 
    phyloseq::distance("euclidean")
  beta_bray_PhyloSeq_gutB  <- gutB_PS_rarefied_oneSpecies %>% 
    microbiome::transform("clr") %>% 
    phyloseq::distance("euclidean")
} else {
  beta_bray_PhyloSeq_endo  <- endoS_PS_rarefied_oneSpecies %>% 
    phyloseq::distance(betaM)
  beta_bray_PhyloSeq_gutB  <- gutB_PS_rarefied_oneSpecies %>% 
    phyloseq::distance(betaM)
}



########################################

# Calculate median beta between all possible pairs of sites 

#####################
## Endo symbionts ###
#####################

# data frame with all pairs of sites 
reshaped_Beta_endo = as.data.frame(melt(as.matrix(beta_bray_PhyloSeq_endo),varnames = c("SampleID_1","SampleID_2"))) %>% 
  subset(!SampleID_1==SampleID_2) %>% 
  mutate(Beta=value) %>% select(-value) %>% 
  left_join(select(Info_depth,sample_alias,sample_site_ID), by=c("SampleID_1"="sample_alias")) %>% 
  left_join(select(Info_depth,sample_alias,sample_site_ID), by=c("SampleID_2"="sample_alias")) %>% 
  group_by(sample_site_ID.x,sample_site_ID.y) %>% 
  summarise(Beta=median(Beta))

# dissimilarity object 
Median_beta_endo = select(reshaped_Beta_endo,sample_site_ID.x,sample_site_ID.y,Beta) %>% 
  spread(sample_site_ID.x, Beta, fill=NA) %>% 
  column_to_rownames(var="sample_site_ID.y") %>% 
  as.matrix


# Ordination
perSite_MDS_endo = Median_beta_endo %>% metaMDS() %>% extract2("points") %>% 
  as_tibble(rownames = "sample_site_ID")  %>% 
  left_join(elevation_sites)

# model endosymbionts
model_endo = adonis_OmegaSq(adonis2(Median_beta_endo ~ Elevation, data=perSite_MDS_endo))
model_endo 
Elevation_effect <- bquote(Omega["Elevation"]^2==.(round(model_endo$parOmegaSq[1],2))~" ; p-value="~.(model_endo$`Pr(>F)`[1]))


Fig2 = perSite_MDS_endo %>% 
  ggplot(aes(y=MDS1,x=MDS2,col=Elevation)) +
  geom_point(shape = 1,size = 2.2,colour = "black")+
  geom_point(alpha=1,size=1.7)  + MyTheme +
  ggtitle(paste0("Endosymbiont community, ", str_replace(species_generalist, "_", " ") )) +
  labs(subtitle=Elevation_effect) +
  xlab("Ordination axis 2") + ylab(("Ordination axis 1"))+
  scale_color_gradientn(colours = terrain.colors(7), na.value = NA) 
#leg <- get_legend(Fig1B)
#Legend = as_ggplot(leg)
#Fig1B = Fig1B + theme(legend.position="none")


#####################
## GUT SYMBIOMNTS ###
#####################

# data frame with all pairs of sites 
reshaped_Beta_gut = as.data.frame(melt(as.matrix(beta_bray_PhyloSeq_gutB),varnames = c("SampleID_1","SampleID_2"))) %>% 
  subset(!SampleID_1==SampleID_2) %>% 
  mutate(Beta=value) %>% select(-value) %>% 
  left_join(select(Info_depth,sample_alias,sample_site_ID), by=c("SampleID_1"="sample_alias")) %>% 
  left_join(select(Info_depth,sample_alias,sample_site_ID), by=c("SampleID_2"="sample_alias")) %>% 
  group_by(sample_site_ID.x,sample_site_ID.y) %>% 
  summarise(Beta=median(Beta))

# dissimilarity object 
Median_beta_gut = select(reshaped_Beta_gut,sample_site_ID.x,sample_site_ID.y,Beta) %>% 
  spread(sample_site_ID.x, Beta, fill=NA) %>% 
  column_to_rownames(var="sample_site_ID.y") %>% 
  as.matrix

# Ordination
perSite_MDS_gut = Median_beta_gut %>% metaMDS() %>% extract2("points") %>% 
  as_tibble(rownames = "sample_site_ID")  %>% 
  left_join(elevation_sites)

# model gut symbionts
model_gut =  adonis_OmegaSq(adonis2(Median_beta_gut ~ Elevation, data=perSite_MDS_gut))
model_gut
Elevation_effect <- bquote(Omega["Elevation"]^2==.(round(model_gut$parOmegaSq[1],2))~" ; p-value="~.(model_gut$`Pr(>F)`[1]))


Fig3 = perSite_MDS_gut %>% 
  ggplot(aes(y=MDS1,x=MDS2,col=Elevation)) +
  geom_point(shape = 1,size = 2.2,colour = "black")+
  geom_point(alpha=1,size=1.7)  + MyTheme +
  labs(title=paste0("Putative gut symbiont community, ", str_replace(species_generalist, "_", " ") ),
        subtitle=Elevation_effect) +
  xlab("Ordination axis 1") + ylab("Ordination axis 2")+
  scale_color_gradientn(colours = terrain.colors(7), na.value = NA)
Fig3

# Group and export plot
export_plot[[species_generalist]] = plot_grid(Fig1, Fig2, Fig3, ncol=1)
}


final_plot = plot_grid(export_plot[[1]], 
                       export_plot[[2]],
                       ncol=2)

ggsave(filename = FigFileName, 
       plot = final_plot, 
       width = 180, height = 210, device = 'pdf',units="mm")










