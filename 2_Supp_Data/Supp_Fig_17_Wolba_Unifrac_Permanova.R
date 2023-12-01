library(phyloseq)
library(vegan)
library(tidyverse)

Wolba_PS

wolba_phylo = read.tree("Wolbachia_ASV.tre")
wolba_phylo

Wolba =  Taxonomy %>%
  subset(Genus=="Wolbachia") %>% 
  mutate(namesPhylo=paste0("Wolba_ASV_",1:length(Wolba$seq)))


wolba_phylo$tip.label = Wolba$seq[match(wolba_phylo$tip.label,Wolba$namesPhylo)]

Wolba_PS2 = merge_phyloseq(Wolba_PS, wolba_phylo)
Wolba_PS2 

################################
#### Fig 1B. Gut symbiont   ####
################################
sample_Wolba_unique=c("MEP1-A10","MEP3-D6","MEP1-D3","MEP4-G5")# to remove

betaM="wunifrac"

bray_Wolba_PS  <- prune_samples(!Wolba_PS2@sam_data$SampleID %in%  sample_Wolba_unique,Wolba_PS2) %>% 
  phyloseq::distance(betaM)
hist(c(bray_Wolba_PS ))

# Plot MDS
MDS_bray_Wolba_PS <- bray_Wolba_PS %>% metaMDS() 

MDS_data_bray_Wolba_PS <-  MDS_bray_Wolba_PS$points %>% as_tibble(rownames = "SampleID") %>% 
  left_join(Info_depth)

hostsepciers_display = expand_grid(color=brewer.pal(n = 8, name="Set1"),shape=c(15,17,18))[1:23,] 

Fig1B = MDS_data_bray_Wolba_PS %>% 
  ggplot(aes(y=MDS1,x=MDS2,col=Species,shape=Species,fill=Suborder)) +#+ ylim(c(-.15,.25))+
  geom_point(alpha=.7,size=.7)  +#geom_convexhull(alpha=0.05,size=.2)+
  MyTheme +  scale_shape_manual(values = hostsepciers_display$shape) +
  scale_color_manual(values = hostsepciers_display$color) + #facet_wrap(vars(Suborder))+
  ggtitle("B. Wolbachia")  + theme(legend.position="none")

Fig1B

model_Wolba = adonis2(bray_Wolba_PS ~ Species + Sex + Elevation ,
                      by="margin", data=MDS_data_bray_Wolba_PS)


model_Wolba


#adonis2(formula = bray_Wolba_PS ~ Species + Sex + Elevation, data = MDS_data_bray_Wolba_PS, by = "margin")
#Df SumOfSqs      R2       F Pr(>F)    
#Species    17   7.3604 0.58697 16.9581  0.001 ***
#  Sex         1   0.0256 0.00204  1.0026  0.304    
#Elevation   1   0.0330 0.00263  1.2934  0.226    
#Residual  198   5.0552 0.40314                   
#Total     217  12.5396 1.00000     
