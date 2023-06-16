
# Load data
set.seed(19)
N_min = 499
betaM="bray"
FigFileName = "Redaction/Submission/ISMEcom_revision/Supp_Data_in_SM_doc/Supp_Fig_8_Spiroplasma_wolbachia_Beta.pdf"

Spiro_samples = Info_depth %>% subset(depth_Spiro>N_min) %>% pull(sample_alias)
Spiro_PS <- Spiro_PS %>% subset_samples(sample_alias %in% Spiro_samples)

Wolba_samples = Info_depth %>% subset(depth_Wolba>N_min) %>% pull(sample_alias)
Wolba_PS <- Wolba_PS %>% subset_samples(sample_alias %in% Wolba_samples)


  ################################
  #### Fig 1A.    Spiro       ####
  ################################
  
  bray_Spiro_PS  <- Spiro_PS  %>% 
    distance(betaM)
  
  hist(c(bray_Spiro_PS ))
  #summary(c(bray_Spiro_PS_rarefied ))
  
  # Plot MDS
  MDS_bray_Spiro_PS <- bray_Spiro_PS %>% metaMDS() 
  
  MDS_data_bray_Spiro_PS <-  MDS_bray_Spiro_PS$points %>% as_tibble(rownames = "sample_alias") %>% 
    left_join(Info_depth)
  
  Fig1A = MDS_data_bray_Spiro_PS %>% 
    ggplot(aes(y=MDS1,x=MDS2,col=host_scientific_name,shape=host_scientific_name,fill=suborder_host)) +
    geom_point(alpha=.7,size=.7)  +#geom_convexhull(alpha=0.05,size=.2)+
    MyTheme +  scale_shape_manual(values = hostsepciers_display$shape) +
    scale_color_manual(values = hostsepciers_display$color) + #facet_wrap(vars(Suborder))+
    ggtitle("A. Spiroplasma") + theme(legend.position="none")
  Fig1A
  
  ################################
  #### Fig 1B. Gut symbiont   ####
  ################################
  sample_Wolba_unique=c("MEP1-A10","MEP3-D6","MEP1-D3","MEP4-G5")# to remove

  bray_Wolba_PS  <- prune_samples(!Wolba_PS@sam_data$sample_alias %in%  sample_Wolba_unique,Wolba_PS) %>% 
    distance(betaM)
  hist(c(bray_Wolba_PS ))
  
  # Plot MDS
  MDS_bray_Wolba_PS <- bray_Wolba_PS %>% metaMDS() 
  
  MDS_data_bray_Wolba_PS <-  MDS_bray_Wolba_PS$points %>% as_tibble(rownames = "sample_alias") %>% 
    left_join(Info_depth)
  
  hostsepciers_display = expand_grid(color=brewer.pal(n = 8, name="Set1"),shape=c(15,17,18))[1:23,] 
  
  Fig1B = MDS_data_bray_Wolba_PS %>% 
    ggplot(aes(y=MDS1,x=MDS2,col=host_scientific_name,shape=host_scientific_name,fill=suborder_host)) +#+ ylim(c(-.15,.25))+
    geom_point(alpha=.7,size=.7)  +#geom_convexhull(alpha=0.05,size=.2)+
    MyTheme +  scale_shape_manual(values = hostsepciers_display$shape) +
    scale_color_manual(values = hostsepciers_display$color) + #facet_wrap(vars(Suborder))+
    ggtitle("B. Wolbachia")  + theme(legend.position="none")
  
  Fig1B
  
  
  ################################
  #### Fig 1C. Model Result   ####
  ################################
  
  # model Spiro
  model_Spiro = adonis2(bray_Spiro_PS ~ host_scientific_name + host_sex + Elevation,
                       by="margin", data=MDS_data_bray_Spiro_PS )
  model_Spiro
  model_Spiro2 = adonis_OmegaSq(model_Spiro,partial = T)
  model_Spiro2
  
  # model gut symbionts
  model_Wolba = adonis2(bray_Wolba_PS ~ host_scientific_name + host_sex + Elevation ,
                      by="margin", data=MDS_data_bray_Wolba_PS)
  model_Wolba2 = adonis_OmegaSq(model_Wolba,partial = T)
  model_Wolba2
  
  # grouping results 
  
  Permanova_res = unnest(model_Wolba2) %>% 
    subset(!is.na(F)) %>% 
    mutate(Predictor=c("Species","Sex" ,"Elevation"),
           Symbiont = rep("Wolbachia",3)) %>% 
    mutate(Signif = ifelse(`Pr(>F)`>.05,"",ifelse(`Pr(>F)`>.01,"*",'**')))
  
  Permanova_res2 = unnest(model_Spiro2) %>% 
    subset(!is.na(F)) %>% 
    mutate(Predictor=c("Species","Sex" ,"Elevation"),
           Symbiont = rep("Spiroplasma",3)) %>% 
    mutate(Signif = ifelse(`Pr(>F)`>.05,"",ifelse(`Pr(>F)`>.01,"*",'**')))
  
  Permanova_res = rbind(Permanova_res,Permanova_res2)
  
  Permanova_res
  
  permanova_plot = Permanova_res %>% 
    subset(!is.na(F)) %>% 
    ggplot(aes(y=parOmegaSq,x=reorder(Predictor,desc(parOmegaSq)),fill=Symbiont))+
    geom_bar(stat='identity',position="dodge")+
    geom_text(aes(label=Signif), position=position_dodge(width=0.9), vjust=.5,size=6)+ ylim(c(0,.6))+
    MyTheme+#scale_fill_manual(values = cols)+#scale_color_manual(values = cols)+
    xlab('Predictor of microbiota beta-diversity')+ylab(expression(Omega^2*' (effect size) of perMANOVA'))+
    ggtitle("C")+
    theme(legend.position="bottom",legend.title = element_blank())
  
  permanova_plot
  
  
  ################################
  ####     FINAL FIGURE       ####
  ################################
  
  model_plot=plot_grid(Fig1A,Fig1B,permanova_plot,nrow = 3)
  ggsave(plot = model_plot ,filename = FigFileName,device="pdf",height = 195,width = 80,units="mm")
  