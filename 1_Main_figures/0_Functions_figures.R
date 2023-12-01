# Function to perform compositional analysis of 16s metabarcoding data, i.e. compute beta-diversity, produce ordination plots and PERMANOVA tests. 
# This helps to produce figure 2 of the manuscript, as well as associated supp figures. 

# input are 
# phyloseq object for putative gut symbionts end endosymbionts
# beta-diveristy metric (unifrac, jaccard, bray curtis or aitchison)
# path to export correspoinding figure in pdf format 

Specificity_figure_stats=function(PhyloSeq_gut,PhyloSeq_endo,betaM,FigFileName){
  
  species_list = unique(Info_depth$host_scientific_name)
  
  aes_code_host = tibble(color=colsHost$Host,
                         shape=rep(c(15,17,18),8),
                         host_scientific_name = names(colsHost$Host))
  
  color_code = aes_code_host$color; names(color_code)=aes_code_host$host_scientific_name
  shape_code = aes_code_host$shape; names(shape_code)=aes_code_host$host_scientific_name
  
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

# Plot 
  Fig1A = MDS_data_bray_PhyloSeq_endo %>% 
    ggplot(aes(y=MDS1,x=MDS2,col=Host_species,shape=Host_species,fill=suborder_host)) + #ylim(-.3,.3)+
    geom_point(alpha=.7,size=.7)  +#geom_convexhull(alpha=0.05,size=.2)+
    MyTheme +  scale_shape_manual(values = shape_code) +
    scale_color_manual(values = color_code) + #facet_wrap(vars(suborder_host))+
    theme(legend.position="none") + xlab("Ordination axis 2") + ylab(("Ordination axis 1"))+
    labs(title = "A. Ordination <span style='color: #FFD92F;'>(endosymbionts)<span>") + 
    theme(plot.title = element_markdown(hjust = 0,face="bold"))
  #theme(plot.margin = unit(c(5.5, 5.5,0, 5.5), "pt")) 
  
  print(Fig1A)

# Extract beta-div values 
flat_beta_endo <- reshape2::melt(as.matrix(beta_bray_PhyloSeq_endo), 
                     varnames=c("Sample_alias_1","Sample_alias_2"),
                     value.name = "beta_diversity")
flat_beta_endo <- subset(flat_beta_endo, !Sample_alias_1==Sample_alias_2) # remove comparison of a sample with itself (beta=0)
  
# Add metadata info to the beta-div flat table 
flat_beta_endo <- left_join(x = flat_beta_endo, 
                         y = MDS_data_bray_PhyloSeq_endo, 
                         by = c("Sample_alias_1"="sample_alias")) # add metadata for sample

flat_beta_endo <- left_join(x = flat_beta_endo, 
                         y = MDS_data_bray_PhyloSeq_endo, 
                         by = c("Sample_alias_2"="sample_alias"),
                         suffix = c("_1","_2")) # add metadata for sample
  
# Create within versus between groups categories 
flat_beta_endo <- mutate(flat_beta_endo,
                      comparisons = ifelse((!host_scientific_name_1==host_scientific_name_2),
                                           "Different \nhost sp.","Same \nhost sp."))

  
comparison_beta_plot_endo <- ggplot(flat_beta_endo, aes(x=reorder(comparisons,desc(beta_diversity))
                                                ,y=beta_diversity)) +
    MyTheme + 
    geom_boxplot(alpha=.005,outlier.size = .2) + 
    theme(axis.text.x = element_text(angle=70),
          axis.title.x = element_blank()) +
  #ggtitle("B") +
    ylab("Dissimilarity") +#+ theme(plot.margin = unit(c(5.5, 5.5, 0, 5.5), "pt")) 
labs(title = "B. Dissimilarity <span style='color: #FFD92F;'>(endosymbionts)<span>") + 
  theme(plot.title = element_markdown(hjust = 0,face="bold"))

  
comparison_beta_plot_endo
    


  
  ###############################
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
  
  
  # Plot MDS
  Fig1B = MDS_data_bray_gut %>% ggplot(aes(y=MDS1,x=MDS2,col=Host_species,shape=Host_species)) +
    geom_point(alpha=.7,size=.7) + 
    scale_shape_manual(values = shape_code)+
    scale_color_manual(values = color_code)+
    ggtitle("C")  +  MyTheme +
    theme(legend.spacing = unit(0, 'cm'),legend.key.size = unit(.3,"cm"),legend.text = element_text(size=4,face="italic")) +
    guides(shape = guide_legend(keywidth=0.2,keyheight=0.3,title="Host species",ncol=1),
           col = guide_legend(keywidth=0.2,keyheight=0.3,title="Host species"))+
    xlab("Ordination axis 2") + ylab(("Ordination axis 1"))+
    labs(title = "C. Ordination <span style='color: #6ab41f;'>(gut symbionts)<span>") + 
    theme(plot.title = element_markdown(hjust = 0,face="bold"))#+  theme(plot.margin = unit(c(5.5, 5.5,0, 5.5), "pt")) 

  print(Fig1B)
  leg <- get_legend(Fig1B)
  Legend = as_ggplot(leg)
  Fig1B = Fig1B + theme(legend.position="none")
  
  # Extract beta-div values 
  flat_beta_gut <- reshape2::melt(as.matrix(beta_bray_gut_PS), 
                                   varnames=c("Sample_alias_1","Sample_alias_2"),
                                   value.name = "beta_diversity")
  flat_beta_gut <- subset(flat_beta_gut, !Sample_alias_1==Sample_alias_2) # remove comparison of a sample with itself (beta=0)
  
  # Add metadata info to the beta-div flat table 
  flat_beta_gut <- left_join(x = flat_beta_gut, 
                              y = MDS_data_bray_gut, 
                              by = c("Sample_alias_1"="sample_alias")) # add metadata for sample
  
  flat_beta_gut <- left_join(x = flat_beta_gut, 
                              y = MDS_data_bray_gut, 
                              by = c("Sample_alias_2"="sample_alias"),
                              suffix = c("_1","_2")) # add metadata for sample
  
  # Create within versus between groups categories 
  flat_beta_gut <- mutate(flat_beta_gut,
                           comparisons = ifelse((!host_scientific_name_1==host_scientific_name_2),
                                                "Different \nhost sp.","Same \nhost sp."))
  
  
  comparison_beta_plot_gut <- ggplot(flat_beta_gut, aes(x=reorder(comparisons,desc(beta_diversity))
                                                          ,y=beta_diversity)) +
    MyTheme + 
    geom_boxplot(alpha=.005,outlier.size = .2) + 
    theme(axis.text.x = element_text(angle=70),
          axis.title.x = element_blank()) +
    ggtitle("D")+
    ylab("Dissimilarity") +#+ theme(plot.margin = unit(c(5.5, 5.5,0, 5.5), "pt")) 
    labs(title = "D. Dissimilarity <span style='color: #6ab41f;'>(gut symbionts)<span>") + 
    theme(plot.title = element_markdown(hjust = 0,face="bold"))#+  theme(plot.margin = unit(c(5.5, 5.5,0, 5.5), "pt")) 
  
  comparison_beta_plot_gut
  
  
  ################################
  #### Fig 1C. Model Result   ####
  ################################
  
  # model endosymbionts
  model_endo = adonis2(beta_bray_PhyloSeq_endo ~ host_scientific_name + host_sex + Elevation,
                       by="margin", data=MDS_data_bray_PhyloSeq_endo )
  #model_endo = adonis2(beta_bray_PhyloSeq_endo ~ host_scientific_name + host_sex + sample_site_ID,
  #                     by="margin", data=MDS_data_bray_PhyloSeq_endo )
  #model_endo
  dim(MDS_data_bray_PhyloSeq_endo)
  model_endo2 = adonis_OmegaSq(model_endo,partial = T)
  model_endo2$parOmegaSq
  
  # model gut symbionts
  model_gut = adonis2(beta_bray_gut_PS ~ host_scientific_name + host_sex + Elevation ,
                      by="margin", data=MDS_data_bray_gut)
  #model_gut = adonis2(beta_bray_gut_PS ~ host_scientific_name + host_sex + sample_site_ID ,
  #                    by="margin", data=MDS_data_bray_gut)
  dim(MDS_data_bray_gut)
  model_gut 
  model_gut2 = adonis_OmegaSq(model_gut)
  model_gut2
  
  # grouping results 
  
  Permanova_res = unnest(model_gut2) %>% 
    subset(!is.na(F)) %>% 
    mutate(Predictor=c("Host species","Sex" ,"Elevation"),
    #mutate(Predictor=c("Host species","Sex" ,"Site"),
           Symbiont = rep("Putative gut symbiont",3)) %>% 
    mutate(Signif = ifelse(`Pr(>F)`>.05,"",ifelse(`Pr(>F)`>.01,"*",'**')))
  
  Permanova_res2 = unnest(model_endo2) %>% 
    subset(!is.na(F)) %>% 
    mutate(Predictor=c("Host species","Sex" ,"Elevation"),
    #mutate(Predictor=c("Host species","Sex" ,"Site"),
           Symbiont = rep("Endosymbiont",3)) %>% 
    mutate(Signif = ifelse(`Pr(>F)`>.05,"",ifelse(`Pr(>F)`>.01,"*",'**')))
  
  Permanova_res = rbind(Permanova_res,Permanova_res2)
  
  Permanova_res
  
  permanova_plot = Permanova_res %>% 
    subset(!is.na(F)) %>% 
    ggplot(aes(y=parOmegaSq,x=reorder(Predictor,dplyr::desc(parOmegaSq)),fill=Symbiont))+
    geom_bar(stat='identity',position="dodge")+
    geom_text(aes(label=Signif), position=position_dodge(width=0.9), vjust=.5,size=3)+ ylim(c(0,0.7))+
    MyTheme+scale_fill_manual(values = cols, name="Microbiota")+#scale_color_manual(values = cols)+
    xlab('Predictor of dissimilarities')+ylab(expression(Omega^2*' (effect size) of perMANOVA'))+
    ggtitle("E. Statistics")+
    #theme(legend.position="bottom",legend.title = element_blank())
    theme(legend.position="right", legend.spacing = unit(0, 'cm'),legend.key.size = unit(.3,"cm"),legend.text = element_text(size=4)) +
    guides(col = guide_legend(keywidth=0.2,keyheight=0.5))
    
  permanova_plot
  
  print("Ratio of omega 2 of host species:  endo/gut symbiont ")
  print(betaM)
  print(Permanova_res %>% subset(Symbiont == "Endosymbiont" & Predictor == "Host species") %>% pull(parOmegaSq) / Permanova_res %>% subset(Symbiont == "Putative gut symbiont" & Predictor == "Host species") %>% pull(parOmegaSq))
  
  print("permanova results")
  print(Permanova_res)
  
  print("model gut ")
  print(model_gut2)
        
  print("model endo")
  print( model_endo2)
  
  
  ################################
  ####     FINAL FIGURE       ####
  ################################
  
# Endo plots 
Fig1AA=plot_grid(Fig1A,comparison_beta_plot_endo,nrow = 1,rel_widths = c(.8,.6), align = "h", axis = 'b')
Fig1AAA=plot_grid(Fig1AA, ggdraw(),nrow = 2,rel_heights = c(.9,.1)) # add padding

# gut plots 
Fig1BB=plot_grid(Fig1B,comparison_beta_plot_gut,nrow = 1,rel_widths = c(.8,.6), align = "h", axis = 'b')
Fig1BBB=plot_grid(Fig1BB, ggdraw(), nrow = 2,rel_heights = c(.9,.1)) # add padding

# Endo+Gut plots 
Fig1AB=plot_grid(Fig1AAA,Fig1BBB,nrow = 2)

# Legends
leg <- get_legend(permanova_plot)
Legend_perma = as_ggplot(leg) + theme(plot.margin = ggplot2::margin(-30, 0, 0, 0)) 
Legend_beta <- Legend + theme(plot.margin = ggplot2::margin(0, 0, 0, 0)) 
title_Legends <- ggdraw() + draw_label("Legends", fontface = 'bold', size = P3) 
Legends <- plot_grid(title_Legends, Legend_beta, Legend_perma,nrow = 3, rel_heights = c(.1,.6,.3))

# Permanova plots 
permanova_plot_noLegend= permanova_plot + theme(legend.position="none")

# Legend and permanova plots
right_plots = plot_grid(Legends, permanova_plot_noLegend , nrow = 2, rel_heights = c(1,1))
right_plots = plot_grid(ggdraw(),right_plots ,nrow = 1, rel_widths = c(.3,.9))

# final plots 
model_plot = plot_grid(Fig1AB,right_plots,nrow = 1, rel_widths = c(.6,.4), align = "h")
ggsave(plot = model_plot ,filename = FigFileName,device="pdf",height = 140,width = 120,units="mm")

}