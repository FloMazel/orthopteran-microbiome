# Function to perform compositional analysis of 16s metabarcoding data, i.e. compute beta-diversity, produce ordination plots and PERMANOVA tests. 
# This helps to produce figure 2 of the manuscript, as well as associated supp figures. 

# input are 
# phyloseq object for putative gut symbionts end endosymbionts
# beta-diveristy metric (unifrac, jaccard, bray curtis or aitchison)
# path to export correspoinding figure in pdf format 

Specificity_figure_stats=function(PhyloSeq_gut,PhyloSeq_endo,betaM,FigFileName){
  
  species_list = unique(Info_depth$host_scientific_name)
  
  aes_code_host = expand_grid(color=brewer.pal(n = 8, name="Set1"),shape=c(15,17,18))[1:24,] %>% 
    mutate(host_scientific_name=species_list,
           host_scientific_name_ok=gsub("_", " ",host_scientific_name))
  
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
    ggplot(aes(y=MDS1,x=MDS2,col=host_scientific_name,shape=host_scientific_name,fill=suborder_host)) + #ylim(-.3,.3)+
    geom_point(alpha=.7,size=.7)  +#geom_convexhull(alpha=0.05,size=.2)+
    MyTheme +  scale_shape_manual(values = hostsepciers_display$shape) +
    scale_color_manual(values = hostsepciers_display$color) + #facet_wrap(vars(suborder_host))+
    ggtitle("A. Endosymbionts community") + theme(legend.position="none") + xlab("Ordination axis 2") + ylab(("Ordination axis 1"))
  
  print(Fig1A)
  
  
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
  
  
  # Plot MDS
  Fig1B = MDS_data_bray_gut %>% ggplot(aes(y=MDS1,x=MDS2,col=host_scientific_name,shape=host_scientific_name)) +
    geom_point(alpha=.7,size=.7) + 
    scale_shape_manual(values = shape_code,labels=aes_code_host$host_scientific_name_ok) + scale_color_manual(values = color_code,labels=aes_code_host$host_scientific_name_ok) + 
    ggtitle("B. Putative gut symbiont community")  +  MyTheme +
    theme(legend.spacing = unit(0, 'cm'),legend.key.size = unit(.3,"cm"),legend.text = element_text(size=4)) +
    guides(shape = guide_legend(keywidth=0.2,keyheight=0.5,title="Host species",ncol=1),
           col = guide_legend(keywidth=0.2,keyheight=0.5,title="Host species"))+
    xlab("Ordination axis 2") + ylab(("Ordination axis 1"))
  
  print(Fig1B)
  leg <- get_legend(Fig1B)
  Legend = as_ggplot(leg)
  Fig1B = Fig1B + theme(legend.position="none")
  
  ################################
  #### Fig 1C. Model Result   ####
  ################################
  
  # model endosymbionts
  model_endo = adonis2(beta_bray_PhyloSeq_endo ~ host_scientific_name + host_sex + Elevation,
                       by="margin", data=MDS_data_bray_PhyloSeq_endo )
  #model_endo = adonis2(beta_bray_PhyloSeq_endo ~ host_scientific_name + host_sex + sample_site_ID,
  #                     by="margin", data=MDS_data_bray_PhyloSeq_endo )
  #model_endo
  model_endo2 = adonis_OmegaSq(model_endo,partial = T)
  model_endo2$parOmegaSq
  
  # model gut symbionts
  model_gut = adonis2(beta_bray_gut_PS ~ host_scientific_name + host_sex + Elevation ,
                      by="margin", data=MDS_data_bray_gut)
  #model_gut = adonis2(beta_bray_gut_PS ~ host_scientific_name + host_sex + sample_site_ID ,
  #                    by="margin", data=MDS_data_bray_gut)
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
    geom_text(aes(label=Signif), position=position_dodge(width=0.9), vjust=.5,size=6)+ ylim(c(0,0.7))+
    MyTheme+scale_fill_manual(values = cols)+#scale_color_manual(values = cols)+
    xlab('Predictor of microbiota beta-diversity')+ylab(expression(Omega^2*' (effect size) of perMANOVA'))+
    ggtitle("C")+
    theme(legend.position="bottom",legend.title = element_blank())
  
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
  
  Fig1AB=plot_grid(Fig1A,Fig1B,nrow = 2)
  Fig1ABlegend = plot_grid(Fig1AB, Legend,ncol=2,rel_widths = c(.8,.4))
  model_plot = plot_grid(Fig1ABlegend,permanova_plot,nrow = 2,rel_heights = c(.6,.3))
  #model_plot
  ggsave(plot = model_plot ,filename = FigFileName,device="pdf",height = 195,width = 80,units="mm")
  
  
  
}

