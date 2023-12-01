
# Load data
source("0_Data_preparation/0_Util_fonctions.R")
source("/1_Main_figures/0_Functions_figures.R")
source("/0_Data_preparation/5_Load_Data.R")

library(gridExtra)


# Main figure 
set.seed(19)
N_raref = 1000
betaM="bray"
# Rarefy
endoS_PS_rarefied <- endoS_PS %>% 
  rarefy_even_depth(N_raref)

gutB_PS_rarefied <-  gutB_PS %>% 
  rarefy_even_depth(N_raref)


  ################################
  ####      Endosymbiont      ####
  ################################
  
  sample_endo_unique=c("MEP2-E7","MEP3-D6","MEP4-G5","MEP1-D7")# to remove
  bray_endoS_PS_rarefied  <- prune_samples(!endoS_PS_rarefied@sam_data$sample_alias %in% sample_endo_unique,endoS_PS_rarefied) %>% 
    distance(betaM)
  
  MDS_bray_endoS_PS_rarefied <- bray_endoS_PS_rarefied %>% metaMDS() 
  
  MDS_data_bray_endoS_PS_rarefied <-  MDS_bray_endoS_PS_rarefied$points %>% as_tibble(rownames = "sample_alias") %>% 
    left_join(Info_depth)
  
  model_endo_disp = betadisper(bray_endoS_PS_rarefied, MDS_data_bray_endoS_PS_rarefied$host_scientific_name,bias.adjust = T)
  anova(model_endo_disp) # significant dispersion
  
  # Re do a model with balanced samples sizes
  
Species_list = MDS_data_bray_endoS_PS_rarefied  %>% 
    group_by(host_scientific_name) %>% 
    summarise(n=n(),depth_endo=min(depth_endo)) %>% 
    subset(n>3) %>% pull(host_scientific_name)

Species_list
# 15 species with 4 individual each with >1000 endosymbiont reads


Info_depth %>% subset(host_scientific_name %in% Species_list) %>% group_by(host_scientific_name) %>%  summarise(depth_endo=mean(depth_endo),depth_gut=median(depth_gut),n=n())
  
Permanova_res = list()

for (i in 1:100)
{
  print(i)
  balanced_sample = MDS_data_bray_endoS_PS_rarefied  %>% 
    subset(host_scientific_name %in% Species_list) %>% 
    group_by(host_scientific_name) %>% 
    sample_n(4)
  
  balanced_sample %>% 
    group_by(host_scientific_name) %>% summarise(n=n())
  
  d = as.matrix(bray_endoS_PS_rarefied)[balanced_sample$sample_alias,balanced_sample$sample_alias]
  
  model_endo_balanced = adonis2(d ~ host_scientific_name + host_sex + Elevation ,
                                by="margin", data=balanced_sample)
  
  
  Permanova_res[[i]] = unnest(model_endo_balanced) %>% 
    subset(!is.na(F)) %>% 
    mutate(Predictor=c("Species","Sex" ,"Elevation"),
           Symbiont = rep("Endosymbiont",3)) %>% 
    mutate(Signif = ifelse(`Pr(>F)`>.05,"",ifelse(`Pr(>F)`>.01,"*",'**')),
           rep=i)
  print(Permanova_res[[i]])
}


Permanova_res_endo <- do.call(rbind,Permanova_res) %>% 
  group_by(Predictor) %>% 
  select(-rep) %>% 
  summarise_if(is.numeric, median, na.rm = TRUE)   %>% mutate(Symbiont="Endosymbiont")

Permanova_res_endo

  ################################
  #### Fig 1B. Gut symbiont   ####
  ################################
sample_gut_unique=c("MEP2-B7","MEP1-G7","MEP2-A4","MEP4-D12")# to remove

  bray_gut_PS_rarefied  <- prune_samples(!gutB_PS_rarefied@sam_data$sample_alias %in% sample_gut_unique,gutB_PS_rarefied) %>% 
    distance(betaM)
  hist(c(bray_gut_PS_rarefied ))
  
  MDS_bray_gut_PS_rarefied <- bray_gut_PS_rarefied %>% metaMDS() 
  MDS_data_bray_gut_rarefied <-  MDS_bray_gut_PS_rarefied$points %>% as_tibble(rownames = "sample_alias") %>% 
    left_join(Info_depth)  
  
  model_gutB_disp = betadisper(bray_gut_PS_rarefied, MDS_data_bray_gut_rarefied$host_scientific_name,bias.adjust = T)
  anova(model_gutB_disp)
  

  # Re do a model with balanced samples sizes
  
  Species_list =   MDS_data_bray_gut_rarefied   %>% 
    group_by(host_scientific_name) %>% 
    summarise(n=n(),depth_endo=min(depth_endo)) %>% 
    subset(n>3) %>% pull(host_scientific_name)
  
  Species_list
  # 15 species with 4 individual each with >1000 endosymbiont reads
  
  
  Info_depth %>% subset(host_scientific_name %in% Species_list) %>% group_by(host_scientific_name) %>%  summarise(depth_endo=mean(depth_endo),depth_gut=median(depth_gut),n=n())
  
  Permanova_resG = list()
  
  for (i in 1:100)
  {
    print(i)
    balanced_sample = MDS_data_bray_gut_rarefied  %>% 
      subset(host_scientific_name %in% Species_list) %>% 
      group_by(host_scientific_name) %>% 
      sample_n(4)
    
    balanced_sample %>% 
      group_by(host_scientific_name) %>% summarise(n=n())
    
    d = as.matrix(bray_gut_PS_rarefied)[balanced_sample$sample_alias,balanced_sample$sample_alias]
    
    model_gut_balanced = adonis2(d ~ host_scientific_name + host_sex + Elevation ,
                                  by="margin", data=balanced_sample)
    
    
    Permanova_resG[[i]] = unnest(model_gut_balanced) %>% 
      subset(!is.na(F)) %>% 
      mutate(Predictor=c("Species","Sex" ,"Elevation"),
             Symbiont = rep("Gut symbiont",3)) %>% 
      mutate(Signif = ifelse(`Pr(>F)`>.05,"",ifelse(`Pr(>F)`>.01,"*",'**')),
             rep=i)
    print(Permanova_resG[[i]])
  }
  
  
  Permanova_res_gut <- do.call(rbind,Permanova_resG) %>% 
    group_by(Predictor) %>% 
    select(-rep) %>% 
    summarise_if(is.numeric, median, na.rm = TRUE)   %>% mutate(Symbiont="Putative gut symbiont")
  
  Permanova_res_gut 
  
  
  final = rbind(Permanova_res_endo,Permanova_res_gut) %>% 
    mutate_if(is.numeric,round,3)
  

pdf("Supp_Data_in_SM_doc/SuppTable_Balanced_Permanova.pdf")       # Export PDF
grid.table(final)
dev.off()
  
  
  
  
