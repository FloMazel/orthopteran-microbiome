# Load script from Figure 1 first 

################################################################################
################################################################################
# SEX and SPECIES  effect (Supp. material)
################################################################################
################################################################################

Summary_Categ_Sample = Taxo_counts_Genus %>% 
  mutate(Endosymbiont = if_else(Genus=="Spiroplasma"|Order=="Rickettsiales"|Order=="Chlamydiales","Endosymbiont","Gut symbiont", missing="Gut symbiont")) %>% 
  group_by(Endosymbiont,Sample) %>% 
  summarise(RRC=sum(Abundance),Suborder=unique(suborder_host),Species=unique(host_scientific_name),Sex=unique(host_sex)) %>% 
  pivot_wider(names_from = Endosymbiont,values_from=RRC) %>% 
  mutate(Ratio_EG = Endosymbiont/`Gut symbiont`)

SexSpeciesSymbiont_plot = Summary_Categ_Sample %>% 
  ggplot(aes(y=Endosymbiont,x=Species, fill=Sex)) + geom_boxplot() + facet_wrap(vars(Suborder),scales = "free_x")+
  MyTheme + ylab("Relative read counts of endosymbionts")

SexSpeciesSymbiont_plot
ggsave("Supp_Data/Supp_Figure_SexSpeciesSymbiont_plot.pdf",SexSpeciesSymbiont_plot,device = 'pdf',
       width = 175,height =80  ,units="mm")


Ensi = Summary_Categ_Sample %>% 
  subset(Suborder=="Ensifera")

EnsiSp = kruskal.test(Ensi$Endosymbiont,Ensi$Species)
EnsiSex = kruskal.test(Ensi$Endosymbiont,Ensi$Sex)

Caeli = Summary_Categ_Sample %>% 
  subset(Suborder=="Caelifera")

CaeliSp = kruskal.test(Caeli$Endosymbiont,Caeli$Species)
CaeliSex = kruskal.test(Caeli$Endosymbiont,Caeli$Sex)


KruskalTests = as.tibble(rbind(EnsiSp,EnsiSex,CaeliSp,CaeliSex)) %>% 
  select(statistic, parameter, p.value) %>% 
  unnest() %>% 
  mutate(statistic=round(statistic,2), 
         p.value=round(p.value,3),
         p.value=ifelse(p.value==0,"<.001",p.value)) %>% 
  mutate(Suborder=c("Ensifera","Ensifera","Caelifera","Caelifera"),
         Factor=c("Host species","Sex","Host species","Sex"))


EnsiSp = lm(Endosymbiont ~ Species*Sex , data=Ensi)
EnsiSp
#plot(EnsiSp)
EnsiMod = anova(EnsiSp) %>% as_tibble(rownames = "Factor") %>%
  mutate(Suborder=rep("Ensifera",4)) %>% mutate_if(is.numeric, round,3)
  

CaeliSp = lm(Endosymbiont ~ Species*Sex , data=Caeli)
CaeliSp 
#plot(CaeliSp )
CaeliMod = anova(CaeliSp) %>% as_tibble(rownames = "Factor") %>%
  mutate(Suborder=rep("Caelifera",4)) %>% mutate_if(is.numeric, round,3)
CaeliMod 


pdf("Supp_Data/SuppTable_4.pdf")       # Export PDF
grid.table(rbind(EnsiMod, CaeliMod))
dev.off()








