#------------------#
#  ASV prevalence. #
#------------------#

###################
## DATA AND CODE ##
###################


rm(list=ls())

# Packages 
library(microbiomeutilities)
library(circlize)
library(microbiome)
library(ComplexHeatmap)
library(randomForest)
source("Scripts/V_submitted/Revision/0_Data_preparation/0_Util_fonctions.R")
source("Scripts/V_submitted/Revision/0_Data_preparation/5_Load_Data.R")



#------------------#
#  ASV prevalence. #
#------------------#

ASV_info <- as.data.frame(tax_table(phyloseq_obj))  %>% 
  mutate(prev = prevalence(phyloseq_obj, count = T),
         type = ifelse(Genus %in% c("Spiroplasma","Wolbachia"),"Endosymbiont","Putative gut symbiont") )

ASV_info %>% group_by(type) %>% 
  summarise(mean_prev = mean(prev),
            sd_prev = sd(prev))

prev_plot <- ASV_info  %>% 
  ggplot(aes(x=prev, fill=type)) + scale_fill_manual(values=cols) + facet_wrap(vars(type)) + 
  geom_histogram() +   xlab("ASV prevalence (sample counts)") + ylab("ASV counts") + #ggtitle("A. ASV prevalence") +
  scale_x_log10() + scale_y_log10() + annotation_logticks() +
  MyTheme +   
  theme(legend.title = element_blank()) + 
  theme(legend.position="none")

prev_plot 
ggsave(plot = prev_plot, filename = "Redaction/Submission/ISMEcom_revision2/Supp_Data_in_SM_doc/Prevalence_plot.pdf",
       device="pdf", width = 120,height =80  ,units="mm")

