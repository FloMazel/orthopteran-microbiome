###################################################################
# Study individual ASV distribution across samples and host species 
###################################################################

rm(list=ls())

# Packages 
library(microbiomeutilities)
library(circlize)
library(microbiome)
library(ComplexHeatmap)
library(randomForest)
source("0_Data_preparation/0_Util_fonctions.R")
source("0_Data_preparation/5_Load_Data.R")

# Random forest = function to compute oob on random data
null_oob = function(i,rf.data, ntree){
  print(i)
  rf.data.random <- rf.data
  rf.data.random$response <- sample(rf.data$response)
  endo.classify <- randomForest(response~., data = rf.data.random, ntree = ntree)
  oob <- endo.classify$err.rate[ntree,1]
  return(1-oob)
}


# Parameters and initial data (phyloseq object loaded from "Load_Data" script)

set.seed(2)
ntree = 100 # trees in randomforest
transformReads="compositional" # transfo of read counts
min_depth = 999 # minimum depth to keep a sample
subset.top.gut = 50 # number of most abundant ASVSs kept for heatmap
subset.top.endo = 50 # number of most abundant ASVSs kept for heatmap

# Gut ASV data
gut_PS <- subset_samples(gutB_PS,sample_sums(gutB_PS)>min_depth)
gut_PS <-  microbiome::transform(gut_PS, transformReads)
gut_PS # 1700 ASVS in 167 samples
otu.mat_gut <- abundances(gut_PS); meta.tab_gut <- meta(gut_PS); taxo.tab_gut <- tax_tibble(gut_PS)

topOTU_gut <- microbiome::top_taxa(gutB_PS, n = subset.top.gut)
gut_PS_HM <- phyloseq::prune_taxa(topOTU_gut,gut_PS)
gut_PS_HM # 100 ASVS in 167 samples

# Endo ASV data
endoS_PS <- subset_samples(endoS_PS,sample_sums(endoS_PS)>min_depth)
endoS_PS <-  microbiome::transform(endoS_PS, transformReads)
endoS_PS # 257 ASVS in 259 samples
otu.mat_endo <- abundances(endoS_PS);meta.tab_endo <- meta(endoS_PS); taxo.tab_end <- tax_tibble(endoS_PS)

topOTU_endoS <- microbiome::top_taxa(endoS_PS, n = subset.top.endo)
endoS_PS_HM <- prune_taxa(topOTU_endoS,endoS_PS)
endoS_PS_HM # 50 ASVS in 259 samples

# Colors
cols_endo <- list(Taxonomy=c('Wolbachia'="#FFD92F",'Spiroplasma'="chocolate3", 'Others endosymbionts'="grey"))
colsGutOrder <- list(Taxonomy=c("Xanthomonadales"="#0072B2","Enterobacterales"="tan4",
                                "Pseudomonadales"="#D55E00","Diplorickettsiales"="#009E73",
                                "Lactobacillales"="darkblue","Peptostreptococcales-Tissierellales"="purple",
                                "Sphingomonadales"="#E69F00","Rhizobiales"='darkgreen',
                                "Mycoplasmatales"="black"))



#------------------------------------------#
#    Random forest.        analysis       #
#------------------------------------------#

# Endosymbionts
#-------------

# RUN MODELS
predictors <- t(otu.mat_endo); dim(predictors)
response <- as.factor(meta.tab_endo$host_scientific_name)
rf.data.endo <- data.frame(response, predictors)

endo.classify <- randomForest(response~., data = rf.data.endo, ntree = ntree)
obs_oob_endo <- endo.classify$err.rate[ntree,1]
print(obs_oob_endo)

null_oob_endo <- sapply(1:99,FUN = null_oob, rf.data.endo, ntree)

# plots

endo_oob_plot <- tibble(null=null_oob_endo) %>%
  ggplot(aes(x=null)) + geom_histogram() + MyTheme + 
  xlim(0,1) + geom_vline(xintercept = 1-obs_oob_endo, colour="#FFD92F") + 
  xlab("Accuracy of random forest model (1-Out of Bag Error)") + ylab("Count") +
  annotate("text",x=.3,y=30, label="Randomized\n data", size=P1/.pt) + 
  annotate("text",x=.85,y=30, label="True data", colour="#FFD92F", size=P1/.pt) + 
  ggtitle("A. Endosymbionts")
endo_oob_plot

pval_endo=sum(null_oob_endo>(1-obs_oob_endo))/100
pval_endo

# Importance scores
imp_endo <- importance(endo.classify)
imp_endo <- data.frame(predictors = rownames(imp_endo), imp_endo)


# Gut symbionts
#-------------

# RUN MODELS
predictors <- t(otu.mat_gut); dim(predictors)
response <- as.factor(meta.tab_gut$host_scientific_name)
rf.data.gut <- data.frame(response, predictors)

gut.classify <- randomForest(response~., data = rf.data.gut, ntree = ntree)
obs_oob_gut <- gut.classify$err.rate[ntree,1]
obs_oob_gut

null_oob_gut <- sapply(1:199,FUN = null_oob, rf.data.gut, ntree)

# Plot
gut_oob_plot <- tibble(null=null_oob_gut) %>%
  ggplot(aes(x=null)) + geom_histogram() + MyTheme + 
  xlim(0,1) + geom_vline(xintercept = 1-obs_oob_gut, colour="#6ab41f") + 
  xlab("Accuracy of random forest model (1-Out of Bag Error)") + ylab("Count") +
  annotate("text",x=.05,y=30, label="Randomized\n data", size=P1/.pt) + 
  annotate("text",x=.5,y=30, label="True data", colour="#6ab41f", size=P1/.pt) + 
  ggtitle("B. Putative gut symbionts")
gut_oob_plot

pval_gut=sum(null_oob_gut>(1-obs_oob_gut))/100
pval_gut

# Importance scores
imp_gut  <- importance(gut.classify)
imp_gut  <- data.frame(predictors = rownames( imp_gut ),  imp_gut )


# Combining plots
#-------------

oobplot <- plot_grid(endo_oob_plot,gut_oob_plot, ncol = 2)
ggsave(plot = oobplot, filename = "Redaction/Submission/ISMEcom_revision2/Supp_Data_in_SM_doc/OOB_plot.pdf",
       device="pdf", width = 140,height =60  ,units="mm")



#------------------#
#    Heat map      #
#------------------#

ht_opt(HEATMAP_LEGEND_PADDING=unit(10,"mm"),
       ANNOTATION_LEGEND_PADDING=unit(10,"mm"))

# Endosymbionts
#-------------

# extract and prepare relevant object for the heatmap
otu.mat_endo <- abundances(endoS_PS_HM); meta.tab_endo <- meta(endoS_PS_HM); 
taxo.tab_endo <- tax_tibble(endoS_PS_HM) %>% 
left_join(imp_endo, by=c("FeatureID"="predictors")) %>% 
mutate(Genus_plot = ifelse(Genus %in% c('Wolbachia','Spiroplasma'),Genus,"Others endosymbionts"))

table(taxo.tab_endo$Genus)

col_fun = colorRamp2(c(min(c(otu.mat_endo)), max(c(otu.mat_endo))), c("white", "#000080"))



# Prepare the heatmap
column_ha = HeatmapAnnotation(Host = meta.tab_endo$Host_species, 
                              annotation_name_side = "right", col=colsHost,
                              annotation_label = "Host \ntaxonomy",
                              annotation_name_gp = gpar(fontsize = 7),
                              annotation_legend_param = list(title_gp = gpar(fontsize = 10), 
                                                             labels_gp = gpar(fontsize = 8, fontface="italic")))

row_ha = rowAnnotation(Taxonomy = taxo.tab_endo$Genus_plot, 
                       col = cols_endo, show_annotation_name = T,
                       Importance = anno_barplot(taxo.tab_endo$MeanDecreaseGini,
                                                 axis_param=list(at=c(0,10,20), labels_rot=45)),
                       annotation_name_rot = 45,
                       annotation_name_gp = gpar(fontsize = 7),
                       annotation_label = c("ASV taxo.", "Importance score \nto classify host sp."), 
                       annotation_name_offset = unit(c(7,7),"mm"),
                       annotation_legend_param = list(title_gp = gpar(fontsize = 10), 
                                                      labels_gp = gpar(fontsize = 8)))


endoHM <- Heatmap(otu.mat_endo, name="Read counts \n(proportion)", col = col_fun,
                  column_dend_side = "bottom",
                  show_row_names = F, show_column_names = F, show_row_dend = F,
                  row_title = "Bacterial ASVs", row_title_gp = gpar(fontsize = 10),
                  column_title = "A. Endosymbionts", column_title_gp = gpar(fontsize = 14, fontface = "bold"),
                  bottom_annotation = column_ha, left_annotation = row_ha,
                  width = unit(6, "cm"),
                  heatmap_legend_param = list(title_gp = gpar(fontsize = 10),
                                              labels_gp = gpar(fontsize = 8)))
endoHM


# Putative gut symbionts
#-------------------------

# extract and prepare relevant object for the heatmap
otu.mat_gut <- abundances(gut_PS_HM); 
meta.tab <- meta(gut_PS_HM) 
taxo.tab_gut <- tax_tibble(gut_PS_HM) %>% 
  mutate(FeatureID = gsub("-",".",FeatureID)) %>% 
  left_join(imp_gut, by=c("FeatureID"="predictors"))

unique(taxo.tab_gut$Order)

col_fun = colorRamp2(c(min(c(otu.mat_gut)), max(c(otu.mat_gut))), c("white", "#000080"))

# Prepare the heatmap
length(unique(unique(meta.tab_gut$Host_genus)))

column_ha = HeatmapAnnotation(Host = meta.tab$Host_species,  show_annotation_name = F, 
                              annotation_label = "Host \ntaxonomy",
                              annotation_name_gp = gpar(fontsize = 7),
                              col=colsHost,
                              annotation_legend_param = list(title_gp = gpar(fontsize = 10), 
                                                             labels_gp = gpar(fontsize = 8, fontface="italic")))

row_ha = rowAnnotation(Taxonomy = taxo.tab_gut$Order, 
                       show_annotation_name = T, col=colsGutOrder,
                       Importance = anno_barplot(taxo.tab_gut$MeanDecreaseGini, 
                                                 axis_param=list(at=c(0,5), labels = c(0,1), labels_rot=45)),
                       annotation_name_rot = 45,
                       annotation_name_gp = gpar(fontsize = 7),
                       annotation_label = c("ASV taxo.", "Importance score \nto classify host sp."), 
                       annotation_name_offset = unit(c(7,7),"mm"),
                       annotation_legend_param = list(title_gp = gpar(fontsize = 10), 
                                                      labels_gp = gpar(fontsize = 8)))


gutHM <- Heatmap(otu.mat_gut, name="Read counts \n(proportion)", col = col_fun,
                 column_dend_side = "bottom",
                 show_row_names = F, show_column_names = F, show_row_dend = F,
                 row_title = "Bacterial ASVs", row_title_gp = gpar(fontsize = 10),
                 column_title = "B. Putative gut symbionts", column_title_gp = gpar(fontsize = 14, fontface = "bold"),
                 bottom_annotation = column_ha, 
                 left_annotation = row_ha,
                 width = unit(6, "cm"),
                 heatmap_legend_param = list(labels_gp = gpar(fontsize = 10)))

gutHM 

# Combining heatmaps
#--------------------

finalHM <- endoHM + gutHM

pdf("Main_Figure/HeatMap.pdf", height = 5, width = 12)
draw(finalHM,ht_gap = unit(2, "cm"))
dev.off()


