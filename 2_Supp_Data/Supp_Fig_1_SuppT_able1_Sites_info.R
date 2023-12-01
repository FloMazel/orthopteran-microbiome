################################################################################
################################################################################
# FSup Table 1
################################################################################
################################################################################


###################
## DATA AND CODE ##
###################


# Packages 
library(tidyverse)
library(proj4)
library(gridExtra)
library(xlsx)

# get a map of the region
library(ggmap)
library(ggsn)

# Tally table for species, site, sex. 
phyloseq_obj = readRDS("Phyloseq_object.RDS")
metaD = as(sample_data(phyloseq_obj), "data.frame")

tally_table = metaD %>% group_by(host_scientific_name,host_sex,sample_site_ID) %>% 
  summarise(n_samples = n())

write_csv(x = tally_table,file = "Tally_table_sampling.csv")

# S

# Plot_info
plot_info <- read.table("3_Metadata/1_site.info.txt",header = T) %>% 
  select(-Elevation_category)

# Projection
proj4string <- "+proj=somerc +lat_0=46.95240555555556 +lon_0=7.439583333333333 +k_0=1 +x_0=600000 +y_0=200000 +ellps=bessel +towgs84=674.4,15.1,405.3,0,0,0,0 +units=m +no_defs "

# Transformed data
pj <- project(select(plot_info,Coordinates.E,Coordinates.N), proj4string, inverse=TRUE)
plot_info$Lat=pj$y; plot_info$Long=pj$x

#Export Table
pdf("Supp_figures/SuppTable_1_Sites_info.pdf",height = 5,width = 10)       # Export PDF
grid.table(plot_info)
dev.off()



# Species and Site infos 
sample_info <-  read.table("sample_info.txt") 

Site_species_Summary <- sample_info  %>%  
  group_by(Species,Site) %>% 
  summarise(n_sample=n(),
            Suborder=unique(Suborder)) %>% 
  left_join(plot_info)

Site_species_plot <- Site_species_Summary  %>%  
  ggplot(aes(y=Elevation, x=reorder(Species,desc(Elevation)),size=n_sample,col=Site,shape=Suborder))  + 
  geom_point() + xlab("Host Species")+
  MyTheme + ggtitle("A")

Site_species_plot



mapImageData1 <- get_map(location = c(left = .99*min(plot_info$Long),
                                      right = 1.01*max(plot_info$Long),
                                      bottom = .999* min(plot_info$Lat),
                                      top = 1.001*max(plot_info$Lat)),
                         color = "color",
                         source = "osm",
                         maptype = "terrain")

map_sites <- ggmap(mapImageData1,
                   extent = "normal",
                   ylab = "Latitude",
                   xlab = "Longitude") +
  geom_point(data = plot_info, aes(x = Long, y = Lat,col=Site)) + 
  labs(x="Longitude", y ="Lattitude")

map_sites <- map_sites +   scalebar(transform = TRUE, dist_unit ="km",dist = 4,
                       x.min = .99*min(plot_info$Long),
                        x.max =  1.01*max(plot_info$Long),
                         y.min = .999* min(plot_info$Lat),
                         y.max = 1.001*max(plot_info$Lat)) + ggtitle("B")
  



fig_sites = plot_grid(Site_species_plot,map_sites,
                 ncol=1)
fig_sites

ggsave(filename = "Supp_figures/Supp_Fig_1_Sites.pdf",fig_sites,device = "pdf",
       width = 175,height =310  ,units="mm")






