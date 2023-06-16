# Define plot features 
# 1 pt = 0.35mm
# Recommended max size=7 , min size =5 in pt 
P1=5;P2=6;P3=7
MyTheme=theme_classic() + 
  theme(plot.title = element_text(size=P3),
        axis.text=element_text(size=P2),
        axis.text.x =element_text(size=P2,angle=45,hjust=1),
        axis.title=element_text(size=P2),
        
        legend.text=element_text(size=P2),
        legend.title=element_text(size=P2),
        
        panel.border = element_blank(),
        
        axis.ticks = element_line(size = 1*.35),
        axis.ticks.length = unit(.5, "mm"),
        
        plot.caption = element_text(hjust = 0), 
        plot.title.position = "plot")