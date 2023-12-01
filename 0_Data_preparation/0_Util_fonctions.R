# Define plot features 
# 1 pt = 0.35mm
# Recommended max size=7 , min size =5 in pt 


colsHost <- list(Host=c("Anonconotus alpinus" = '#DDAD4B',
                        "Chorthippus biguttulus"= 'forestgreen', 
                        "Chorthippus parallelus"= 'yellowgreen', 
                        "Arcyptera fusca" = 'darkgray', 
                        "Decticus verrucivorus"   = 'cornflowerblue', 
                        "Euthystira brachyptera" = 'darkolivegreen4', 
                        "Gomphocerippus rufus" = 'indianred1', 
                        "Mecostethus parapleurus" = 'tan4', 
                        "Metrioptera roeselii" = 'darkblue', 
                        "Metrioptera saussuriana" = 'cadetblue1',  
                        "Oedipoda caerulescens" = 'lightsalmon',
                        "Omocestus rufipes" = "tan1",
                        "Omocestus viridulus" = 'orange', 
                        "Phaneroptera falcata" = 'wheat4', 
                        "Pholidoptera griseoaptera"= 'black', 
                        "Platycleis albopunctata"= 'moccasin', 
                        "Podisma pedestris"= 'mediumvioletred', 
                        "Polysarcus denticauda" = 'seagreen',
                        "Psophus stridulus" = 'firebrick4',
                        "Ruspolia nitidula"  = "tomato3" ,
                        "Stauroderus scalaris" = "gainsboro",
                        "Stenobothrus lineatus" = "purple",
                        "Tettigonia cantans" = "chocolate3",
                        "Miramella alpina"= "red"))

P1=5;P2=6;P3=7
MyTheme=theme_classic() + 
  theme(plot.title = element_text(size=P3, face="bold"),
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
