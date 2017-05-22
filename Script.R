#################################################
#                                               #
#         NEW Mimulus sp. MAP                   #
#                                               #
# Researchers                                   #
#   Leandro Mendonca                            #
#   Giovanni Galli                              #
#   Julia Morosini                              #
#   Humberto Fanelli                            #
#                                               #
#################################################

#Package
library(onemap)

#Taking data
data <- read_mapmaker(dir="C:\\biometriademarcadores\\",
                      file="m_feb06.raw")

# Defining parameters (grid)

combinations <- matrix(NA, length(3:12), length(seq(0.30, 0.5, by = 0.01)))
rownames(combinations) <- as.character(3:12)
colnames(combinations) <- as.character(seq(0.30, 0.5, by = 0.01))

# First for()
# Uses LOD from 3 to 12
# Second for()
# Uses r from 0.3 to 0.5
# Number of gruoups are estimated for every LOD x r combination within the specified range

for(lod in 3:12){
  for(r in seq(0.30, 0.5, by = 0.01)){
    TpT <- rf_2pts(data, LOD=lod, max.rf=r)
    Seq <- make_seq(TpT, arg="all")
    combinations[as.character(lod),as.character(r)] <- (group(Seq))$n.groups
  }
}

# number of groups (defined by LOD and r)
frame <- data.frame(LOD = as.numeric(rep(rownames(combinations),dim(combinations)[2])),
                    r = as.numeric(rep(colnames(combinations),each = dim(combinations)[1])),
                    Groups = factor(as.vector(combinations)))

ggplot(frame, aes(x=r,y=LOD,group=(Groups))) +
  geom_point(aes(color=Groups),show.legend = T)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  annotate("rect", xmin = 0.305, xmax = 0.315, ymin = 6.8, ymax = 7.2, alpha=.2)

suggest_lod(data)

#Defining parameters
n.mar <- data$n.mar
rmax <- 0.31 
LOD=7

#Two point test, considering all markers
TpT <- rf_2pts(data, LOD=LOD, max.rf=rmax)

#Get sequence of markers
Seq <- make_seq(TpT, arg="all")

#Formation od linkage groups
(LGall <- group(Seq))

# Kosambi map function
set_map_fun(type="kosambi") 

###############################################################
#        Markers orderin on each linkage group               #
###############################################################

###############################################################
#                    Linkage Group 1                          #
###############################################################

#take only markers of LG1
G1 <- make_seq(LGall, 1) 

#ordening
LG1_ord <- order_seq(input.seq=G1, n.init = 5,
                       subset.search = "twopt",
                       twopt.alg = "rcd", THRES = 3,
                       draw.try = F, wait = 1,
                       touchdown=TRUE)
LG1_ord

#Bilding final order and plot
LG1.final<-make_seq(LG1_ord, "force")
rf_graph_table(LG1.final)

###############################################################
#                    Linkage Group 2                          #
###############################################################

#take sequence
G2 <- make_seq(LGall,2)

#ordening
LG2_ord <- order_seq(input.seq=G2, n.init = 5,
                     subset.search = "twopt",
                     twopt.alg = "rcd", THRES = 3,
                     draw.try = FALSE, wait = 1,
                     touchdown=TRUE)

#Bilding final order and plot
LG2.final<-make_seq(LG2_ord, "force")
rf_graph_table(LG2.final)

###############################################################
#                    Linkage Group 3                          #
###############################################################

#take sequence
G3 <- make_seq(LGall,3)

#ordening
LG3_ord <- order_seq(input.seq=G3, n.init = 5,
                     subset.search = "twopt",
                     twopt.alg = "rcd", THRES = 3,
                     draw.try = FALSE, wait = 1,
                     touchdown=TRUE)

#Bilding final order and plot
LG3.final<-make_seq(LG3_ord, "force")
rf_graph_table(LG3.final)

###############################################################
#                    Linkage Group 4                          #
###############################################################

#take only markers of LG4
G4 <- make_seq(LGall, 4) 

#ordening
LG4_ord <- order_seq(input.seq=G4, n.init = 5,
                     subset.search = "twopt",
                     twopt.alg = "ug", THRES = 3,
                     draw.try = F, wait = 1,
                     touchdown=TRUE)

#Bilding final order and plot
LG4.final<-make_seq(LG4_ord, "force")
rf_graph_table(LG4.final) # too bed :`(

#Take out 77 and 126 
G4 <- drop_marker(G4, c(77,126))

#ordening
LG4_ord <- order_seq(input.seq=G4, n.init = 5,
                     subset.search = "twopt",
                     twopt.alg = "ug", THRES = 3,
                     draw.try = F, wait = 1,
                     touchdown=TRUE)

#Bilding final order and plot
LG4.final<-make_seq(LG4_ord, "force")
rf_graph_table(LG4.final) # better :) 

###############################################################
#                    Linkage Group 5                          #
###############################################################

#take sequence
G5 <- make_seq(LGall,5)

#ordening
LG5_ord <- order_seq(input.seq=G5, n.init = 5,
                     subset.search = "twopt",
                     twopt.alg = "rcd", THRES = 3,
                     draw.try = FALSE, wait = 1,
                     touchdown=TRUE)

#Bilding final order and plot
LG5.final<-make_seq(LG5_ord, "force")
rf_graph_table(LG5.final)

###############################################################
#                    Linkage Group 6                          #
###############################################################

#take only markers of LG6
G6 <- make_seq(LGall, 6) 

#ordening
LG6_ord <- order_seq(input.seq=G6, n.init = 5,
                     subset.search = "twopt",
                     twopt.alg = "ser", THRES = 3,
                     draw.try = F, wait = 1,
                     touchdown=TRUE)

#Bilding final order and plot
LG6.final<-make_seq(LG6_ord, "force")
rf_graph_table(LG6.final, scale=1.5) #not good :(

#Take out 150, 213, 312 and 204
G6 <- drop_marker(G6, c(150, 213, 312, 204))

#ordening
LG6_ord <- order_seq(input.seq=G6, n.init = 5,
                     subset.search = "twopt",
                     twopt.alg = "ug", THRES = 3,
                     draw.try = F, wait = 1,
                     touchdown=TRUE)
LG6_ord

#Bilding final order and plot
LG6.final<-make_seq(LG6_ord, "force")
rf_graph_table(LG6.final, scale = 1.5) # better :)

###############################################################
#                    Linkage Group 7                          #
###############################################################

#We have 60 markers for this linkage group, lest try to devide tham using a high LOD

#take only markers of LG7
G7 <- make_seq(LGall, 7) 

#ordening
LG7_ord <- order_seq(input.seq=G7, n.init = 5,
                     subset.search = "twopt",
                     twopt.alg = "ug", THRES = 3,
                     draw.try = F, wait = 1,
                     touchdown=TRUE)

#Bilding final order and plot
LG7.final<-make_seq(LG7_ord, "force")
rf_graph_table(LG7.final ,scale=1.5) #too bed :( 

#New senquence and group formation without 405 and 265
G7 <- drop_marker(G7, c(405,256))
(G7 <- group(G7))

#Now 2 groups was formed !!!
#Defining markers order

#take sequences
G7.1 <- make_seq(G7, 1) 
G7.2 <- make_seq(G7, 2)

#ordening
LG7.1_ord <- order_seq(input.seq=G7.1, n.init = 5,
                       subset.search = "twopt",
                       twopt.alg = "ug", THRES = 3,
                       draw.try = F, wait = 1,
                       touchdown=TRUE)

LG7.2_ord <- order_seq(input.seq=G7.2, n.init = 5,
                       subset.search = "twopt",
                       twopt.alg = "ser", THRES = 3,
                       draw.try = F, wait = 1,
                       touchdown=TRUE)


#Bilding final order and plot
LG7.1.final<-make_seq(LG7.1_ord, "force")
rf_graph_table(LG7.1.final, scale = 1.5)

LG7.2.final<-make_seq(LG7.2_ord,"force")
rf_graph_table(LG7.2.final, scale = 1.5)

###############################################################
#                    Linkage Group 8                          #
###############################################################

#take sequence
G8 <- make_seq(LGall,8)

#ordening
LG8_ord <- order_seq(input.seq=G8, n.init = 5,
                     subset.search = "twopt",
                     twopt.alg = "rcd", THRES = 3,
                     draw.try = FALSE, wait = 1,
                     touchdown=TRUE)

#Bilding final order and plot
LG8.final<-make_seq(LG8_ord, "force")
rf_graph_table(LG8.final)

###############################################################
#                    Linkage Group 9                          #
###############################################################

#take sequence
G9 <- make_seq(LGall,9)

#ordening
LG9_ord <- order_seq(input.seq=G9, n.init = 5,
                     subset.search = "twopt",
                     twopt.alg = "ug", THRES = 3,
                     draw.try = FALSE, wait = 1,
                     touchdown=TRUE)

#Bilding final order and plot
LG9.final<-make_seq(LG9_ord, "force")
rf_graph_table(LG9.final)

###############################################################
#                    Linkage Group 10                        #
###############################################################

#take sequence
G10 <- make_seq(LGall,10)

#ordening
LG10_ord <- order_seq(input.seq=G10, n.init = 5,
                      subset.search = "twopt",
                      twopt.alg = "ser", THRES = 3,
                      draw.try = FALSE, wait = 1,
                      touchdown=TRUE)

#Bilding final order and plot
LG10.final<-make_seq(LG10_ord, "force")
rf_graph_table(LG10.final)

###############################################################
#                    Linkage Group 11                         #
###############################################################

#take sequence
G11 <- make_seq(LGall,11)

#ordening
LG11_ord <- order_seq(input.seq=G11, n.init = 5,
                      subset.search = "twopt",
                      twopt.alg = "ser", THRES = 3,
                      draw.try = FALSE, wait = 1,
                      touchdown=TRUE)

#Bilding final order and plot
LG11.final<-make_seq(LG11_ord, "force")
rf_graph_table(LG11.final)

###############################################################
#                    Linkage Group 12                        #
###############################################################

#take sequence
G12 <- make_seq(LGall,12)
G12<- drop_marker(G12, c(176, 182, 194, 203, 205, 242))

#ordening
LG12_ord <- order_seq(input.seq=G12, n.init = 5,
                      subset.search = "twopt",
                      twopt.alg = "rcd", THRES = 3,
                      draw.try = FALSE, wait = 1,
                      touchdown=TRUE)

#Bilding final order and plot
LG12.final<-make_seq(LG12_ord, "force")
rf_graph_table(LG12.final, scale = 1.5)#

###############################################################
#                    Linkage Group 13                         #
###############################################################

#take sequence
G13 <- make_seq(LGall,13)

#ordening
LG13_ord <- order_seq(input.seq=G13, n.init = 5,
                      subset.search = "twopt",
                      twopt.alg = "rcd", THRES = 3,
                      draw.try = FALSE, wait = 1,
                      touchdown=TRUE)

#Bilding final order and plot
LG13.final<-make_seq(LG13_ord, "force")
rf_graph_table(LG13.final)

############## Heatmaps #########

rf_graph_table(LG1.final, scale=2, main ="LG1", inter=F)   #LG1
rf_graph_table(LG2.final, scale=2, main ="LG2",inter=F)    #LG2
rf_graph_table(LG3.final, scale=2, main ="LG3", inter=F)   #LG3
rf_graph_table(LG4.final, scale=2, main ="LG4", inter=F)   #LG4
rf_graph_table(LG5.final, scale=2, main ="LG5", inter=F)   #LG5
rf_graph_table(LG6.final, scale=2, main ="LG6", inter=F)   #LG6
rf_graph_table(LG7.1.final, scale=2, main ="LG7", inter=F) #LG7
rf_graph_table(LG8.final, scale=2, main ="LG8", inter=F) #LG8
rf_graph_table(LG9.final, scale=2, main ="LG9", inter=F)   #LG9
rf_graph_table(LG7.2.final, scale=2, main ="LG10", inter=F)  #LG10
rf_graph_table(LG11.final, scale=2, main ="LG11", inter=F) #LG11
rf_graph_table(LG12.final, scale=2, main ="LG12", inter=F) #LG12
rf_graph_table(LG13.final, scale=2, main ="LG13", inter=F) #LG13
rf_graph_table(LG10.final, scale=2, main ="LG14", inter=F) #LG14

### MapChart ###

map.list = list(LG1.final, LG2.final, LG3.final, LG4.final,
                LG5.final, LG6.final, LG7.1.final, LG8.final,
                LG9.final, LG7.2.final, LG11.final, LG12.final,
                LG13.final, LG10.final)
write_map(map.list, file.out = "LGs.txt")
(map.chart = read.table("LGs.txt"))
attach(map.chart)
map.chart$V3 = round(map.chart$V3, 2)
map.chart = subset(map.chart, select = c(2,3))
head(map.chart)
write.table(map.chart, "map.chart.txt")

### Export Heatmaps ###

tiff("LG1.tif", width=1500, height=1500, units="px", res=250)
rf_graph_table(LG1.final, scale=2, main ="LG1", inter=F)  
dev.off()

tiff("LG2.tif", width=1500, height=1500, units="px", res=250)
rf_graph_table(LG2.final, scale=2, main ="LG2", inter=F)  
dev.off()

tiff("LG3.tif", width=1500, height=1500, units="px", res=250)
rf_graph_table(LG3.final, scale=2, main ="LG3", inter=F)  
dev.off()

tiff("LG4.tif", width=1500, height=1500, units="px", res=250)
rf_graph_table(LG4.final, scale=2, main ="LG4", inter=F)  
dev.off()

tiff("LG5.tif", width=1500, height=1500, units="px", res=250)
rf_graph_table(LG5.final, scale=2, main ="LG5", inter=F)  
dev.off()

tiff("LG6.tif", width=1500, height=1500, units="px", res=250)
rf_graph_table(LG6.final, scale=2, main ="LG6", inter=F)  
dev.off()

tiff("LG7.tif", width=1500, height=1500, units="px", res=250)
rf_graph_table(LG7.1.final, scale=2, main ="LG7", inter=F)  
dev.off()

tiff("LG8.tif", width=1500, height=1500, units="px", res=250)
rf_graph_table(LG8.final, scale=2, main ="LG8", inter=F)  
dev.off()

tiff("LG9.tif", width=1500, height=1500, units="px", res=250)
rf_graph_table(LG9.final, scale=2, main ="LG9", inter=F)  
dev.off()

tiff("LG10.tif", width=1500, height=1500, units="px", res=250)
rf_graph_table(LG7.2.final, scale=2, main ="LG10", inter=F)  
dev.off()

tiff("LG11.tif", width=1500, height=1500, units="px", res=250)
rf_graph_table(LG11.final, scale=2, main ="LG11", inter=F)  
dev.off()

tiff("LG12.tif", width=1500, height=1500, units="px", res=250)
rf_graph_table(LG12.final, scale=2, main ="LG12", inter=F)  
dev.off()

tiff("LG13.tif", width=1500, height=1500, units="px", res=250)
rf_graph_table(LG13.final, scale=2, main ="LG13", inter=F)  
dev.off()

tiff("LG14.tif", width=1500, height=1500, units="px", res=250)
rf_graph_table(LG10.final, scale=2, main ="LG14", inter=F)  
dev.off()
