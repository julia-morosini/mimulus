##################################
##                              ##
##   QTL Mapping - IM and CIM   ##
##                              ##
##################################

library(qtl)
library(onemap)
library(knitr)
library(kableExtra)
library(reshape2)

# Importing
data <- read_mapmaker(dir=getwd(), file="m_feb06.raw")

mapaMim <- read.cross(format = "mm", file = "m_feb06.raw",
                        mapfile = "LGs.txt",map.function = "kosambi") # esse LG foi do mapchart

# Conditional genotype probabilities
probMim <- calc.genoprob(mapaMim, step = 1, map.function = "kosambi")
simMim <- sim.geno(probMim, n.draws = 5, map.function = "kosambi")

###############################################################################################
###                                 INTERVAL MAPPING (IM)                                   ###
###############################################################################################

#### Threshold ####

#em

perIM_em <- scanone(simMim, pheno.col = 8, method= "em", model = "normal",n.perm = 10000, verbose = T)
summary(perIM_em,alpha=.05)

ggplot(as.data.frame(perIM_em), aes(x=perIM_em[,1])) + 
  xlim(1,7) + xlab("Maximum LOD") + ylab("Frequency") +
  ggtitle("IM - em") + theme(plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5, size=13)) +
  geom_histogram(binwidth = 0.1, color="white", fill="#000099") +
  geom_vline(xintercept = quantile(perIM_em[,1], 0.95), linetype="solid", color="red", size=0.7) +
  annotate("text", x = 5.8, y = 150, label = "Threshold = 3.86", size = 3.5)

#hk

perIM_hk <- scanone(simMim, pheno.col = 8, method= "hk", model = "normal",n.perm = 10000, verbose = F)
summary(perIM_hk,alpha=.05)
#plot(perIM_hk, col = "orange", breaks = 50)
#abline(v=summary(perIM_hk,alpha=.05)[1], col = "black", lwd = 2)

ggplot(as.data.frame(perIM_hk), aes(x=perIM_hk[,1])) + 
  xlim(1,7) + xlab("Maximum LOD") + ylab("Frequency") +
  ggtitle("IM - Halley-Knott") + theme(plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5, size=13)) +
  geom_histogram(binwidth = 0.1, color="white", fill="#009E73") +
  geom_vline(xintercept = quantile(perIM_hk[,1], 0.95), linetype="solid", color="red", size=0.6) +
  annotate("text", x = 6, y = 150, label = "Threshold = 3.83", size = 3.5)

#### Mapping ####

#em
IMscan_em <- scanone(simMim, pheno.col = 8, method = "em", model = "normal")
summary.IMscan.em <- summary(IMscan_em, perms=perIM_em, alpha = 0.05) # sign QTL

sumpos <- data.frame(chr=IMscan_em[,1], pos=cumsum(IMscan_em[,2]))
a <- as.vector(tapply(sumpos$pos, sumpos$chr, mean))
bpLabels <- paste(1:14)

ggplot(as.data.frame(IMscan_em), aes(cumsum(IMscan_em[,2]), IMscan_em[,3])) + 
  geom_line() + geom_area(fill=alpha('blue',0.2)) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  geom_hline(yintercept = 3.86, linetype="solid", color="red", size=0.6) +
  ggtitle("IM - em") + theme(plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5, size=13)) +
  scale_x_continuous(breaks=a, labels=bpLabels) + xlab("Chromosome") + ylab("LOD")

#hk
IMscan_hk <- scanone(simMim, pheno.col = 8, method = "hk", model = "normal")
summary.IMscan.hk <- summary(IMscan_hk, perms=perIM_hk, alpha = 0.05) # sign QTL

ggplot(as.data.frame(IMscan_hk), aes(cumsum(IMscan_hk[,2]), IMscan_hk[,3])) + 
  geom_line() + geom_area(fill=alpha('red',0.2)) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  geom_hline(yintercept = 3.83, linetype="solid", color="red", size=0.6) +
  ggtitle("IM - Halley-Knott") + theme(plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5, size=13)) +
  scale_x_continuous(breaks=a, labels=bpLabels) + xlab("Chromosome") + ylab("LOD")

#plot em x hk
ggplot(NULL, aes(cumsum(pos), lod, fill=cumsum(pos))) + 
  geom_line(data = as.data.frame(IMscan_em), colour="blue", size=1.5) + ylim(0,9) +
  geom_step(data = as.data.frame(IMscan_hk), colour="magenta1", size=0.8) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  scale_x_continuous(breaks=a, labels=bpLabels) + xlab("Chromosome") + ylab("LOD") +
  geom_hline(yintercept = 3.86, linetype="solid", color="red", size=0.6) +
  geom_hline(yintercept = 3.83, linetype="dashed", color="black", size=0.6) +
  ggtitle("IM - em x hk") + theme(plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5, size=15))

# plot difference LOD between methods
plot(IMscan_em - IMscan_hk, col = c("blue"))

ggplot(as.data.frame(IMscan_em[,3] - IMscan_hk[,3]), aes(cumsum(IMscan_hk[,2]), IMscan_em[,3] - IMscan_hk[,3])) + 
  geom_line() + geom_area(fill=alpha('red',0.2)) + ylim(-0.25,0.9) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  ggtitle("Difference LOD in IM: em - hk") + theme(plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5, size=13)) +
  scale_x_continuous(breaks=a, labels=bpLabels) + xlab("Chromosome") + ylab("LOD")

# plot principal QTLs
plot(IMscan_em, IMscan_hk, chr=c(6,13,14), col=c("blue", "magenta"), incl.markers = T, show.marker.names = F, bandcol="lightgray")

# QTLs effect
max(IMscan_hk)
mar <- find.marker(simMim, chr=c(13), pos=c(58))
plotPXG(simMim, pheno.col = 8, marker=mar) # pheno x geno do QTL de maior efeito
effectplot(simMim, pheno.col = 8, mname1=mar)

# Tabelas 1.1("hk")

tabela.im = data.frame("Linkage Group"=c(summary.IMscan.hk$chr[1],summary.IMscan.hk$chr[2],summary.IMscan.hk$chr[3]),
                       "Position"=c(summary.IMscan.hk$pos[1],summary.IMscan.hk$pos[2],summary.IMscan.hk$pos[3]),
                       "LOD"=c(round(summary.IMscan.hk$lod[1],2),round(summary.IMscan.hk$lod[2],2),round(summary.IMscan.hk$lod[3],2)))
                       
kable(tabela.im)

###############################################################################################
###                            COMPOSITE INTERVAL MAPPING (CIM)                             ###
###############################################################################################

#### Threshold ####

#em

perCIM_em <- cim(cross = simMim, pheno.col = 8, window = 10,n.marcovar = 3, method = "em", n.perm = 1000)
summary(perCIM_em,alpha=.05)

ggplot(as.data.frame(perCIM_em), aes(x=perCIM_em[,1])) + 
  xlim(1,8) + xlab("Maximum LOD") + ylab("Frequency") +
  ggtitle("CIM - em") + theme(plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5, size=13)) +
  geom_histogram(binwidth = 0.1, color="white", fill="#000099") +
  geom_vline(xintercept = quantile(perCIM_em[,1], 0.95), linetype="solid", color="red", size=0.7) +
  annotate("text", x = 6.5, y = 100, label = "Threshold = 4.52", size = 3.5)

#hk

perCIM_hk <- cim(cross = simMim, pheno.col = 8, window = 10, n.marcovar = 3, method = "hk", n.perm = 1000)
summary(perCIM_hk, alpha=.05)

ggplot(as.data.frame(perCIM_hk), aes(x=perCIM_hk[,1])) + 
  xlim(1,7) + xlab("Maximum LOD") + ylab("Frequency") +
  ggtitle("CIM - Halley-Knott") + theme(plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5, size=13)) +
  geom_histogram(binwidth = 0.1, color="white", fill="#009E73") +
  geom_vline(xintercept = quantile(perCIM_hk[,1], 0.95), linetype="solid", color="red", size=0.6) +
  annotate("text", x = 6, y = 150, label = "Threshold = 4.43", size = 3.5)

# Mapping

#em
CIMscan_em <- cim(cross = simMim, pheno.col = 8, window = 10, n.marcovar = 3, method = "em")
summary.CIMscan.em <- summary(CIMscan_em,perms = perCIM_em, alpha = 0.05, pvalues = T)

sumpos_cim <- data.frame(chr=CIMscan_em[,1], pos=cumsum(CIMscan_em[,2]))
a.cim <- as.vector(tapply(sumpos_cim$pos, sumpos_cim$chr, mean))
bpLabels <- paste(1:14)

ggplot(as.data.frame(CIMscan_em), aes(cumsum(CIMscan_em[,2]), CIMscan_em[,3])) + 
  geom_line() + geom_area(fill=alpha('blue',0.2)) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  geom_hline(yintercept = 4.52, linetype="solid", color="red", size=0.6) +
  ggtitle("CIM - em") + theme(plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5, size=13)) +
  scale_x_continuous(breaks=a.cim, labels=bpLabels) + xlab("Chromosome") + ylab("LOD")

#hk
CIMscan_hk <- cim(cross = simMim, pheno.col = 8, window = 10, n.marcovar = 3, method = "hk")
summary.CIMscan.hk <- summary(CIMscan_hk, perms = perCIM_hk, alpha = 0.05, pvalues = T)

ggplot(as.data.frame(CIMscan_hk), aes(cumsum(CIMscan_hk[,2]), CIMscan_hk[,3])) + 
  geom_line() + geom_area(fill=alpha('red',0.2)) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  geom_hline(yintercept = 4.43, linetype="solid", color="red", size=0.6) +
  ggtitle("CIM - Halley-Knott") + theme(plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5, size=13)) +
  scale_x_continuous(breaks=a, labels=bpLabels) + xlab("Chromosome") + ylab("LOD")

#plot em x hk
plot(CIMscan_em, CIMscan_hk, chr=c(1:14), col=c("blue", "orange"), incl.markers = T, show.marker.names = F, bandcol="gray70")
add.threshold(out=CIMscan_em, perm=perCIM_em,col="black",lwd=2)
add.threshold(out=CIMscan_hk, perm=perCIM_hk,col="red",lwd=2)

ggplot(NULL, aes(cumsum(pos), lod)) + 
  geom_line(data = as.data.frame(CIMscan_em), colour="blue", size=0.8) + ylim(0,9) +
  geom_step(data = as.data.frame(CIMscan_hk), colour="magenta1", size=0.8) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  scale_x_continuous(breaks=a, labels=bpLabels) + xlab("Chromosome") + ylab("LOD") +
  geom_hline(yintercept = 4.52, linetype="solid", color="red", size=0.6) +
  geom_hline(yintercept = 4.43, linetype="dashed", color="black", size=0.6) +
  ggtitle("CIM - em x hk") + theme(plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5, size=15))

# plot difference LOD between methods
plot(IMscan_em - IMscan_hk, col = c("blue"))

# rodar a função CIMscan_em mais uma vez se a diferença de lod desse gráfico der perto de 10
ggplot(as.data.frame(CIMscan_em[,3] - CIMscan_hk[,3]), aes(cumsum(CIMscan_hk[,2]), CIMscan_em[,3] - CIMscan_hk[,3])) + 
  geom_line() + geom_area(fill=alpha('red',0.2)) + #ylim(-0.3,0.9) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  ggtitle("Difference LOD in CIM: em - hk") + theme(plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5, size=13)) +
  scale_x_continuous(breaks=a, labels=bpLabels) + xlab("Chromosome") + ylab("LOD")

# plot principal QTLs
plot(CIMscan_em, CIMscan_hk, chr=c(6,13), col=c("blue", "magenta"), incl.markers = T, show.marker.names = F, bandcol="lightgray")

# tabela 1.2
tabela.cim = data.frame("Linkage Group"=c(summary.CIMscan.hk$chr[1],summary.CIMscan.hk$chr[2]),
                       "Position"=c(summary.CIMscan.hk$pos[1],summary.CIMscan.hk$pos[2]),
                       "LOD"=c(round(summary.CIMscan.hk$lod[1],2),round(summary.CIMscan.hk$lod[2],2)))
kable(tabela.cim)

###############################################################################################
###                                   COMPARING IM x CIM                                    ### 
###############################################################################################

# em
plot(IMscan_em, CIMscan_em, chr=c(1:14), col=c("blue", "orange"))
add.threshold(out=IMscan_em, perm=perIM_em,col="red",lwd=2)
add.threshold(out=CIMscan_em, perm=perCIM_em,col="green",lwd=2)

ggplot(NULL, aes(cumsum(pos), lod)) + 
  geom_line(data = as.data.frame(IMscan_em), colour="green", size=0.8) + ylim(0,9) +
  geom_step(data = as.data.frame(CIMscan_em), colour="blue", size=0.8) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  scale_x_continuous(breaks=a, labels=bpLabels) + xlab("Chromosome") + ylab("LOD") +
  geom_hline(yintercept = 4.52, linetype="solid", color="red", size=0.6) +
  geom_hline(yintercept = 3.86, linetype="dashed", color="black", size=0.6) +
  ggtitle("CIM x IM (em)") + theme(plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5, size=15))

# hk
plot(IMscan_hk, CIMscan_hk, chr=c(1:14), col=c("blue", "orange"))
add.threshold(out=IMscan_hk, perm=perIM_hk,col="red",lwd=2)
add.threshold(out=CIMscan_hk, perm=perCIM_hk,col="green",lwd=2)

ggplot(NULL, aes(cumsum(pos), lod)) + 
  geom_line(data = as.data.frame(IMscan_hk), colour="green1", size=0.8) + ylim(0,9) +
  geom_step(data = as.data.frame(CIMscan_hk), colour="blue", size=0.8) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  scale_x_continuous(breaks=a, labels=bpLabels) + xlab("Chromosome") + ylab("LOD") +
  geom_hline(yintercept = 4.43, linetype="solid", color="red", size=0.6) +
  geom_hline(yintercept = 3.83, linetype="dashed", color="black", size=0.6) +
  ggtitle("CIM x IM (hk)") + theme(plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5, size=15))

# dentro das regiões QTL (6 e 13)
plot(IMscan_em, CIMscan_em, chr=c(6), col=c("green", "blue"))
plot(IMscan_em, CIMscan_em, chr=c(13), col=c("green", "blue"))


###############################################################################################
###                           MULTIPLE INTERVAL MAPPING (MIM)                               ###
###############################################################################################

# Estimação das penalizações

tmpMim <- calc.genoprob(simMim, step=5, map.function = "kosambi") #objeto temporário
scannedTwo <- scantwo(tmpMim, pheno.col=8, method="hk", n.perm=1000, verbose = T)
penalties <- calc.penalties(scannedTwo,alpha=0.05)

# Determinando os QTL
qtl <- makeqtl(simMim, chr=c(6,13), pos=c(186,56), what="prob")

rqtl <- refineqtl(simMim, qtl=qtl, method="hk", pheno.col = 8)

# Função stepwiseqtl para determinação do modelo mais provavel e interações epistatica
stepQtl <- stepwiseqtl(simMim, max.qtl=5, method="hk", penalties = penalties,
                       verbose=T, pheno.col = 8, qtl = rqtl,
                       refine.locations = T)

# Marcadores significativos
summary(stepQtl)

# segundo a stepwise o modelo é:
formula <- y ~ Q1 + Q2 + Q3 + Q1:Q3

# Determinando os QTL
qtl2 <- makeqtl(simMim, chr=c(6,8,13), pos=c(186,87,56), what="prob")

# Refinando a posição dos QTL
rqtl2 <- refineqtl(simMim, qtl=qtl2, method="hk", pheno.col = 8)

out.fqr <- fitqtl(simMim, qtl=rqtl2, method="hk",pheno.col = 8,
                  formula = formula, get.ests = T)
summary(out.fqr) # obs: todos foram significativos

# análise de variancia
out.fqr$result.full

# efeitos
out.fqr$ests$ests

# Plot dos QTLs
plot(stepQtl)
plotModel(rqtl2, formula, chronly = F, cex.name = 1)

# Plot dos LOD
plotLodProfile(stepQtl)


###############################################################################################
###                                   COMPARING CIM x MIM                                    ### 
###############################################################################################

# Comparando CIM e MIM

# CIM: não considera o 8 como significativo (verificar pelo threshold)
plot(CIMscan_em, chr=c(6,8,13), col=c("blue"), incl.markers = T, show.marker.names = F, bandcol="lightgray")
add.threshold(out=CIMscan_em, perm=perCIM_em,col="red",lwd=2)

cxm <- data.frame("CIM"=c("Chr 6 - pos 185", "-", "Chr 13 - pos 58"),
                  "MIM"=c("Chr 6 - pos 186", "Chr 8 - pos 87", "Chr 13 - pos 56"))

#kable(cxm)

plot(CIMscan_em, chr=c(6,8,13), col=c("black"), ylim = c(0,20)) # CIM: apenas 6 e 13 signif
plotLodProfile(stepQtl, lty = c("solid", "solid", "solid"), 
               col = c("blue", "red", "green"), add = T, mtick = "line") # MIM: 6, 8 e 13