library(cowplot)
library(dplyr)
library(data.table)
library(tidyr)
library(rtracklayer)
# par(mfrow=c(1,1))
par(mfrow = c(2, 2))
setwd("~/Documents/UCDavis/ECL243/git/ECL243/Assignment4/all_Admixture")
all.Admixture.2<-read.table("all.2.qopt")
all.Admixture.3<-read.table("all.3.qopt")
all.Admixture.4<-read.table("all.4.qopt")
all.Admixture.5<-read.table("all.5.qopt")

x.lab=c("Maize_BKN009","Maize_BKN011","Maize_BKN014","Maize_BKN015",
        "Maize_BKN019","Maize_BKN022","Maize_BKN025","Maize_BKN026","Maize_BKN027","Maize_BKN033",
        "Maize_BKN035","Teosinte_TIL01","Teosinte_TIL03","Teosinte_TIL04",
        "Teosinte_TIL07","Teosinte_TIL09","Teosinte_TIL11",
        "Teosinte_TIL15","Teosinte_TIL16","Teosinte_TIL17",
        "Mexicana_TIL08","Mexicana_TIL25")
# 
barplot(t(as.matrix(all.Admixture.2)), col=rainbow(4), ylab="Admixture Proportion",
         names.arg=x.lab, las=2, cex.names = 0.65, border=NA, main="K=2") 
barplot(t(as.matrix(all.Admixture.3)), col=rainbow(4), ylab="Admixture Proportion", 
        names.arg=x.lab, las=2, cex.names = 0.65, border=NA, main="K=3") 
barplot(t(as.matrix(all.Admixture.4)), col=rainbow(4), ylab="Admixture Proportion", 
         names.arg=x.lab, cex.names = 0.65, las=2, border=NA, main="K=4") 
barplot(t(as.matrix(all.Admixture.5)), col=rainbow(4), ylab="Admixture Proportion", 
        names.arg=x.lab, las=2, cex.names = 0.65,border=NA, main="K=5") 

