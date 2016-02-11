library(cowplot)
library(dplyr)
library(data.table)
library(tidyr)
library(rtracklayer)
# par(mfrow=c(1,1))
par(mfrow = c(2, 2))
setwd("~/Documents/UCDavis/ECL243/git/ECL243/Assignment4/Maize_Admixture")
Maize.Admixture.2<-read.table("Maize.2.qopt")
Maize.Admixture.3<-read.table("Maize.3.qopt")
Maize.Admixture.4<-read.table("Maize.4.qopt")
Maize.Admixture.5<-read.table("Maize.5.qopt")

barplot(t(as.matrix(Maize.Admixture.2)), col=rainbow(4), xlab="Individuals", 
        ylab="Admixture Proportion", border=NA, main="K=2") 
barplot(t(as.matrix(Maize.Admixture.3)), col=rainbow(4), xlab="Individuals", 
        ylab="Admixture Proportion", border=NA, main="K=3") 
barplot(t(as.matrix(Maize.Admixture.4)), col=rainbow(4), xlab="Individuals", 
        ylab="Admixture Proportion", border=NA, main="K=4") 
barplot(t(as.matrix(Maize.Admixture.5)), col=rainbow(4), xlab="Individuals", 
        ylab="Admixture Proportion", border=NA, main="K=5") 
mtext("Maize", outer=TRUE,  cex=1, line=-2)
setwd("~/Documents/UCDavis/ECL243/git/ECL243/Assignment4/Teosinte_Admixture")
Teosinte.Admixture.2<-read.table("Teosinte.2.qopt")
Teosinte.Admixture.3<-read.table("Teosinte.3.qopt")
Teosinte.Admixture.4<-read.table("Teosinte.4.qopt")
Teosinte.Admixture.5<-read.table("Teosinte.5.qopt")
# par(mfrow=c(1,1))
par(mfrow = c(2, 2))
barplot(t(as.matrix(Teosinte.Admixture.2)), col=rainbow(4), xlab="Individuals", 
        ylab="Admixture Proportion", border=NA, main="K=2") 
barplot(t(as.matrix(Teosinte.Admixture.3)), col=rainbow(4), xlab="Individuals", 
        ylab="Admixture Proportion", border=NA, main="K=3") 
barplot(t(as.matrix(Teosinte.Admixture.4)), col=rainbow(4), xlab="Individuals", 
        ylab="Admixture Proportion", border=NA, main="K=4") 
barplot(t(as.matrix(Teosinte.Admixture.5)), col=rainbow(4), xlab="Individuals", 
        ylab="Admixture Proportion", border=NA, main="K=5") 
mtext("Teosinte", outer=TRUE,  cex=1, line=-2)
setwd("~/Documents/UCDavis/ECL243/git/ECL243/Assignment4/Mexicana_Admixture")
Mexicana.Admixture.2<-read.table("Mexicana.2.qopt")
Mexicana.Admixture.3<-read.table("Mexicana.3.qopt")
Mexicana.Admixture.4<-read.table("Mexicana.4.qopt")
Mexicana.Admixture.5<-read.table("Mexicana.5.qopt")
# par(mfrow=c(1,1))
par(mfrow = c(2, 2))
barplot(t(as.matrix(Mexicana.Admixture.2)), col=rainbow(4), xlab="Individuals", 
        ylab="Admixture Proportion", border=NA, main="K=2") 
barplot(t(as.matrix(Mexicana.Admixture.3)), col=rainbow(4), xlab="Individuals", 
        ylab="Admixture Proportion", border=NA, main="K=3") 
barplot(t(as.matrix(Mexicana.Admixture.4)), col=rainbow(4), xlab="Individuals", 
        ylab="Admixture Proportion", border=NA, main="K=4") 
barplot(t(as.matrix(Mexicana.Admixture.5)), col=rainbow(4), xlab="Individuals", 
        ylab="Admixture Proportion", border=NA, main="K=5") 
mtext("Mexicana", outer=TRUE,  cex=1, line=-2)
