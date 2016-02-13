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

# 
barplot(t(as.matrix(all.Admixture.2)), col=rainbow(4), xlab="Individuals", 
        ylab="Admixture Proportion", border=NA, main="K=2") 
barplot(t(as.matrix(all.Admixture.3)), col=rainbow(4), xlab="Individuals", 
        ylab="Admixture Proportion", border=NA, main="K=3") 
barplot(t(as.matrix(all.Admixture.4)), col=rainbow(4), xlab="Individuals", 
        ylab="Admixture Proportion", border=NA, main="K=4") 
barplot(t(as.matrix(all.Admixture.5)), col=rainbow(4), xlab="Individuals", 
        ylab="Admixture Proportion", border=NA, main="K=5") 

