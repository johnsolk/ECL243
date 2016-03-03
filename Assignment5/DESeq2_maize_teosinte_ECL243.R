library(DESeq2)
library("lattice")
library(biomaRt)
library("RColorBrewer")
library("gplots")
library(dplyr)
library(magrittr)
 
source('~/Documents/scripts/plotPCAWithSampleNames.R')
setwd("~/Documents/UCDavis/ECL243/git/ECL243/Assignment5")

teo<-read.table("teosintereadcount-2.csv",header=T,row.names=1,sep=",")
maize<-read.table("Maizereadcount.csv",header=T,row.names=1,sep=",")

cols <-data.frame(genotype=c("RIMMA1","RIMMA1","RIMMA19","RIMMA19","RIMMA140","RIMMA140","RIMMA809",
"RIMMA809","RIMMA1","RIMMA1","RIMMA19","RIMMA19","RIMMA140","RIMMA140","RIMMA809","RIMMA809"), 
sample=c("RIMMA1-1","RIMMA1-1","RIMMA19-1","RIMMA19-1",
"RIMMA140-1","RIMMA140-1","RIMMA809-1","RIMMA809-1","RIMMA1-2","RIMMA1-2",
"RIMMA19-2","RIMMA19-2","RIMMA140-2","RIMMA140-2","RIMMA809-2","RIMMA809-2"), 
condition=as.factor(c(rep("265ppm",8), rep("400ppm",8))))
designModel = formula(~ genotype + condition)


dds.maize <- DESeqDataSetFromMatrix(as.matrix(maize), colData=cols, design = designModel)
ddsColl.maize<-collapseReplicates(dds.maize, groupby= dds.maize$sample)
#ddsColl.maize<- DESeq(dds.maize)
ddsColl.maize<- DESeq(ddsColl.maize)
#ddsColl<- DESeq(ddsColl,betaPrior=FALSE)

teosinte.data<-read.csv("teosintereadcount-2.csv",header=T,row.names=1,sep=",")
# 400ppm is control (modern)
# 265ppm (ancient)
cols<-data.frame(genotype=c("Pop1","Pop1","Pop1","Pop2","Pop2","Pop2","Pop3","Pop3",
"Pop3","Pop4","Pop4","Pop4","Pop1","Pop1","Pop1","Pop2","Pop2","Pop2","Pop3","Pop3",
"Pop3","Pop4","Pop4"),condition=as.factor(c(rep("265ppm",12),rep("400ppm",11))))
designModel = formula(~ genotype + condition)
dds.teo <- DESeqDataSetFromMatrix(as.matrix(teosinte.data), colData=cols, design =designModel)
ddsColl.teo<- DESeq(dds.teo)

# log2 transformation for PCA plot
log_maize<-rlog(ddsColl.maize)
log_teo<-rlog(ddsColl.teo)

teo.names<-c("Pop1-1-265","Pop1-2-265","Pop1-3-265","Pop2-1-265","Pop2-2-265","Pop2-3-265","Pop3-1-265","Pop3-2-265",
                     "Pop3-3-265","Pop4-1-265","Pop4-2-265","Pop4-3-265","Pop1-1-400","Pop1-2-400","Pop1-3-400","Pop2-1-400",
                     "Pop2-2-400","Pop2-3-400","Pop3-1-400","Pop3-2-400","Pop3-3-400","Pop4-1-400","Pop4-2-400")

colnames(log_teo)<-teo.names  
  
plotPCAWithSampleNames(log_maize, intgroup="condition", ntop=40000)
plotPCAWithSampleNames(log_teo, intgroup="condition", ntop=40000)

res.teosinte<-results(ddsColl.teo,contrast=c("condition","265ppm","400ppm"))
res.maize<-results(ddsColl.maize,contrast=c("condition","265ppm","400ppm"))

resultsNames(res.teosinte)
resultsNames(res.maize)
#plotMA(res2,alpha=0.01,main="Teosinte, 265 vs. 400",ylim=c(-6,7))
plot(log2(res.teosinte$baseMean), res.teosinte$log2FoldChange, 
     col=ifelse(res.teosinte$padj < 0.05, "red","gray67"),
     main="Teosinte (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")


plot(log2(res.maize$baseMean), res.maize$log2FoldChange, 
     col=ifelse(res.maize$padj < 0.05, "red","gray67"),
     main="Maize (padj<0.05)",xlim=c(0,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")

###
norm_counts.teo<-counts(ddsColl.teo,normalized=TRUE)
norm_counts.maize<-counts(ddsColl.maize,normalized=TRUE)
###
# Teosinte
# get gene names from Ensembl gene ID
#mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org")
mart = useMart(biomart = "plants_mart", host="plants.ensembl.org")
#listDatasets(mart)
ensembl=useMart("plants_mart",dataset="zmays_eg_gene", host="plants.ensembl.org")
# old command:
#ensembl=useMart("ensembl")
ensembl = useDataset("zmays_eg_gene",mart=ensembl)
data_table<-norm_counts.teo
ensembl_id<-rownames(data_table)
ensembl_id_fixed<-gsub("^gene:","",ensembl_id)
#listAttributes(ensembl)
query<-getBM(attributes=c('ensembl_gene_id','description','gene_biotype'), filters = 'ensembl_gene_id', values = ensembl_id_fixed, mart=ensembl)
#query<-getBM(attributes=c('description','ensembl_gene_id','start_position'),filters=c('start'),values=c(130000),mart=ensembl)
new_res<-as.data.frame(data_table)
data_table<-cbind(new_res,ensembl_id_fixed)
col.names<-c("ensembl_id_fixed",'description',"gene_biotype")
colnames(query)<-col.names
merge_biomart_res_counts <- merge(data_table,query,by="ensembl_id_fixed")
head(merge_biomart_res_counts)
temp_data_merged_counts<-merge_biomart_res_counts
###
# merge gene names with Teosinte expression data
head(res.teosinte)
ensembl_id<-rownames(res.teosinte)
ensembl_id_fixed<-gsub("^gene:","",ensembl_id)
new_res<-as.data.frame(res.teosinte)
new_res<-cbind(new_res,ensembl_id_fixed)
merge_biomart_res <- merge(new_res,query,by="ensembl_id_fixed")
head(merge_biomart_res)
merge_biomart_res_all<-merge(merge_biomart_res,temp_data_merged_counts,by="ensembl_id_fixed")
merge_biomart_res_all<-merge_biomart_res_all[order(merge_biomart_res_all$padj),]
merge_biomart_res_all<-subset(merge_biomart_res_all,merge_biomart_res_all$padj!="NA")
head(merge_biomart_res_all)
#write.csv(merge_biomart_res_all,file="Teosinte_all.csv")
res_merged_cutoff<-subset(merge_biomart_res_all,merge_biomart_res_all$padj<0.05)
dim(res_merged_cutoff)
#write.csv(res_merged_cutoff,file="Teosinte_padj0.05.csv")
up_1.5FC<-subset(res_merged_cutoff,res_merged_cutoff$log2FoldChange>1.5)
down_1.5FC<-subset(res_merged_cutoff,res_merged_cutoff$log2FoldChange< -1.5)
dim(up_1.5FC)
dim(down_1.5FC)
up_down_1.5FC<-subset(res_merged_cutoff,res_merged_cutoff$log2FoldChange>1.5 |res_merged_cutoff$log2FoldChange< -1.5)
dim(up_down_1.5FC)
head(up_down_1.5FC)
colnames(up_down_1.5FC)
#write.csv(up_down_1.5FC,file="Teosinte_padj0.05_log2FC1.5.csv")

# heatmap Teosinte
head(up_down_1.5FC)
colnames(up_down_1.5FC)
d.1.5FC <- as.matrix(up_down_1.5FC[,c(10:17)])
rownames(d.1.5FC) <- up_down_1.5FC[,1]
d.1.5FC<-na.omit(d.1.5FC)
#colnames(d.1.5FC)<-teo.names
hr <- hclust(as.dist(1-cor(t(d.1.5FC), method="pearson")), method="complete")
mycl <- cutree(hr, h=max(hr$height/1.5))
clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
myheatcol <- greenred(75)
#png("Teosinte_heatmap.png", width = 7*300,height = 7*300,res = 1200,pointsize = 2) 
heatmap.2(d.1.5FC, main="Teosinte, padj<0.05, log2FC +-1.5", 
          Rowv=as.dendrogram(hr),
          cexRow=0.15,cexCol=0.5,srtCol= 90,
          adjCol = c(NA,0),offsetCol=1, 
          Colv=NA, dendrogram="row", 
          scale="row", col=myheatcol, 
          density.info="none", 
          trace="none")
#dev.off()



###
# Maize
# get gene names from Ensembl gene ID
data_table<-norm_counts.maize
ensembl_id<-rownames(data_table)
ensembl_id_fixed<-gsub("^gene:","",ensembl_id)
#listAttributes(ensembl)
query<-getBM(attributes=c('ensembl_gene_id','description','gene_biotype'), filters = 'ensembl_gene_id', values = ensembl_id_fixed, mart=ensembl)
#query<-getBM(attributes=c('description','ensembl_gene_id','start_position'),filters=c('start'),values=c(130000),mart=ensembl)
new_res<-as.data.frame(data_table)
data_table<-cbind(new_res,ensembl_id_fixed)
col.names<-c("ensembl_id_fixed",'description',"gene_biotype")
colnames(query)<-col.names
merge_biomart_res_counts <- merge(data_table,query,by="ensembl_id_fixed")
head(merge_biomart_res_counts)
temp_data_merged_counts<-merge_biomart_res_counts
###
# merge gene names with Maize expression data
head(res.maize)
ensembl_id<-rownames(res.maize)
ensembl_id_fixed<-gsub("^gene:","",ensembl_id)
new_res<-as.data.frame(res.maize)
new_res<-cbind(new_res,ensembl_id_fixed)
merge_biomart_res <- merge(new_res,query,by="ensembl_id_fixed")
head(merge_biomart_res)
merge_biomart_res_all<-merge(merge_biomart_res,temp_data_merged_counts,by="ensembl_id_fixed")
merge_biomart_res_all<-merge_biomart_res_all[order(merge_biomart_res_all$padj),]
merge_biomart_res_all<-subset(merge_biomart_res_all,merge_biomart_res_all$padj!="NA")
head(merge_biomart_res_all)
#write.csv(merge_biomart_res_all,file="Maize_all.csv")
res_merged_cutoff<-subset(merge_biomart_res_all,merge_biomart_res_all$padj<0.05)
dim(res_merged_cutoff)
#write.csv(res_merged_cutoff,file="Maize_padj0.05.csv")
up_1.5FC<-subset(res_merged_cutoff,res_merged_cutoff$log2FoldChange>1.5)
down_1.5FC<-subset(res_merged_cutoff,res_merged_cutoff$log2FoldChange< -1.5)
dim(up_1.5FC)
dim(down_1.5FC)
up_down_1.5FC<-subset(res_merged_cutoff,res_merged_cutoff$log2FoldChange>1.5 |res_merged_cutoff$log2FoldChange< -1.5)
dim(up_down_1.5FC)
head(up_down_1.5FC)
colnames(up_down_1.5FC)
#write.csv(up_down_1.5FC,file="Maize_padj0.05_log2FC1.5.csv")

# heatmap Maize
head(up_down_1.5FC)
colnames(up_down_1.5FC)
d.1.5FC <- as.matrix(up_down_1.5FC[,c(10:32)])
rownames(d.1.5FC) <- up_down_1.5FC[,1]
d<-na.omit(d.1.5FC)
hr <- hclust(as.dist(1-cor(t(d.1.5FC), method="pearson")), method="complete")
mycl <- cutree(hr, h=max(hr$height/1.5))
clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
myheatcol <- greenred(75)
#png("Teosinte_heatmap.png", width = 7*300,height = 7*300,res = 1200,pointsize = 2) 
heatmap.2(d.1.5FC, main="Maize, padj<0.05, log2FC +-1.5", 
          Rowv=as.dendrogram(hr),
          cexRow=0.15,cexCol=0.5,srtCol= 90,
          adjCol = c(NA,0),offsetCol=1, 
          Colv=NA, dendrogram="row", 
          scale="row", col=myheatcol, 
          density.info="none", 
          trace="none")
#dev.off()

