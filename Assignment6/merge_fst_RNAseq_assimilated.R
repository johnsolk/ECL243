# For enrichment analysis:
# make bed file 
# https://www.biostars.org/p/12107/
# https://www.biostars.org/p/12283/
# with chr, winstart, winend, fst_teo_mz
# use bedtools to compare with gff file
# what is list of genes in high Fst windows?
# what % of assimilation genes in high Fst windows?

setwd("~/Documents/UCDavis/ECL243/git/ECL243/Assignment6")
# SLR (column 28) = number of SNP in land-raised maize
# filter greater than 500 SLR

# commands:
# awk -F '\t' '$28 > 500' Hufford_et_al._2012_10kb_statistics.txt > Hufford_et_al._2012_10kb_statistics_filtered.txt
# cut -f 1,9,10,15 Hufford_et_al._2012_10kb_statistics_filtered.txt > teo_maize_fst.bed
# sed '/NA/d' teo_maize_fst.bed > teo_maize_fst_fixed.bed
# sed '/chr/d' teo_maize_fst_fixed.bed > teo_maize_fst_fixed_noheader.bed
# sed '/UNMAPPED/d' ZmB73_4a.53_FGS_mssgd_fixed.gff3 > ZmB73_4a.53_FGS_mssgd_fixed2.gff3
# bedtools intersect -a teo_maize_fst_fixed_noheader.bed -b ZmB73_4a.53_FGS_mssgd_fixed_unmapped.gff3 -wa -wb > regions

fst=read.table("regions",sep="\t")
assim<-read.csv("Teosinte_subset_notsigmaize_padj0.05.csv")

# need to merge fst with assim based on ID,
# problem is ENSEMBL that "ID=ENSEMBL" is within fst$V13

fst_ID<-strsplit(as.character(fst$V13),";")
assim_ID<-assim$ensembl_id_fixed

length(assim_ID)
# [1] 2110
length(fst_ID)
# [1] 3820

# From hints, use sample to test for enrichment:
# sample takes a sample of the specified size from the elements 
# of x using either with or without replacement.
