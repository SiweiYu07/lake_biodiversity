## alpha diversity
rm(list=ls())#clear Global Environment

#install packages
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggprism)
library(vegan)
library(picante)
library(dplyr)
library(RColorBrewer)

1.Filter data
otu<-read.csv("CHabundance.csv")
otu = read.table('CHabundance.txt', header=T, sep="\t", quote = "", row.names=1, comment.char="", stringsAsFactors = FALSE)

#calculate sum
colSums(otu)
otu_Flattening =as.data.frame(t(rrarefy(t(otu),min(colSums(otu)))))

#check out_Flattening
colSums(otu_Flattening)

#save the file and named as “otu_Results.csv”
write.table(otu_Flattening, file="choupingtable.csv",sep =",", quote=FALSE)

2.Alpha diversity (ACE, Chao,Shannon,Simpson)

##import data, column name was the OTU data, row name are the sample name
df <- read.table("CHabundance.txt",header = T, row.names = 1, check.names = F)

## select csv, delete species names of the first column 
df<-read.csv("CHabundance.csv")

#calculate shannon, simpson,richness with vegan package
Shannon <- diversity(df, index = "shannon", MARGIN = 2, base = exp(1))
Simpson <- diversity(df, index = "simpson", MARGIN = 2, base =  exp(1))
Richness <- specnumber(df, MARGIN = 2)

###tabulate diversity indices
index <- as.data.frame(cbind(Shannon, Simpson, Richness))

#transpose matrix
tdf <- t(df)
tdf<-ceiling(as.data.frame(t(df)))

#calculate obs，chao，ace indices
obs_chao_ace <- t(estimateR(tdf))
obs_chao_ace <- obs_chao_ace[rownames(index),]#unify the row name

#merge obs，chao，ace indices with above indices
index$Chao <- obs_chao_ace[,2]
index$Ace <- obs_chao_ace[,4]
index$obs <- obs_chao_ace[,1]

#calculate Pielou and coverage
index$Pielou <- Shannon / log(Richness, 2)
index$Goods_coverage <- 1 - colSums(df ==1) / colSums(df)

#export table
write.table(cbind(sample=c(rownames(index)),index),'CH_diversity.index.txt', row.names = F, sep = '\t', quote = F)



-----------------------------------------------------------------
