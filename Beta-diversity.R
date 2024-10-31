######Beta diversity (dissimilarity, turnover, nestedness)

library(betapart)
library(vegan)

1. Abundance-based pair-wise dissmilarities（Bray-curties）
abund_phyto<-read.table("CH_abundn.txt",header=T,row.names = 1)
## .txt file first row name is otu name, column name is sample name
abund_phyto<-read.csv("CH_abund.csv")
abund_total<-beta.multi.abund(abund_phyto,index.family = "bray")
abund_total
abund_phyto_detail<-beta.pair.abund(abund_phyto,index.family = "bray")
abund_phyto_detail
write.table(abund_phyto_detail,'abund_CH_result.csv')

library(dplyr)
incid_phyto %>% mutate(across(where(is.numeric), ~ +as.logical(.x)))
incid_phyto %>% mutate_if(is.numeric, ~1 * (. != 0))
write.table(incid_phyt,'incid_CH.csv')

2. Incidence-based pair-wise dissimilarities (jaccard)
incid_phyto<-read.table("CH_incid.txt",header=T,row.names = 1)
## .txt file first row name is otu name, first column name is sample name 
incid_phyto<-read.csv("CH_incid.csv")
incid_total<-beta.multi(incid_phyto, index.family = "jaccard")
incid_total
incid_phyto<-beta.pair(incid_phyto,index.family = "jac")
incid_phyto
write.table(incid_phyto_detail,'incid_CH.csv')


