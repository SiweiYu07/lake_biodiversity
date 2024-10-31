Triangle Correlation with significance 

library("pheatmap")
library("corrplot")
getwd()

#import data
al<-read.csv("CHEF.csv")
head(al)

tdc1 = cor (al, method="spearman")
tdc1

#choose color
addcol <- colorRampPalette(c("CornflowerBlue","White","tomato"))
testRes = cor.mtest(al, menthod = "spearman", conf.level = 0.95)
testRes$uppCI

corrplot(tdc1, method = "circle", 
         type="upper",
         col = addcol(100), 
         tl.col = "black", 
         tl.cex = 1, 
         tl.srt = 90,
         #tl.pos = "lt",
         cl.pos = "r",
         p.mat = testRes$p, 
         diag = F,
         sig.level = c(0.001, 0.01, 0.05), 
         pch.cex = 1,
         insig = 'label_sig',
         family="serif")

