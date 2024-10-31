###Variance Partitioning Analysis

library("ecodist")
library("vegan")
spe<-read.csv("otu22.csv")
env<-read.csv("env.csv")

# 2.1.1 RDA-environmental variables
RDA=rda(decostand(spe, "hel"),env)
RDA
sumr=summary(RDA)

# 2.1.2 RDA-WC、pH and TC
## R Partial effects of RDA-WC, pH, and TC: variance explained solely by WC, pH, and TC individually
RDA2=rda(decostand(spe, "hel"),env[c(1,2,5)],env[c(3,4,6)]) 
RDA2
sumr2=summary(RDA2)
anova(RDA2) # Significance test.

#permutest(RDA2,permu=999) #permutest(RDA2, permu=999) # same function as anova()
RsquareAdj(RDA2) # Specify the adjustment of environmental factor correlation results. When there is more than one environmental factor in the classification, the R^2 results need to be adjusted

## RMarginal effects of RDA-WC, pH, and TC: variance explained solely by WC, pH, and TC + variance jointly explained by WC, pH, TC, and other environmental factors
RDA4=rda(decostand(spe, "hel"),env[c(1,2,5)]) 
RDA4
sumr4=summary(RDA4)
