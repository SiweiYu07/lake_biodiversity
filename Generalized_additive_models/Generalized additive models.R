###Generalized additive models

library("mgcv")
#> Loading required package: nlme
#> This is mgcv 1.8-24. For overview type 'help("mgcv-package")'.
library("scam")
library("ggplot2")
library("cowplot")
library("tidyr")
## gratia is not on CRAN, install from github:
install.packages("devtools")
devtools::install_github("gavinsimpson/gratia")
library("gratia") # need to change the name of the package
## Default ggplot theme
theme_set(theme_bw())
## source Small Water data
TNTurnover <- read.csv("TNTurnover.csv",sep = ',',header = TRUE)
head(TNTurnover)
## load braya so data set
TPTurnover <- read.csv("TPTurnover.csv",sep = ',',header = TRUE)
## clean up variable names
names(TNTurnover) <- c("Depth", "DepthUpper", "DepthLower", "Year", "YearYoung",
                  "YearOld", "Turnover")
## add a variable for the amount of time per sediment sample
braya <- transform(TNTurnover, sampleInterval = YearYoung - YearOld)
head(braya)

#GAM线性拟合：

small<-read.csv("TNTurnover.csv")
Sim<-expression(Sim)
small_plt <- ggplot(small, aes(x = TP, y = Sim)) +
   geom_point(size=10) +
   labs(y = Sim, x = "TP")
plot(small_plt, ncol = 1, labels = "auto", align = "hv",     axis = "lr")
m <- gam(Sim ~ s(TP, k = 15), data = small, method = "REML")
mod <- gamm(Sim ~ s(TP, k = 15), data = small,  correlation = corCAR1(form = ~ Depth), method = "REML")
smallPhi <- intervals(mod$lme, which = "var-cov")$corStruct
summary(mod$gam)
S <- seq(0, 50, length = 150)
car1 <- setNames(as.data.frame(t(outer(smallPhi, S, FUN = `^`)[1, , ])), c("Lower","Correlation","Upper"))
car1 <- transform(car1, S = S)
car1Plt <- ggplot(car1, aes(x = S, y = Correlation)) +geom_ribbon(aes(ymax = Upper, ymin = Lower),fill = "black", alpha = 0.2) +geom_line() 

newYear <- with(small, data.frame(TP = seq(min(TP), max(TP),length.out = 200)))
newYear <- cbind(newYear,data.frame(predict(mod$gam, newYear, se.fit = TRUE)))
crit.t <- qt(0.975, df = df.residual(mod$gam))
newYear <- transform(newYear,upper = fit + (crit.t * se.fit),lower = fit - (crit.t * se.fit))
sfittt <- ggplot(newYear, aes(x = TP, y = fit))+
  geom_ribbon(aes(ymin = lower, ymax = upper, x = TP), alpha = 0.4,inherit.aes = FALSE, fill = "#66CCCC") +
  geom_point(data = small, size=7,color="#3366CC", mapping = aes(x = TP, y = Sim),inherit.aes = FALSE) +
  geom_line(linewidth=2,color="#FF3366") +
  labs(y = Sim, x = "TP")+
  theme(axis.text.x=element_text(vjust=1, hjust=1,size=15),
   axis.text.y=element_text(size=15))
sfittt

library("patchwork")
sfitt+sfittt
