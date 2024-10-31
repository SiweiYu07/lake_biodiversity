###Piecewise linear regresion

library(SiZer)
library(ggplot2)
seg<-read.csv("Robustness.csv")

#Exploring the relationship through smoothing fit reveals that there seems to be a threshold time
ggplot(seg,aes(Temperature,Robustness))+geom_point(size=10,color="#33CCCC")+geom_smooth(size=3,color="blue",fill="#FF9999")
  +xlab("Temperature")+ylab("Robustness")
+theme(axis.text.x=element_text(size=27),axis.text.y = element_text(size=27),axis.title.x=element_text(size=27),axis.title.y=element_text(size=27))
+theme(panel.border = element_rect(color = "black",fill = NA,linewidth = 1))
model <- piecewise.linear(x = seg$Temperature, y = seg$Robustness,CI = TRUE, bootstrap.samples = 100, sig.level = 0.05)

#It suggests that the relationship can potentially be described by piecewise linear regression, where the threshold time is the breakpoint
#Perform piecewise linear regression
#The breakpoint position is automatically determined, and a 95% confidence interval is estimated using 1000 bootstraps
model <- piecewise.linear(x = seg$Temperature, y = seg$Robustness, 
                          CI = TRUE, bootstrap.samples = 1000, sig.level = 0.05)
model

#Create a simple plot to view the piecewise linear regression graph
plot(model, xlab = 'Temperature', ylab = 'Robustness',pch=19,col="#33CCCC",cex=2,cex.axis=2,cex.lab=2,cex.sub=2,cex.main=2,lwd=4)
lines(model,col="blue",lwd=4)
#SSre：the sum of the squares of the distances of the points from the fit
#SStot：the sum of the squares of the distances of the points from a horizontal line through the mean of all Y values
fit1 <- predict(model, seg$Temperature)
SSre <- sum((seg$Robustness-fit1)^2)
SStot <- sum((seg$Robustness-mean(seg$Robustness))^2)
R2 <- round(1 - SSre/SStot, 3)
R2

library(segmented)

#fit a simple linear regression model
fit_lm <- lm(Robustness~Temperature, data = seg)
summary(fit_lm)

#Using the known linear regression model above, identify possible breakpoints with the segmented() function
#For details, see ?segmented. Two methods are available

#First method: specify the number of breakpoints using npsi
#For example, with 1 breakpoint, it will automatically search for the possible breakpoint position across the entire range
lm_seg1 <- segmented(fit_lm, seg.Z = ~TP, npsi = 2)
summary(lm_seg1)
plot(lm_seg1, xlab = 'Temperature', ylab = 'Robustness',col="blue",lwd=8)
points(Robustness~Temperature, data = seg,pch=19,col="#33CCCC",cex=2,cex.axis=2,cex.lab=2,cex.sub=2,cex.main=2,lwd=20)


#Second method: specify an approximate initial position for the breakpoint using psi
#For example, if x=1997 is a potential breakpoint position, setting it will automatically search for the optimal breakpoint position near x=1997
lm_seg2 <- segmented(fit_lm, seg.Z = ~TP, psi = 1997) 
summary(lm_seg2)

plot(lm_seg2, xlab = 'TP', ylab = 'Richness')

points(Richness~TP, data = seg)

#Note: The two methods in segmented() may yield different breakpoint positions or parameter estimates
#The breakpoint positions or parameter estimates obtained may also differ between methods from different R packages

