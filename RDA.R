### RDA analysis 

rm(list=ls())#clear Global Environment


#import packages
library(vegan)
library(ggplot2)

#import OTU table
df <- read.table("CHabundance.txt",sep="\t",header = T,row.names = 1,check.names = F)
head(df)

#import environment data
env <- read.table("CHenv.txt",sep="\t",header = T,row.names = 1,check.names = F)
head(env)

print(decorana(t(df)))
#Select the analysis method based on the Axis Lengths value of DCA1: if it is >4.0, choose CCA; if it is between 3.0 and 4.0, either RDA or CCA can be used; if it is <3.0, choose RDA analysis.

RDA <- rda(t(df),env,scale = T)
df_rda <- data.frame(RDA$CCA$u[,1:2],rownames(env))
colnames(df_rda)=c("RDA1","RDA2","samples")

#extract species scores
df_rda_score <- data.frame(RDA$CCA$v[,1:2])

#Calculate axis label data (= axis eigenvalue / sum of eigenvalues of all axes)
RDA1 =round(RDA$CCA$eig[1]/sum(RDA$CCA$eig)*100,2)
RDA2 =round(RDA$CCA$eig[2]/sum(RDA$CCA$eig)*100,2)

p1<-ggplot(data=df_rda,aes(x=RDA1,y=RDA2))+#Specify data, X-axis, Y-axis, and color
  theme_bw()+#Theme settings
  geom_point(size=4,shape=16,color = "#339933")+#Plot a scatter plot and set the size
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+#Dashed line in the plot
  geom_text(aes(label=samples, y=RDA2+0.03,x=RDA1+0.03,  vjust=0),size=4)+#Add labels to data points
  # guides(color=guide_legend(title=NULL))+#Remove legend title.
  labs(x=paste0("RDA1 (",RDA1,"%)"),
       y=paste0("RDA2 (",RDA2,"%)"))+#Change the x and y-axis titles to Contribution
  stat_ellipse(data=df_rda,
               level=0.95,
               linetype = 2,linewidth =0.8,
               show.legend = T)+
  scale_fill_manual(values = c("#339933"))+
  theme(axis.title.x=element_text(size=16),#Modify the X-axis title text.
        axis.title.y=element_text(size=16,angle=90),#Modify the y-axis title text.
        axis.text.y=element_text(size=14),#Modify the y-axis tick label text.
        axis.text.x=element_text(size=14),#Modify the X-axis tick label text.
        panel.grid=element_blank())#Hide grid lines.

p1

#Extract environmental factor scores
df_rda_env <- RDA$CCA$biplot[,1:2]
df_rda_env <- as.data.frame(df_rda_env)
head(df_rda_env)
# Add environmental factor data.
p2<-p1+
  geom_segment(data=df_rda_env,aes(x=0,y=0,xend=df_rda_env[,1],yend=df_rda_env[,2]),
               color="#3366FF",size=1,
               arrow=arrow(angle = 35,length=unit(0.4,"cm")))+
  geom_text(data=df_rda_env,aes(x=df_rda_env[,1],y=df_rda_env[,2],
                                label=rownames(df_rda_env)),size=7,
            color="#3366FF", fontface ="bold",
            hjust="inward",
            vjust=0.5*(1-sign(df_rda_env[,2])))+
  theme(legend.position = "top")
p2

p3<-p2+
  geom_point(data=df_rda_score,size=3,shape=2,color="Red")+
p3




