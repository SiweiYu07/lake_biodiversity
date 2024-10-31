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
setwd("/Users/siweiyu/SW/Alpha_diversity")
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
write.table(cbind(sample=c(rownames(index)),index),'LG_protozoa.index.txt', row.names = F, sep = '\t', quote = F)



-----------------------------------------------------------------
3、差异性计算及绘图1）读入数据及分组文件#读入文件
index <- read.delim('diversity.index_pd.txt', header = T, row.names = 1)
##figure:take shannon for example
index$samples <- rownames(index)#将样本名写到文件中
#读入分组文件
groups <- read.delim('group.txt',header = T, stringsAsFactors = F)
colnames(groups)[1:2] <- c('samples','group')#改列名
#合并分组信息与多样性指数
df2 <- merge(index,groups,by = 'samples')
2）绘图--Shannon#Shannon
p1 <- ggplot(df2,aes(x=group,y=Shannon))+#指定数据
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+#添加误差线,注意位置，放到最后则这条线不会被箱体覆盖
  geom_boxplot(aes(fill=group), #绘制箱线图函数
               outlier.colour="white",size=0.8)+#异常点去除
  theme(panel.background =element_blank(), #背景
        axis.line=element_line(),#坐标轴的线设为显示
        plot.title = element_text(size=14))+#图例位置
  # scale_fill_manual(values=c("#ffc000","#a68dc8","blue"))+#指定颜色
  geom_jitter(width = 0.2)+#添加抖动点
  geom_signif(comparisons = list(c("A","B"),
                                 c("A","C"),
                                 c("B","C")),# 设置需要比较的组
              map_signif_level = T, #是否使用星号显示
              test = t.test, ##计算方法
              y_position = c(3,3.5,3.25),#图中横线位置 设置
              tip_length = c(c(0,0),
                             c(0,0),
                             c(0,0)),#横线下方的竖线设置
              size=0.8,color="black")+
  theme_prism(palette = "candy_bright",
              base_fontface = "plain", # 字体样式，可选 bold, plain, italic
              base_family = "serif", # 字体格式，可选 serif, sans, mono, Arial等
              base_size = 16,  # 图形的字体大小
              base_line_size = 0.8, # 坐标轴的粗细
              axis_text_angle = 45)+ # 可选值有 0，45，90，270
  scale_fill_prism(palette = "candy_bright")+
  theme(legend.position = 'none')#去除图例
p1
<img src="https://picx.zhimg.com/50/v2-a13c2bee7e69f87bdae65f21ec8d5b51_720w.jpg?source=1940ef5c" data-caption="" data-size="normal" data-rawwidth="640" data-rawheight="439" data-original-token="v2-a683f2a5d5339b32c53cccb70f839c4f" class="origin_image zh-lightbox-thumb" width="640" data-original="https://pic1.zhimg.com/v2-a13c2bee7e69f87bdae65f21ec8d5b51_r.jpg?source=1940ef5c"/>--Simpson#Simpson
p2 <- ggplot(df2,aes(x=group,y=Simpson))+#指定数据
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+#添加误差线,注意位置，放到最后则这条线不会被箱体覆盖
  geom_boxplot(aes(fill=group), #绘制箱线图函数
               outlier.colour="white",size=0.8)+#异常点去除
  theme(panel.background =element_blank(), #背景
        axis.line=element_line(),#坐标轴的线设为显示
        plot.title = element_text(size=14))+#图例位置
  # scale_fill_manual(values=c("#ffc000","#a68dc8","blue"))+#指定颜色
  geom_jitter(width = 0.2)+#添加抖动点
  geom_signif(comparisons = list(c("A","B"),
                                 c("A","C"),
                                 c("B","C")),# 设置需要比较的组
              map_signif_level = T, #是否使用星号显示
              test = t.test, ##计算方法
              y_position = c(1,1.1,1.05),#图中横线位置 设置
              tip_length = c(c(0,0),
                             c(0,0),
                             c(0,0)),#横线下方的竖线设置
              size=0.8,color="black")+
  theme_prism(palette = "candy_bright",
              base_fontface = "plain", # 字体样式，可选 bold, plain, italic
              base_family = "serif", # 字体格式，可选 serif, sans, mono, Arial等
              base_size = 16,  # 图形的字体大小
              base_line_size = 0.8, # 坐标轴的粗细
              axis_text_angle = 45)+ # 可选值有 0，45，90，270
  scale_fill_prism(palette = "candy_bright")+
  theme(legend.position = 'none')#去除图例
p2
<img src="https://pica.zhimg.com/50/v2-05af65bcafb816aa8b690e7ca8440801_720w.jpg?source=1940ef5c" data-caption="" data-size="normal" data-rawwidth="640" data-rawheight="427" data-original-token="v2-a18a7f9774069677f05b965d8c2a51f6" class="origin_image zh-lightbox-thumb" width="640" data-original="https://picx.zhimg.com/v2-05af65bcafb816aa8b690e7ca8440801_r.jpg?source=1940ef5c"/>--Ace#Ace
p3 <- ggplot(df2,aes(x=group,y=Ace))+#指定数据
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+#添加误差线,注意位置，放到最后则这条线不会被箱体覆盖
  geom_boxplot(aes(fill=group), #绘制箱线图函数
               outlier.colour="white",size=0.8)+#异常点去除
  theme(panel.background =element_blank(), #背景
        axis.line=element_line(),#坐标轴的线设为显示
        plot.title = element_text(size=14))+#图例位置
  # scale_fill_manual(values=c("#ffc000","#a68dc8","blue"))+#指定颜色
  geom_jitter(width = 0.2)+#添加抖动点
  geom_signif(comparisons = list(c("A","B"),
                                 c("A","C"),
                                 c("B","C")),# 设置需要比较的组
              map_signif_level = T, #是否使用星号显示
              test = t.test, ##计算方法
              y_position = c(55,65,60),#图中横线位置 设置
              tip_length = c(c(0,0),
                             c(0,0),
                             c(0,0)),#横线下方的竖线设置
              size=0.8,color="black")+
  theme_prism(palette = "candy_bright",
              base_fontface = "plain", # 字体样式，可选 bold, plain, italic
              base_family = "serif", # 字体格式，可选 serif, sans, mono, Arial等
              base_size = 16,  # 图形的字体大小
              base_line_size = 0.8, # 坐标轴的粗细
              axis_text_angle = 45)+ # 可选值有 0，45，90，270
  scale_fill_prism(palette = "candy_bright")+
  theme(legend.position = 'none')#去除图例
p3
<img src="https://picx.zhimg.com/50/v2-8189faa4ca9f82baca30b1eb6516f79d_720w.jpg?source=1940ef5c" data-caption="" data-size="normal" data-rawwidth="640" data-rawheight="439" data-original-token="v2-0fdf014e58d324550f7d2b32d993b407" class="origin_image zh-lightbox-thumb" width="640" data-original="https://picx.zhimg.com/v2-8189faa4ca9f82baca30b1eb6516f79d_r.jpg?source=1940ef5c"/>--Chao#Chao
p4 <- ggplot(df2,aes(x=group,y=Chao))+#指定数据
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+#添加误差线,注意位置，放到最后则这条线不会被箱体覆盖
  geom_boxplot(aes(fill=group), #绘制箱线图函数
               outlier.colour="white",size=0.8)+#异常点去除
  theme(panel.background =element_blank(), #背景
        axis.line=element_line(),#坐标轴的线设为显示
        plot.title = element_text(size=14))+#图例位置
  # scale_fill_manual(values=c("#ffc000","#a68dc8","blue"))+#指定颜色
  geom_jitter(width = 0.2)+#添加抖动点
  geom_signif(comparisons = list(c("A","B"),
                                 c("A","C"),
                                 c("B","C")),# 设置需要比较的组
              map_signif_level = T, #是否使用星号显示
              test = t.test, ##计算方法
              y_position = c(55,65,60),#图中横线位置 设置
              tip_length = c(c(0,0),
                             c(0,0),
                             c(0,0)),#横线下方的竖线设置
              size=0.8,color="black")+
  theme_prism(palette = "candy_bright",
              base_fontface = "plain", # 字体样式，可选 bold, plain, italic
              base_family = "serif", # 字体格式，可选 serif, sans, mono, Arial等
              base_size = 16,  # 图形的字体大小
              base_line_size = 0.8, # 坐标轴的粗细
              axis_text_angle = 45)+ # 可选值有 0，45，90，270
  scale_fill_prism(palette = "candy_bright")+
  theme(legend.position = 'none')#去除图例
p4
<img src="https://pic1.zhimg.com/50/v2-4e0e0ea5c258811ae43f166a01b3b676_720w.jpg?source=1940ef5c" data-caption="" data-size="normal" data-rawwidth="640" data-rawheight="439" data-original-token="v2-4de6e8b655de5cbbaedf39b314b857da" class="origin_image zh-lightbox-thumb" width="640" data-original="https://pic1.zhimg.com/v2-4e0e0ea5c258811ae43f166a01b3b676_r.jpg?source=1940ef5c"/>4、拼图library("gridExtra")
library("cowplot")
plot_grid(p1,p2,p3,p4, labels=c('A','B','C','D'), ncol=2, nrow=2)#拼图及标注
ggsave('xxx.pdf',width=12,height = 4)
<img src="https://picx.zhimg.com/50/v2-d03b5760eefd893cc4360e5df0977944_720w.jpg?source=1940ef5c" data-caption="" data-size="normal" data-rawwidth="640" data-rawheight="397" data-original-token="v2-c91f116c143681482ad246b142497ce2" class="origin_image zh-lightbox-thumb" width="640" data-original="https://pic1.zhimg.com/v2-d03b5760eefd893cc4360e5df0977944_r.jpg?source=1940ef5c"/>