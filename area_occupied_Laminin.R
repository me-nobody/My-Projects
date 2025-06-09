library(tidyverse)
library(ggplot2)
library(corrplot)
library(ggplot2)
library(uwot)
library(scales)
library(FactoMineR)
library(factoextra)
library(fpc)
# add t-test to plot
library(ggpubr)
library(rstatix)
library(ggprism)
library(ggsignif)

# first checking how many of the 86 variables correlate with each other
# removing the following features
#  - label,slice,moments_normalized.0.0,moments_normalized.0.1,moments_normalized.1.0 : non-numerical data
# !(names(df) %in% cols_to_drop)
cols_to_drop <- c('label','slice','moments_normalized.0.0','moments_normalized.0.1','moments_normalized.1.0')
infil_roi1_area_laminin<-infil_roi1_area_laminin[,!(names(infil_roi1_area_laminin) %in% cols_to_drop)]
infil_roi2_area_laminin<-infil_roi2_area_laminin[,!(names(infil_roi2_area_laminin) %in% cols_to_drop)]
infil_roi3_area_laminin<-infil_roi3_area_laminin[,!(names(infil_roi3_area_laminin) %in% cols_to_drop)]
infil_roi4_area_laminin<-infil_roi4_area_laminin[,!(names(infil_roi4_area_laminin) %in% cols_to_drop)]
infil_roi5_area_laminin<-infil_roi5_area_laminin[,!(names(infil_roi5_area_laminin) %in% cols_to_drop)]
tumor_roi1_area_laminin<-tumor_roi1_area_laminin[,!(names(tumor_roi1_area_laminin) %in% cols_to_drop)]
tumor_roi2_area_laminin<-tumor_roi2_area_laminin[,!(names(tumor_roi2_area_laminin) %in% cols_to_drop)]
tumor_roi3_area_laminin<-tumor_roi3_area_laminin[,!(names(tumor_roi3_area_laminin) %in% cols_to_drop)]
tumor_roi4_area_laminin<-tumor_roi4_area_laminin[,!(names(tumor_roi4_area_laminin) %in% cols_to_drop)]
tumor_roi5_area_laminin<-tumor_roi5_area_laminin[,!(names(tumor_roi5_area_laminin) %in% cols_to_drop)]
tumor_roi6_area_laminin<-tumor_roi6_area_laminin[,!(names(tumor_roi6_area_laminin) %in% cols_to_drop)]
tumor_roi7_area_laminin<-tumor_roi7_area_laminin[,!(names(tumor_roi7_area_laminin) %in% cols_to_drop)]
# scatter plot to check correlation between variables
pairs(infil_roi1_area_laminin)
# figure margins too large message
par(mfcol=c(5,3),mai=c(0.5,0.5,0.5,0))
pairs(infil_roi1_area_laminin)
# opening a separate window did not help
dev.new(noRStudioGD = TRUE)
pairs(infil_roi1_area_laminin)
dev.off()
# now checking subsets of the dataframe
pairs(infil_roi1_area_laminin[,c(20:25)]) # inertia
# inertia_tensor.0.1=inertia_tensor.1.0
pairs(infil_roi1_area_laminin[,c(42:57)]) # moments_central
pairs(infil_roi1_area_laminin[,79:ncol(infil_roi1_area_laminin)]) # perimeter onwards
# perimeter = perimeter_crofton
pairs(infil_roi1_area_laminin[,c(1:10)]) # area
# new set of cols to drop
infil_roi1_area_laminin<-infil_roi1_area_laminin[,-c(2,3,4,6,22,80)]
infil_roi2_area_laminin<-infil_roi2_area_laminin[,-c(2,3,4,6,22,80)]
infil_roi3_area_laminin<-infil_roi3_area_laminin[,-c(2,3,4,6,22,80)]
infil_roi4_area_laminin<-infil_roi4_area_laminin[,-c(2,3,4,6,22,80)]
infil_roi5_area_laminin<-infil_roi5_area_laminin[,-c(2,3,4,6,22,80)]
tumor_roi1_area_laminin<-tumor_roi1_area_laminin[,-c(2,3,4,6,22,80)]
tumor_roi2_area_laminin<-tumor_roi2_area_laminin[,-c(2,3,4,6,22,80)]
tumor_roi3_area_laminin<-tumor_roi3_area_laminin[,-c(2,3,4,6,22,80)]
tumor_roi4_area_laminin<-tumor_roi4_area_laminin[,-c(2,3,4,6,22,80)]
tumor_roi5_area_laminin<-tumor_roi5_area_laminin[,-c(2,3,4,6,22,80)]
tumor_roi6_area_laminin<-tumor_roi6_area_laminin[,-c(2,3,4,6,22,80)]
tumor_roi7_area_laminin<-tumor_roi7_area_laminin[,-c(2,3,4,6,22,80)]
pairs(infil_roi1_area_laminin[,c(3:6)]) # bounding boxes
# bbox.0 = bbox.2, bbox.1=bobx.3
pairs(infil_roi1_area_laminin[,c(7:10)]) # centroids
# let us build up the argument slowly, beginning with area
boxplot(infil_roi1_area_laminin$area,tumor_roi1_area_laminin$area,ylim=c(0,10000))
# add labels
infil_roi1_area_laminin$sample<-"roi1"
infil_roi1_area_laminin$class<-"infil"

infil_roi2_area_laminin$sample<-"roi2"
infil_roi2_area_laminin$class<-"infil"

infil_roi3_area_laminin$sample<-"roi3"
infil_roi3_area_laminin$class<-"infil"

infil_roi4_area_laminin$sample<-"roi4"
infil_roi4_area_laminin$class<-"infil"

infil_roi5_area_laminin$sample<-"roi5"
infil_roi5_area_laminin$class<-"infil"

tumor_roi1_area_laminin$sample<-"roi1"
tumor_roi1_area_laminin$class<-"tumor"

tumor_roi2_area_laminin$sample<-"roi2"
tumor_roi2_area_laminin$class<-"tumor"

tumor_roi3_area_laminin$sample<-"roi3"
tumor_roi3_area_laminin$class<-"tumor"

tumor_roi4_area_laminin$sample<-"roi4"
tumor_roi4_area_laminin$class<-"tumor"

tumor_roi5_area_laminin$sample<-"roi5"
tumor_roi5_area_laminin$class<-"tumor"

tumor_roi6_area_laminin$sample<-"roi6"
tumor_roi6_area_laminin$class<-"tumor"

tumor_roi7_area_laminin$sample<-"roi7"
tumor_roi7_area_laminin$class<-"tumor"

# now attach all the area columns
area_infil_1 <- infil_roi1_area_laminin[,c(1,76,77)]
area_infil_2 <- infil_roi2_area_laminin[,c(1,76,77)]
area_infil_3 <- infil_roi3_area_laminin[,c(1,76,77)]
area_infil_4 <- infil_roi4_area_laminin[,c(1,76,77)]
area_infil_5 <- infil_roi5_area_laminin[,c(1,76,77)]
area_tumor_1 <- tumor_roi1_area_laminin[,c(1,76,77)]
area_tumor_2 <- tumor_roi2_area_laminin[,c(1,76,77)]
area_tumor_3 <- tumor_roi3_area_laminin[,c(1,76,77)]
area_tumor_4 <- tumor_roi4_area_laminin[,c(1,76,77)]
area_tumor_5 <- tumor_roi5_area_laminin[,c(1,76,77)]
area_tumor_6 <- tumor_roi6_area_laminin[,c(1,76,77)]
area_tumor_7 <- tumor_roi7_area_laminin[,c(1,76,77)]
area <- rbind(area_infil_1,area_infil_2,area_infil_3,area_infil_4,
              area_infil_5,area_tumor_1,area_tumor_2,area_tumor_3,
              area_tumor_4,area_tumor_5,area_tumor_5,area_tumor_6,
              area_tumor_7)
area$sample <- as.factor(area$sample)
area$class <- as.factor(area$class)

# plot the area
ggplot(data=area,aes(x=class,y=area,col=sample)) + geom_boxplot()
area_laminin <-area
names(area_laminin)<-c("area","ROI","region")

ggplot(data=area_laminin,aes(x=region,y=area,color=ROI)) + geom_boxplot()+
  ylim(0,10000)+theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
                      axis.title.x = element_text(color="black", size=22, face="bold"),
                      axis.title.y = element_text(color="black", size=22, face="bold"))+
  theme(axis.text = element_text(color="black",size=20,face="bold",angle=0))+
  ggtitle("area occupied by laminin") + theme(panel.background = element_blank())
# combining all the rois
laminin_all <- rbind(infil_roi1_area_laminin,infil_roi2_area_laminin,
                       infil_roi3_area_laminin,infil_roi4_area_laminin,
                      infil_roi5_area_laminin, tumor_roi1_area_laminin,
                       tumor_roi2_area_laminin,tumor_roi3_area_laminin,
                       tumor_roi4_area_laminin,tumor_roi5_area_laminin,
                       tumor_roi6_area_laminin,tumor_roi7_area_laminin)
laminin_all$marker<-"laminin"
write.csv(laminin_all,"laminin_all.csv",row.names = FALSE)
# update infil to margin
area_laminin<-area_laminin %>%
  mutate(across('region', str_replace, 'infil', 'margin'))
# updated plots
ggplot(data=area_laminin,aes(x=region,y=area,color=ROI)) + geom_boxplot()+
  ylim(0,10000)+theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
                      axis.title.x = element_text(color="black", size=22, face="bold"),
                      axis.title.y = element_text(color="black", size=22, face="bold"))+
  theme(axis.text = element_text(color="black",size=20,face="bold",angle=0))+
  ggtitle("area occupied by laminin") + theme(panel.background = element_blank())
# calculate means of margin and tumor
area_laminin%>%group_by(region)%>%summarise(mean=mean(area),sd=sd(area))
area_laminin_mean<-area_laminin%>%group_by(ROI,region)%>%summarise(mean=mean(area),sd=sd(area))
area_laminin_mean<-as.data.frame(area_laminin_mean)
t.test(mean ~ region, data = area_laminin_mean)
# create annotated plot
ggplot(data=area_laminin,aes(x=region,y=sqrt(area),color=ROI)) + geom_boxplot()+
  theme(plot.title = element_text(color="black", size=22, face="bold.italic"),
        axis.title.x = element_text(color="black", size=22, face="bold"),
        axis.title.y = element_text(color="black", size=22, face="bold"))+
  theme(axis.text = element_text(color="black",size=20,face="bold",angle=0))+
  ggtitle("area occupied by Laminin") + theme(panel.background = element_blank())+
  annotate("text", x = 1, y = 220, label = "Wilcoxon rank sum test")+ 
  annotate("text", x = 1, y = 206, label = "                  mean")+
  annotate("text", x = 1, y = 190, label = "        margin 1288.44")+ 
  annotate("text", x = 1, y = 175, label = "        tumor  1024.42")+
  annotate("text", x = 1, y = 160, label = "        (p-value=0.202)")
  