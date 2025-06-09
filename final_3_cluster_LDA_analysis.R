library(tidyverse)
library(ggplot2)
library(caret)
library(MASS)
library(FactoMineR)
library(factoextra)
library(ggord)
library(regclass)
library(pheatmap)
library(ggpubr)
library(GGally)
library(corrplot)
theme_set(theme_classic())
# LDA is a dimensionality reduction technique
set.seed(42)

old.LDA.model <- readRDS("all_df_lda_model.rds")
old.predictions<-old.LDA.model%>%predict(comp_size_unscaled_df)
# we shall create a LDA model with the new classification
# and only from the laminin,claudin and shape features
composition_core<- rbind(
  ub18_66043_6_antigen_full,
  ub18_66043_7_antigen_full,
  ub18_66043_8A_set3_512px_mask_antigen_full,
  ub19_49455_2_antigen_full,
  ub19_49455_3_antigen_full,
  ub19_49455_4_antigen_full,
  ub19_49455_5_antigen_full,
  ub19_52388_1B_antigen_full,
  ub19_52388_5_antigen_full)
  
shape_core <- rbind(
  ub18_66043_6_set3_512px_shape_full,
  ub18_66043_7_set3_512px_shape_full,
  ub18_66043_8A_set3_512px_shape_full,
  ub19_49455_2_set3_512px_shape_full,
  ub19_49455_3_set3_512px_shape_full,
  ub19_49455_4_set3_512px_shape_full,
  ub19_49455_5_set3_512px_shape_full,
  ub19_52388_1B_set3_512px_shape_full,
  ub19_52388_5_set3_512px_shape_full
)
composition_margin<-ub19_49455_1_antigen_full
composition_normal<-ub19_52388_2_antigen_full

shape_margin <-ub19_49455_1_set3_512px_shape_full
shape_normal<- ub19_52388_2_set3_512px_shape_full

# model building with set1 cluster
row.names(new_clusters_all_df) <- new_clusters_all_df[,1]
new_clusters_all_df <- new_clusters_all_df[,-1]
new_clusters_all_df$newcluster <- as.factor(new_clusters_all_df$newcluster)

endothelium<-new_clusters_all_df%>%dplyr::select(!starts_with("red"))
# data partition
#create ID column
endothelium$id <- 1:nrow(endothelium)
#use 70% of data set as training set and 30% as test set 
train <- endothelium %>% dplyr::sample_frac(0.70)
train_id <- train$id


test  <- dplyr::anti_join(endothelium, train, by = 'id')
test_id <- test$id

train_id <- train$id
test_id <- test$id

train <- train%>%dplyr::select(-id)
test <- test%>%dplyr::select(-id)

# build LDA on the new train set
# create the initial model
set1.lda.model <- lda(newcluster~.,train)  
# Make predictions
set1.predictions <- set1.lda.model %>% predict(test[,-42])
# confusion matrix
confusionMatrix(set1.predictions$class,test$newcluster)
table(test$newcluster)
# plot the biplot
ggord(set1.lda.model,
      as.factor(train$newcluster),repel=TRUE,size=2,
      alpha_el=0.1,alpha=0.8)
# plot the predictions
x<-set1.predictions$x%>%as.data.frame()
x$predictions <- as.factor(set1.predictions$class)
ggplot(x,aes(x=LD1,y=LD2,col=predictions))+
  geom_point()+
  xlim(-8.0,5.0)+ylim(-8.0,10)

# common data set for composition and shape for set3
composition_core$region <- "tumor"
shape_core$region <- "tumor"

composition_margin$region <- "margin"
shape_margin$region <- "margin"

composition_normal$region <- "normal"
shape_normal$region <- "normal"

set3_composition <-rbind(composition_core,composition_margin,composition_normal)
set3_shape <- rbind(shape_core,shape_margin,shape_normal)
set3_composition <- set3_composition[,-1]
set3_shape <- set3_shape[,-c(1:2)]
set3_shape$sample<-sub("set3_512px_mask_","",set3_shape$sample)
set3_composition$sample <- change_name(set3_composition$sample)

set3_all<-merge(set3_composition,set3_shape,by=c("sample","region"))
# now from set3 we will remove red (PDGFRb)

set3_all$region <- as.factor(set3_all$region)
set3_endothelium<- set3_all%>%dplyr::select(!starts_with("red"))
set1_endothelium<- comp_size_unscaled_df%>%dplyr::select(!starts_with("red"))
# now we will compare set3 with set 1 endothelium
set1_endothelium <- set1_endothelium[,c(1,5:46)]
# choose the same subset of endothelium for set1 and set3
set3_endothelium<- set3_endothelium[,colnames(set3_endothelium) %in% colnames(set1_endothelium)]

set1_endothelium$set <- as.factor("set1")
set3_endothelium$set <- as.factor("set3")
# joint the 2 sets
endothelium.unscaled <- rbind(set1_endothelium,set3_endothelium)
endothelium.dab <- endothelium.unscaled[,c(3:10,44)]
endothelium.green <- endothelium.unscaled[,c(11:18,44)] 
endothelium.shape <- endothelium.unscaled[,c(19:27,44)] 

# modify the dab dataset for boxplot
endothelium.dab<-endothelium.dab%>%gather(!set,key="features",value="value")
ggplot(endothelium.dab,aes(x=features,y=value,col=set))+
  geom_boxplot()+
  facet_wrap(~features,scales = "free_y")+
  theme(axis.text.x = element_blank())+
  ggtitle("DAB")+theme(panel.background = element_blank())+
  theme(axis.text.y = element_text(color="black",size=9,face="bold",angle=0),
        axis.title.x = element_text(color="black", size=12, face="bold"),
        axis.title.y = element_text(color="black", size=12, face="bold"))+
  theme(legend.position = "none")+
  theme(
    strip.text.x = element_text(
      size = 9, color = "black", face = "bold.italic"
    ),
    strip.text.y = element_text(
      size = 9, color = "black", face = "bold.italic"
    )
  )

endothelium.green<-endothelium.green%>%gather(!set,key="features",value="value")
ggplot(endothelium.green,aes(x=features,y=value,col=set))+
  geom_boxplot()+
  facet_wrap(~features,scales = "free_y")+
  theme(axis.text.x = element_blank())+
  ggtitle("Green")+theme(panel.background = element_blank())+
  theme(axis.text.y = element_text(color="black",size=9,face="bold",angle=0),
        axis.title.x = element_text(color="black", size=12, face="bold"),
        axis.title.y = element_text(color="black", size=12, face="bold"))+
  theme(legend.position = "none")+
  theme(
    strip.text.x = element_text(
      size = 9, color = "black", face = "bold.italic"
    ),
    strip.text.y = element_text(
      size = 9, color = "black", face = "bold.italic"
    )
  )

endothelium.shape<-endothelium.shape%>%gather(!set,key="features",value="value")

ggplot(endothelium.shape,aes(x=features,y=value,col=set))+
  geom_boxplot(position = position_dodge2(preserve = "single"),alpha=0.6)+
  facet_wrap(~features,scales = "free_y")+
  theme(axis.text.x = element_blank())+
  ggtitle("Shape")+theme(panel.background = element_blank())+
  theme(axis.text.y = element_text(color="black",size=9,face="bold",angle=0),
        axis.title.x = element_text(color="black", size=12, face="bold"),
        axis.title.y = element_text(color="black", size=12, face="bold"))+
  theme(legend.position = "none")+
  theme(
    strip.text.x = element_text(
      size = 9, color = "black", face = "bold.italic"
    ),
    strip.text.y = element_text(
      size = 9, color = "black", face = "bold.italic"
    )
  )
# Building the LDA model based on Endothelium and BM
endothelium_metadata <- endothelium.unscaled[,c(1,2,44)]
endothelium_numdata <- endothelium.unscaled[,-c(1,2,44)]
endothelium.num.scaled <- scale(endothelium_numdata)
endothelium.num.scaled <- as.data.frame(endothelium.num.scaled)

# create the scaled data set
endothelium.scaled <- cbind(endothelium_metadata,endothelium.num.scaled)

# divide data set into test and train

#create ID column
endothelium.scaled$id <- 1:nrow(endothelium.scaled)
endothelium.scaled$region <- as.factor(endothelium.scaled$region)
endothelium.scaled$set <- as.factor(endothelium.scaled$set)
#use 70% of data set as training set and 30% as test set 
train <- endothelium.scaled %>% dplyr::sample_frac(0.70)
train_id <- train$id


test  <- dplyr::anti_join(endothelium.scaled, train, by = 'id')
test_id <- test$id

train.lda <- train%>%dplyr::select(-c(id,sample,set))
test.lda <- test%>%dplyr::select(-c(id,sample,set))
table(train.lda$region)

# build LDA on the new train set
# Calculate group sizes
group_sizes <- table(train.lda$region)

# Calculate prior probabilities
prior_probabilities <- group_sizes / sum(group_sizes)
prior_probabilities <- as.vector(prior_probabilities)
# create the LDA model
lda.model <- lda(region~.,train.lda,prior=prior_probabilities)  
lda.model
# Make predictions
predictions <- lda.model %>% predict(train.lda[,-1])
# confusion matrix
original_conf_matrix <- confusionMatrix(predictions$class,train.lda$region)
original_conf_matrix
table(train.lda$region)
# based on the initial model new classification 
train.lda$predictions <- predictions$class
tumor <- train.lda%>%filter(train.lda$region == "tumor" & train.lda$predictions == "tumor")
margin <- train.lda%>%filter(train.lda$region == "margin" & train.lda$predictions == "margin")
normal <- train.lda%>%filter(train.lda$region == "normal" & train.lda$predictions == "normal")
# revised train dataset for LDA
revised.train <- rbind(tumor,margin,normal)
revised.train <- revised.train%>%dplyr::select(-predictions)
revised.lda.model <- lda(region ~.,revised.train)
# predictions revised train
revised.predictions <- revised.lda.model%>%predict(revised.train)
# confusion matrix for revised data set
revised_conf_matrix <- confusionMatrix(revised.predictions$class,revised.train$region)
revised_conf_matrix
# based on above criteria, we shall classify the entire dataset
endothelium.scaled.lda <- endothelium.scaled%>%dplyr::select(-c(sample,id,set))
str(endothelium.scaled.lda)
endothelium.scaled.model <- lda(region ~.,endothelium.scaled.lda)
endothelium.scaled.predictions <- endothelium.scaled.model%>%predict(endothelium.scaled.lda)
endothelium.scaled.lda$predictions <- endothelium.scaled.predictions$class
# full data confusion matrix
endothelium.conf.matrix <- confusionMatrix(endothelium.scaled.lda$predictions,endothelium.scaled.lda$region)
endothelium.conf.matrix

# reconstruct the data sets
endothelium.scaled.lda$sample <- endothelium.scaled$sample
endothelium.scaled.lda$set <- endothelium.scaled$set
endothelium.scaled.lda <- endothelium.scaled.lda%>%relocate(sample,region,set,predictions)

tumor <- endothelium.scaled.lda%>%dplyr::filter(region == "tumor" & predictions == "tumor")
margin <- endothelium.scaled.lda%>%dplyr::filter(region == "margin" & predictions == "margin")
normal <- endothelium.scaled.lda%>%dplyr::filter(region == "normal" & predictions == "normal")

confirmed.class.df <- rbind(tumor,margin,normal)
common <- dplyr::anti_join(endothelium.scaled.lda, confirmed.class.df, by = 'sample')

# final classification
tumor <- tumor%>%dplyr::select(-region,-predictions)
tumor$region <- as.factor("tumor")



normal <- normal%>%dplyr::select(-region,-predictions)
normal$region <- as.factor("normal")

margin <- margin%>%dplyr::select(-region,-predictions)
margin$region <- as.factor("margin")

common <- common%>%dplyr::select(-region,-predictions)
common$region <- as.factor("common")

# rejoining reclassified datasets

reclassifed.df <- rbind(tumor,normal,margin,common)
reclassifed.df <- reclassifed.df%>%relocate(sample,set,region)
reclassifed.df$set <- as.factor(reclassifed.df$set)
reclassifed.df$region <- as.factor(reclassifed.df$region)
final.lda.model <- lda(region ~.,reclassifed.df[,c(3:44)])
final.lda.model

# plot the model
# plot the biplot
ggord(final.lda.model,
      as.factor(reclassifed.df$region),repel=TRUE,size=2,
      alpha_el=0.1,alpha=0.8)
# ordination plot with reduced tumor samples
resample.tumor<-sample(1:nrow(tumor),2000)

resample.tumor <- tumor[resample.tumor,]

names(common)[!names(common)%in%names(resample.tumor)]
resample.common <- common%>%dplyr::select(-source,-image,-vessel,-predictions)
resample.common$region <- "common"
resample.df <- rbind(resample.tumor,normal,margin,resample.common)
resample.df <- resample.df%>%relocate(region,.after=set)
resample.lda <- lda(region ~.,resample.df[,c(3:44)])
# create a category for common vs non-common
facet_on_common <- ifelse(resample.df$region=="common",1,0)
resample.mat <- as.matrix(resample.df[,4:44])
resample.coords <- as.data.frame(resample.mat%*%resample.lda$scaling)
resample.coords$region <- as.factor(resample.df$region)
# plotting only the defined regions
ggplot(resample.coords[resample.coords$region!="common",],aes(x=LD1,y=LD2))+ 
  geom_point(aes(colour=region),size=0.5,show.legend = FALSE)+scale_color_manual(values=c('#ff33ff', 'green', 'blue'))+
  scale_shape_manual('Groups', values = c(16,17,18))+xlim(-4,5)+
  theme(axis.text = element_text(color="black",size=20,face="bold",angle=0),
        axis.title.x = element_text(color="black", size=22, face="bold"),
        axis.title.y = element_text(color="black", size=22, face="bold"))+
  theme(legend.title = element_text(size=16),
        legend.text = element_text(size=16))+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 2))

 

ggplot(resample.coords[resample.coords$region =="common",],aes(x=LD1,y=LD2))+ 
  geom_point(colour="#B4B4B4FF",size=0.5,show.legend = F)+xlim(-4,5)+ylim(-4,4)+
  scale_shape_manual('Groups', values = c(16,17,18))+
  theme(axis.text = element_text(color="black",size=20,face="bold",angle=0),
        axis.title.x = element_text(color="black", size=22, face="bold"),
        axis.title.y = element_text(color="black", size=22, face="bold"))+
  theme(legend.title = element_text(size=16),
        legend.text = element_text(size=16))+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 2))


ggord(resample.lda,grp_in=resample.df$region,
      repel=TRUE,size=1,facet= FALSE,
      alpha_el=0.2,alpha=0.6,shape=2,
      cols = c('#ff33ff', 'green', 'blue','#f1e8e7'))
  
  
# scalings
scalings <- as.data.frame(final.lda.model$scaling)
# plot the biplot
ggord(final.lda.model,
      as.factor(reclassifed.df$region),repel=TRUE,size=1.6,
      alpha_el=0.08,alpha=0.9,
      cols = c('#ff33ff', 'green', 'blue','#B4B4B4FF'))+
    scale_shape_manual('Groups', values = c(16,17,18,19))+
  theme(axis.text = element_text(color="black",size=20,face="bold",angle=0),
        axis.title.x = element_text(color="black", size=22, face="bold"),
        axis.title.y = element_text(color="black", size=22, face="bold"))+
  theme(legend.title = element_text(size=16),
        legend.text = element_text(size=16))


# plot the biplot
ggord(final.lda.model,
      as.factor(reclassifed.df$region),repel=TRUE,size=1.6,
      alpha_el=0.08,alpha=0.9,
      cols = c('#ff33ff', 'green', 'blue','#B4B4B4FF'))+
  scale_shape_manual('marvel', values = c(16,17,18,19))+
  theme(axis.text = element_text(color="black",size=20,face="bold",angle=0),
        axis.title.x = element_text(color="black", size=22, face="bold"),
        axis.title.y = element_text(color="black", size=22, face="bold"))+
  theme(legend.title = element_text(size=16),
        legend.text = element_text(size=16))

# model stats
plot(final.lda.model)
final.predictions <- final.lda.model%>%predict(reclassifed.df[,c(4:44)])
conf.mat.lda<-confusionMatrix(final.predictions$class,reclassifed.df$region)
conf.mat.lda
# Based on conf.mat.lda we have the distribution of the 
# common vessels in the various regions
# based on values of conf.mat.lda, we have
# either common recognized as tumor/normal/margin
# or tumor/normal/margin recognized as common
# tumor normal margin LDA model
defined_df <- rbind(tumor,normal,margin)
defined_df <- defined_df%>%relocate(sample,set,region)
defined.lda.model <- lda(region ~.,defined_df[,c(3:44)])

defined.predictions <- defined.lda.model%>%predict(defined_df[,c(3:44)])
conf.defined.lda<-confusionMatrix(defined.predictions$class,defined_df$region)
conf.defined.lda
table(defined_df$region)

# measure the CD31 and Laminin for set3
reclassified_set3 <- reclassifed.df%>%dplyr::filter(set=="set3")
reclassified_set3_metadata <- reclassified_set3[,c(1:3)]
names(reclassified_set3_metadata) <- c("sample","set","predictions")

# merge set3_all and reclassified_metadata
set3_all_reclassified <- merge(set3_all,reclassified_set3_metadata,by="sample")
set3_all_reclassified <- set3_all_reclassified%>%relocate(sample,region,set,predictions)
# analyze Laminin
my_comparisons <- list( c("tumor", "normal"), c("tumor", "margin"), c("normal", "margin"),c("tumor","common") )
ggplot(data=set3_all_reclassified,aes(x= predictions,y=green_Intensity.Mean,fill=predictions))+ 
  geom_boxplot()+
  theme(panel.background = element_blank())+ylim(0,220)+
  theme(title = element_text(color="black", size=18, face="bold"),
        axis.text.x = element_text(color="black",size=18,face="bold",angle=45,vjust=.8, hjust=0.8),
        axis.text.y = element_text(color="black",size=18,face="bold",angle=0),
        axis.title.x = element_text(color="black", size=22, face="bold"),
        axis.title.y = element_text(color="black", size=22, face="bold"))+
  labs(title="Laminin",x = "region", y = "mean intensity")+theme(legend.position = "none")+
  stat_compare_means(method="t.test",label.x=1.0,
                     label = "p.signif",
                     label.sep=" ",comparisons = my_comparisons,
                     label.y= c(130,150,190,220 ),
                     size=6,ref.group = "tumor")

# analyze CD31
my_comparisons <- list( c("tumor", "normal"), c("tumor", "margin"), c("normal", "margin"),c("tumor","common") )
ggplot(data=set3_all_reclassified,aes(x= predictions,y=dab_Intensity.Mean,fill=predictions))+ 
  geom_boxplot()+
  theme(panel.background = element_blank())+ylim(0,220)+
  theme(title = element_text(color="black", size=18, face="bold"),
        axis.text.x = element_text(color="black",size=18,face="bold",angle=45,vjust=.8, hjust=0.8),
        axis.text.y = element_text(color="black",size=18,face="bold",angle=0),
        axis.title.x = element_text(color="black", size=22, face="bold"),
        axis.title.y = element_text(color="black", size=22, face="bold"))+
  labs(title= "CD31", x = "region", y = "mean intensity")+theme(legend.position = "none")+
  stat_compare_means(method="t.test",label.x=1.0,
                     label = "p.signif",
                     label.sep=" ",comparisons = my_comparisons,
                     label.y= c(150,170,190,210 ),
                     size=6,ref.group = "tumor")

# analyze PDGFRb
my_comparisons <- list( c("tumor", "normal"), c("tumor", "margin"), c("normal", "margin"),c("tumor","common") )

set3_all_reclassified$predictions <- factor(set3_all_reclassified$predictions,levels=c("common","margin","normal","tumor"))

ggplot(data=set3_all_reclassified,aes(x= predictions,y=red_Intensity.Mean,fill=predictions))+ 
  geom_boxplot()+
  theme(panel.background = element_blank())+ylim(0,220)+
  theme(title = element_text(color="black", size=18, face="bold"),
        axis.text.x = element_text(color="black",size=18,face="bold",angle=45,vjust=.8, hjust=0.8),
        axis.text.y = element_text(color="black",size=18,face="bold",angle=0),
        axis.title.x = element_text(color="black", size=22, face="bold"),
        axis.title.y = element_text(color="black", size=22, face="bold"))+
  labs(title= "PDGFRb", x = "region", y = "mean intensity")+theme(legend.position = "none")+
  stat_compare_means(method="t.test",label.x=1.5,
                     label = "p.signif",
                     label.sep=" ",comparisons = my_comparisons,
                     label.y= c(150,170,190,210 ),
                     size=6,ref.group = "tumor")

# set3 new data
metadata <- set3_new_metadata
metadata <- metadata%>%group_by(predictions,sample)%>%summarize(count=n())%>%arrange(desc(count))
metadata.tumor <- metadata%>%filter(predictions=="tumor")
metadata.common <- metadata%>%filter(predictions == "common")
metadata.margin <- metadata%>%filter(predictions == "margin")
metadata.normal <- metadata%>%filter(predictions == "normal")
metadata.tumor[metadata.tumor$sample%>%str_detect("ub19_49455_3"),]%>%print(n=20)
metadata.tumor[metadata.tumor$sample%>%str_detect("ub19_49455_2"),]%>%print(n=20)
# analysis of correlations
set3_all.tumor <- set3_all%>%filter(region =="tumor")
set3_all.normal <- set3_all%>%filter(region =="normal")
set3_all.margin <- set3_all%>%filter(region == "margin")

set3_all.tumor <- set3_all.tumor[,3:165]
set3_all.normal <- set3_all.normal[,3:165]
set3_all.margin <- set3_all.margin[,3:165]

ggpairs(set3_all.tumor) # not effective due to memory
sample.x <- 1:nrow(set3_all.tumor)
sample.y <- sample(sample.x,800,replace = F)
set3_all.tumor<-set3_all.tumor[sample.y,]
pheatmap(set3_all.tumor,show_rownames = F,show_colnames = F)
ggpairs(set3_all.tumor)
# all data
sample.all <- 1:nrow(set3_all)
sample.all <- sample(sample.all,800,replace = F)
set3_all.subset <-set3_all[sample.all,]
set3_all.subset <- set3_all.subset[,3:165]

pheatmap(set3_all.subset,show_rownames = F,show_colnames = F)
# enhance color of heat map
color.divisions <- 100
pheatmap(set3_all.subset,
         breaks= seq(-1000,1000,length.out=(color.divisions + 1)), # this readjusts the color palette to the enhanced settings
         show_rownames = F,
         show_colnames = F)

pheatmap(set3_all.tumor,
         breaks= seq(-1000,1000,length.out=(color.divisions + 1)), # this readjusts the color palette to the enhanced settings
         show_rownames = F,
         show_colnames = F)

pheatmap(set3_all.normal,
         breaks= seq(-1000,1000,length.out=(color.divisions + 1)), # this readjusts the color palette to the enhanced settings
         show_rownames = F,
         show_colnames = F)

pheatmap(set3_all.margin,
         breaks= seq(-1000,1000,length.out=(color.divisions + 1)), # this readjusts the color palette to the enhanced settings
         show_rownames = F,
         show_colnames = F)
# min max scaling the datasets
process <- preProcess(as.data.frame(set3_all.subset), method=c("range"))
set3_all.subset <- predict(process, as.data.frame(set3_all.subset))
pheatmap(set3_all.subset,
         breaks= seq(-1,1,length.out=(color.divisions + 1)), # this readjusts the color palette to the enhanced settings
         show_rownames = F,
         show_colnames = T,
         fontsize = 6)
corrplot(cor(set3_all.subset),tl.cex=0.2,tl.pos='ld',
         type='lower',tl.srt=30)
# we will start rationalizing and reducing variables
set3_all.subset%>%dplyr::select(-ends_with("Range"))%>%dim()
set3_all.subset <- set3_all.subset%>%dplyr::select(-ends_with("Range"))
corrplot(cor(set3_all.subset),tl.cex=0.2,tl.pos='ld',
         type='lower',tl.srt=30)
set3_all.subset%>%dplyr::select(-contains("HuMoments"))%>%dim()
set3_all.subset <- set3_all.subset%>%dplyr::select(-contains("HuMoments"))
corrplot(cor(set3_all.subset),tl.cex=0.2,tl.pos='ld',
         type='lower',tl.srt=30)
set3_all.subset%>%dplyr::select(-contains("FSD"))%>%dim()
set3_all.subset <- set3_all.subset%>%dplyr::select(-contains("FSD"))
set3_all.subset <- set3_all.subset%>%dplyr::select(-contains("Min"))
set3_all.subset <- set3_all.subset%>%dplyr::select(-contains("IQR"))
set3_all.subset <- set3_all.subset%>%dplyr::select(-green_Gradient.Mag.Std,-green_Gradient.Mag.Skewness,-green_Gradient.Mag.Kurtosis)
set3_all.subset <- set3_all.subset%>%dplyr::select(-red_Haralick.IMC1.Mean,-red_Haralick.IMC2.Mean)

set3_all.subset <- set3_all.subset%>%dplyr::select(-contains("HuMoments"))
set3_all.subset <- set3_all.subset%>%dplyr::select(-ends_with("Range"))
# heatmap with the refined variables
pheatmap(set3_all.subset,
         breaks= seq(0,1,length.out=(color.divisions + 1)), # this readjusts the color palette to the enhanced settings
         show_rownames = F,
         show_colnames = T,
         fontsize = 6)
# create the annotation dataframe
sample.all <- 1:nrow(set3_all)
sample.all <- sample(sample.all,800,replace = F)
set3_all.subset <-set3_all[sample.all,]
set3_all.subset <- set3_all.subset%>%dplyr::arrange(region)
row.names(set3_all.subset) <- set3_all.subset$sample
annot_names <- as.data.frame(set3_all.subset$region)
names(annot_names) <- "region"
row.names(annot_names) <- set3_all.subset$sample
all(row.names(set3_all.subset)==row.names(annot_names))
set3_all.subset <- set3_all.subset[,3:165]
set3_all.subset <- set3_all.subset%>%dplyr::select(-contains("FSD"))
set3_all.subset <- set3_all.subset%>%dplyr::select(-contains("Min"))
set3_all.subset <- set3_all.subset%>%dplyr::select(-contains("IQR"))
set3_all.subset <- set3_all.subset%>%dplyr::select(-green_Gradient.Mag.Std,-green_Gradient.Mag.Skewness,-green_Gradient.Mag.Kurtosis)
set3_all.subset <- set3_all.subset%>%dplyr::select(-red_Haralick.IMC1.Mean,-red_Haralick.IMC2.Mean)

# min max scaling the datasets
process <- preProcess(as.data.frame(set3_all.subset), method=c("range"))
set3_all.subset <- predict(process, as.data.frame(set3_all.subset))
#plot
pheatmap(set3_all.subset,
         annotation_row = annot_names,
         annotation_names_col = FALSE,
         annotation_names_row = FALSE,
         breaks= seq(0,1,length.out=(color.divisions + 1)), # this readjusts the color palette to the enhanced settings
         show_rownames = F,
         show_colnames = T,
         fontsize = 5,
         cluster_cols = T,cluster_rows = F)

# find the top features of the model and then do
# the heatmap
s3_rec.tumor <- set3_all_reclassified%>%filter(predictions == "tumor")
s3_rec.normal <- set3_all_reclassified%>%filter(predictions == "normal")
s3_rec.margin <- set3_all_reclassified%>%filter(predictions == "margin")
s3_rec.common <- set3_all_reclassified%>%filter(predictions == "common")

# sample the tumor subset
sample.tumor <- 1:nrow(s3_rec.tumor)
sample.tumor <- sample(sample.tumor,400,replace = F)
s3_rec.tumor <-s3_rec.tumor[sample.tumor,]
# I am using the defined LDA model to find the key
# features on the basis which the bv are classified
# get the loadings LD1
loadings <- defined.lda.model$scaling
loadings <- as.data.frame(loadings)
loadings.ld1<-loadings%>%arrange(desc(abs(LD1)))
loadings.ld1 <- loadings.ld1%>%mutate(abs_LD1=abs(LD1))
loadings.ld1<-loadings.ld1%>%arrange(desc(abs_LD1))
plot(loadings.ld1$abs_LD1,ylab="LD1",main="LD1 loadings")

loadings.ld1$features<-row.names(loadings.ld1)
ggplot(data=loadings.ld1,aes(x=fct_reorder(features, abs_LD1, .desc = TRUE),y=abs_LD1))+geom_col()+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))+
  labs(x="features",y="LD1")+
  theme(panel.background = element_blank())


loadings.ld2<-loadings%>%arrange(desc(abs(LD2)))
loadings.ld2 <- loadings.ld2%>%mutate(abs_LD2=abs(LD2))
loadings.ld2<-loadings.ld2%>%arrange(desc(abs_LD2))
plot(loadings.ld2$abs_LD2,ylab="LD2",main="LD2 loadings")
loadings.ld2$features<-row.names(loadings.ld2)
ggplot(data=loadings.ld2,aes(x=fct_reorder(features, abs_LD2, .desc = TRUE),y=abs_LD2))+
  geom_col()+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))+
  labs(x="features",y="LD2")+
  theme(panel.background = element_blank())

# LDA based on final model including common
# get the loadings LD1
loadings <- final.lda.model$scaling
loadings <- as.data.frame(loadings)
loadings.ld1<-loadings%>%arrange(desc(abs(LD1)))
loadings.ld1 <- loadings.ld1%>%mutate(abs_LD1=abs(LD1))
loadings.ld1<-loadings.ld1%>%arrange(desc(abs_LD1))
plot(loadings.ld1$abs_LD1,ylab="LD1",main="LD1 loadings")

loadings.ld1$features<-row.names(loadings.ld1)
ggplot(data=loadings.ld1,aes(x=fct_reorder(features, abs_LD1, .desc = TRUE),y=abs_LD1))+geom_col()+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))+
  labs(x="features",y="LD1")+
  theme(panel.background = element_blank())

loadings.ld2<-loadings%>%arrange(desc(abs(LD2)))
loadings.ld2 <- loadings.ld2%>%mutate(abs_LD2=abs(LD2))
loadings.ld2 <-loadings.ld2%>%arrange(desc(abs_LD2))
plot(loadings.ld2$abs_LD2,ylab="LD2",main="LD2 loadings")
loadings.ld2$features<-row.names(loadings.ld2)
ggplot(data=loadings.ld2,aes(x=fct_reorder(features, abs_LD2, .desc = TRUE),y=abs_LD2))+
  geom_col()+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))+
  labs(x="features",y="LD2")+
  theme(panel.background = element_blank())
# find the top features
x<-merge(loadings.ld1,loadings.ld2,by="features")
x<-x%>%arrange(desc(abs_LD1))
features <- x$features
# create dataframe for heatmap
s3_rec <- rbind(s3_rec.tumor,s3_rec.margin,s3_rec.normal,s3_rec.common)
row.names(s3_rec) <- s3_rec$sample
annot_names <- as.data.frame(s3_rec$predictions)
names(annot_names) <- "region"
row.names(annot_names) <- row.names(s3_rec)
s3_rec <- s3_rec[,5:167]
s3_rec <- s3_rec[,features]
# min max scaling the datasets
process <- preProcess(as.data.frame(s3_rec), method=c("range"))
s3_rec <- predict(process, as.data.frame(s3_rec))
s3_rec <- s3_rec%>%select(-Shape.HuMoments6,-Shape.HuMoments5)
par(oma=c(3,3,3,3))
par(mar = c(5, 5, 5, 5))
# now build the features
jpeg(file="heatmap features.jpeg")
pheatmap(s3_rec,
         annotation_row = annot_names,
         annotation_names_col = FALSE,
         annotation_names_row = FALSE,
         scale = "column",
         breaks= seq(-1,1,length.out=(color.divisions + 1)), # this readjusts the color palette to the enhanced settings
         show_rownames = F,
         show_colnames = T,
         fontsize = 9,
         cluster_cols = T,cluster_rows = F,
         width=1,height =1,
         annotation_legend = T,
         legend = T)
dev.off() 
# final dataframe for plotting key differences
s3_rec_complete <- merge(s3_rec,annot_names,by=0)
write.csv(s3_rec_complete,"blood_vessel_key_properties.csv",row.names = F)
# features common to tumor and margin
reclassed_tumor_margin  <- reclassifed.df%>%filter(region=="tumor"|region=="margin")
reclassed_tumor_common  <- reclassifed.df%>%filter(region=="tumor"|region=="common")
reclassed_margin_common  <- reclassifed.df%>%filter(region=="margin"|region=="common")
reclassed_normal_common  <- reclassifed.df%>%filter(region=="normal"|region=="common")
reclassed_normal_tumor  <- reclassifed.df%>%filter(region=="normal"|region=="tumor")
reclassed_margin_normal  <- reclassifed.df%>%filter(region=="margin"|region=="normal")

x <- scale(reclassed_margin_normal[,c(4:44)])
x <- as.data.frame(x)
scaled_reclassed_tumor_margin <- cbind(reclassed_tumor_margin[,c(1:3)],x)
scaled_reclassed_tumor_common <- cbind(reclassed_tumor_common[,c(1:3)],x)
scaled_reclassed_margin_common <- cbind(reclassed_margin_common[,c(1:3)],x)
scaled_reclassed_normal_common <- cbind(reclassed_normal_common[,c(1:3)],x)
scaled_reclassed_margin_normal <- cbind(reclassed_margin_normal[,c(1:3)],x)

corrplot(cor(x),tl.cex = 0.4)
# comparision between sites
df.sites <- s3_rec_complete
df.sites.region <- df.sites$region
row.names(df.sites) <- df.sites$Row.names
df.sites <- df.sites%>%dplyr::select(-set,-region,-sample)
df.sites <- df.sites[,-1]
df.sites <- df.sites[,-40]
df.sites <- as.data.frame(t(df.sites))
df.sites <-df.sites[,1:70]
cor.df.sites <- cor(df.sites)
corrplot(cor.df.sites,tl.cex = 0.4)
# plot the key features
ggplot(data=comp_size_unscaled_df,aes(x=newcluster,y = Shape.HuMoments3))+ 
  geom_boxplot()+labs(x="region",y="Hu Moments 3")+
  theme(panel.background = element_blank())+
  theme(axis.text.x = element_text(color="black",size=20,face="bold",angle=45,vjust=1,hjust=1),
        axis.text.y = element_text(color="black",size=16,face="bold",angle=0),
        axis.title.x = element_text(color="black", size=22, face="bold"),
        axis.title.y = element_text(color="black", size=22, face="bold"))+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 2))

ggplot(data=comp_size_unscaled_df,aes(x=newcluster,y = green_Haralick.Entropy.Mean))+ 
  geom_boxplot()+labs(x="region",y="Laminin entropy")+
  theme(panel.background = element_blank())+
  theme(axis.text.x = element_text(color="black",size = 20,face="bold",angle=45,vjust=1,hjust=1),
        axis.text.y = element_text(color="black",size = 16,face="bold",angle=0),
        axis.title.x = element_text(color="black", size = 22, face="bold"),
        axis.title.y = element_text(color="black", size = 22, face="bold"))+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 2))

ggplot(data=comp_size_unscaled_df,aes(x=newcluster,y = dab_Haralick.Entropy.Mean))+ 
  geom_boxplot()+labs(x="region",y="CD31 entropy")+
  theme(panel.background = element_blank())+
  theme(axis.text.x = element_text(color="black",size=20,face="bold",angle=45,vjust=1,hjust=1),
        axis.text.y = element_text(color="black",size=16,face="bold",angle=0),
        axis.title.x = element_text(color="black", size=22, face="bold"),
        axis.title.y = element_text(color="black", size=22, face="bold"))+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 2))


ggplot(data=comp_size_unscaled_df,aes(x=newcluster,y = Size.Area))+ 
  geom_boxplot()+labs(x="region",y="Area")+
  theme(panel.background = element_blank())+
  theme(axis.text.x = element_text(color="black",size=20,face="bold",angle=45,vjust=1,hjust=1),
        axis.text.y = element_text(color="black",size=16,face="bold",angle=0),
        axis.title.x = element_text(color="black", size=22, face="bold"),
        axis.title.y = element_text(color="black", size=22, face="bold"))+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 2))
