library(tidyverse)
library(ggplot2)
library(pheatmap)


core<-rbind(core_ub18_66043_6_full,
core_ub18_66043_7_full,
            core_ub18_66043_8A_full,
            core_ub19_49455_2_full,
            core_ub19_49455_3_full,
            core_ub19_49455_4_full,
            core_ub19_49455_5_full,
            core_ub19_52388_5_full)

margin<-rbind(margin_ub19_49455_1_full,
              margin_ub19_52388_2_full)
# scale the datasets
core$region <- "tumor"
margin$region <- "margin"
core$region <- as.factor(core$region)
margin$region <- as.factor(margin$region)

core_sample<-core$sample
margin_sample<-margin$sample


# core numerics

core_nums <- core%>%select(-index,-region,-sample)
margin_nums <- margin%>%select(-index,-region,-sample)

core_nums<- as.data.frame(scale(core_nums))
margin_nums <- as.data.frame(scale(margin_nums))

apply(core_nums,2,function(x) sum(is.na(x)))
apply(margin_nums,2,function(x) is.na(x))

core_nums<- core_nums%>%select(-red_Intensity.Min)
margin_nums<- margin_nums%>%select(-red_Intensity.Min)

row.names(core_nums) <-core$sample
row.names(margin_nums) <- margin$sample
#core_nums_t<-t(core_nums)
#colnames(core_nums_t) <- core$region
#core_nums_t<- as.data.frame(core_nums_t)
combined_df <- rbind(core_nums,margin_nums) 
combined_df <- as.data.frame(t(combined_df))
annot_names<- c(core$region,margin$region)
annot_names <- as.data.frame(annot_names)

row.names(annot_names) <-colnames(combined_df)

pheatmap(combined_df,annotation_col = annot_names,
         annotation_names_col = FALSE,
         annotation_names_row = FALSE,
         show_rownames = T,show_colnames = F,
         cluster_cols = T,cluster_rows = F)


# employ WGCNA principles
high_order_df<-(combined_df)*100000 # enhance the heatmap signal. multiplying works better than exponentiation

# the color palette has to be readjusted to take the enhanced values
makeColorRampPalette <- function(colors, cutoff.fraction, num.colors.in.palette) {
  stopifnot(length(colors) == 4)
  ramp1 <- colorRampPalette(colors[1:2])(num.colors.in.palette * cutoff.fraction)
  ramp2 <- colorRampPalette(colors[3:4])(num.colors.in.palette * (1 - cutoff.fraction))
  return(c(ramp1, ramp2))
}
color.divisions <- 100
color <- makeColorRampPalette(c("navy", "white","white", "red"), 
                              cutoff.distance / max(distmat), color.divisions)
pheatmap(high_order_df,annotation_col = annot_names,
         annotation_names_col = FALSE,
         annotation_names_row = FALSE,
         breaks= seq(-100000,100000,length.out=(color.divisions + 1)), # this readjusts the color palette to the enhanced settings
         show_rownames = T,show_colnames = F,
         cluster_cols = T,cluster_rows = T)

contrastive_df <- combined_df-high_order_df

pheatmap(contrastive_df,annotation_col = annot_names,
         annotation_names_col = FALSE,
         annotation_names_row = FALSE,
         show_rownames = T,show_colnames = F,
         cluster_cols = T,cluster_rows = F)
