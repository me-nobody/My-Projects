# plot ideogram
library(tidyverse)
library(RIdeogram)
# load standard data
data(human_karyotype, package="RIdeogram")
data(gene_density, package="RIdeogram")
data(Random_RNAs_500, package="RIdeogram")
heatmap.png <- read_csv("heatmap.png.csv", show_col_types = FALSE)
# refactor the dataset
data <- heatmap.png%>%separate(chromosome,into=c("chrom","bp"),sep=":")
data <- data%>%separate(bp,into=c("Start","End"))
data <- data%>%select(-End)
data$Start <- as.integer(data$Start)
data$End <- data$Start+100
data <- data[,c("chrom","Start","End","colorectal","lung","prostate","ovarian")]
# add the Type column for the markers
data <- data%>%mutate(Type=case_when((data$colorectal!=0) &(data$lung==0)&(data$prostate==0)~"colorectal",
                                (data$colorectal==0) &(data$lung!=0)&(data$prostate==0)~"lung",
                                (data$colorectal==0) &(data$lung==0)&(data$ovarian!=0)~"ovarian",
                                TRUE~"multiple"))
# add the Color column
data <- data%>%mutate(color=case_when((data$colorectal!=0) &(data$lung==0)&(data$prostate==0)~"6a3d9a",
                                     (data$colorectal==0) &(data$lung!=0)&(data$prostate==0)~"33a02c",
                                     (data$colorectal==0) &(data$lung==0)&(data$ovarian!=0)~"ff7f00",
                                     TRUE~"36ff33"))

# add shape column
data$Shape <- "triangle"
# subset the necessary columns
plot_data <- data[,c("Type","Shape","chrom","Start","End","color")]
names(plot_data) <-c("Type","Shape","Chr","Start","End","color") 
# create the standard dataset
plot_data$Type <- as.character(plot_data$Type)
plot_data$Start <- as.integer(plot_data$Start)
plot_data$End <- as.integer(plot_data$End)
plot_data$Chr <- substr(plot_data$Chr, 4, nchar(plot_data$Chr))

ideogram(human_karyotype, label = as.data.frame(plot_data), label_type = "marker", output = "chromosome.svg")
convertSVG("chromosome.svg", device = "png")

