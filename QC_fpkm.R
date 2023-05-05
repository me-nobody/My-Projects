library(dplyr)
library(DESeq2)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(sva)
par(mfrow=c(3,1))
tumor <- sapply(fpkm_gse156451,sum)
cellline <- sapply(fpkm_gse185055,sum)
cetuximab <- sapply(fpkm_gse196576,sum)
hist(tumor,main="fpkm counts tumor",cex=1.2,cex.lab=1.8,xlab="fpkm counts")
hist(cellline,main="fpkm counts cell line",cex=1.2,cex.lab=1.8,xlab="fpkm counts")
hist(cetuximab,main="fpkm counts cetuximab",cex=1.2,cex.lab=1.8,xlab="fpkm counts")
dev.off()
# we need to obtain a common subset of fpkm reads 
a <- rownames(fpkm_gse156451)
b <- rownames(fpkm_gse185055)
c <- rownames(fpkm_gse196576)

common_gene_subset <- intersect(intersect(a,b),c)
# subset all dataframes with the common subset of genes
fpkm_gse156451 <- fpkm_gse156451[rownames(fpkm_gse156451)%in%common_gene_subset,]
fpkm_gse185055 <- fpkm_gse185055[rownames(fpkm_gse185055)%in%common_gene_subset,]
fpkm_gse196576 <- fpkm_gse196576[rownames(fpkm_gse196576)%in%common_gene_subset,]

# create a common merged dataframe
merged_fpkm <- cbind(fpkm_gse156451,fpkm_gse185055,fpkm_gse196576)

# plot PCA to study batch effects
sample_names <- names(merged_fpkm)
condition <- c(rep("tumor",50),rep("cell_line",50),rep("cetuximab",50))
colData <- as.data.frame(cbind(sample_names,condition))
colData$condition <- as.factor(colData$condition)
dds <- DESeqDataSetFromMatrix(countData = merged_fpkm,colData = colData,design = ~ 1)

plotPCA(merged_fpkm) # won't work with ratios

# using ComBat from sva package to correct for batch effects
# create a batch
batch <- c(rep(1,50),rep(2,50),rep(3,50))
# run the ComBat function
batch_corr_merged_fpkm = ComBat(dat=merged_fpkm,
                          batch=batch,
                          mod=NULL,
                          par.prior=TRUE,
                          prior.plots=TRUE)

# final dataframe for analysis obtained!



