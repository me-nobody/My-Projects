library(dplyr)
library(DESeq2)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
par(mfrow=c(3,1))
tumor <- sapply(GSE156451,sum)
cellline <- sapply(GSE185055,sum)
cetuximab <- sapply(GSE196576,sum)
hist(sapply(tumor,log2),main="raw counts tumor",cex=1.2,cex.lab=1.8,xlab="log2(counts)")
hist(sapply(cellline,log2),main="raw counts cell line",cex=1.2,cex.lab=1.8,xlab="log2(counts)")
hist(sapply(cetuximab,log2),main="raw counts cetuximab",cex=1.2,cex.lab=1.8,xlab="log2(counts)")
dev.off()

# we have again lost the rownames, so proceedign w/o them
rownames(GSE156451) <- paste("tumor",1:50,sep = "_")
rownames(GSE185055) <- paste("cellline",1:50,sep = "_")
rownames(GSE196576) <- paste("cetuximab",1:50,sep = "_")

# we can't combine right away as GSE156451 is FPKM while GSE185055,
# and GSE196576 are count data
# we shall use DESeq2 to calculate FPKM

# first transpose GSE185055
GSE185055 <- as.data.frame(t(GSE185055))
sampleNames <- names(GSE185055)
condition <- as.factor(rep("celline",50))
colData <- as.data.frame(cbind(sampleNames,condition))
colData$sampleNames==names(GSE185055)
# create the deseq object
dds_gse185055 <- DESeqDataSetFromMatrix(countData = GSE185055,colData = colData, design = ~1 )
# here design=~1 has to be used as this is not a differential expression issue.
dds_gse185055 <- DESeq(dds_gse185055)
# get the list of known genes. txdb is obscure. i have simply gone to bioMart website and downloaded
# the necessary table
gene_lengths <- as.data.frame(cbind(gene_names_transcript_lengths$Gene.name,
                                    gene_names_transcript_lengths$Transcript.length..including.UTRs.and.CDS.))
gene_lengths$basepairs <- as.integer(gene_lengths$basepairs)
# above will not work as it is transcripts and there are many transripts
# now getting gene lengths imported from biomart
# issue with duplicate gene names
# due to multiple transcripts, it is very difficult to get row numbers of dds and row numbers 
# of gene_lengths equal. I had to use dplyr distinct
gene_lengths <- gene_name_gene_lengths
names(gene_lengths) <- c("genes","start","end")
gene_lengths$basepairs <- gene_lengths$end-gene_lengths$start
# use dplyr filter to get the genes common
gene_lengths_gse185055 <- gene_lengths%>%filter(genes%in%dds_genes_gse185055)
# there are duplicates and hence use dplyr distinct to get unique entries
gene_lengths_gse185055 <- gene_lengths_gse185055%>%distinct(genes,.keep_all = TRUE)
# now we have to reorder the gene names as per the dds object to get the gene lengths
# first create a dataframe from dds_genes_gse185055
dds_genes_gse185055 <- as.data.frame(dds_genes_gse185055)
names(dds_genes_gse185055) <- "genes"
# now use dplyr left_join to merge gene_lengths_gse185055 with this dataframe
dds_genes_gse185055 <- dds_genes_gse185055%>%left_join(gene_lengths_gse185055,by="genes")
# now from the above dataframe we  will get the genelengths to calculate FPKM
mcols(dds_gse185055)$basepairs <- dds_genes_gse185055$basepairs
# now calculation of FPKM by DESeq2
fpkm_gse185055 <- as.data.frame(fpkm(dds_gse185055))
