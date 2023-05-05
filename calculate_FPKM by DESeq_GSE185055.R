library(dplyr)
library(DESeq2)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# we shall use DESeq2 to calculate FPKM
#--------calculate for GSE185055------------------------------------------------

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
# filter the low count genes to remove noise
keep_gse185055 <- rowSums(counts(dds_gse185055))>=10
dds_55 <- dds_gse185055[keep_gse185055,]
# the dataframe gene_name_gene_lengths was obtained from ENSEMBL biomart.
# we have to derive the basepairs column to provide as input for FPKM
gene_lengths <- gene_name_gene_lengths
names(gene_lengths) <- c("genes","start","end")
gene_lengths$basepairs <- gene_lengths$end-gene_lengths$start

# get the list and order of genes in the dds object
dds_genes_gse185055 <- rownames(dds_55)
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
# we will assign these gene lengths to GSE185055
mcols(dds_55)$basepairs <- dds_genes_gse185055$basepairs
# now calculation of FPKM by DESeq2
fpkm_gse185055 <- as.data.frame(fpkm(dds_55))







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
#----------calculate for GSE196576---------------------------------
