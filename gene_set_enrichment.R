library(AnnotationHub)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
library(GenomicRanges)
library(annotatr)
library(Homo.sapiens)
# create the mart object
ensembl <-useEnsembl(biomart="genes", dataset = "hsapiens_gene_ensembl")
# get the list of important genes
enriched_genes <- df_dm_annotated_genes$annot.symbol
enriched_genes <- unique(enriched_genes)
# get the full list of genes
binned_peaks <- merge(binned_scc104_narrow,binned_scc152_narrow,by=c("chr","start","end"))
# name the columns 
names(binned_peaks) <- c("chrom","chromStart","chromEnd","scc104","scc152")
# create a GRanges object from the binned dataframe
bins.GRanges_all <- makeGRangesFromDataFrame(binned_peaks)
# merge the binned dataframe with TxDb annotation
annotated.bins <- subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), bins.GRanges_all)
# list the set of annotations needed
annots_genes="hg19_genes_cds"
# Build the annotations (a single GRanges object)
annotations = build_annotations(genome = 'hg19', annotations = annots)
# Intersect the regions we read in with the annotations
all_bins_annotated = annotate_regions(
  regions = bins.GRanges_all,
  annotations = annotations_genes, # only genes annotated
  ignore.strand = TRUE,
  quiet = FALSE)
# Coerce to a data.frame
all_bins_annotated_df <-  data.frame(all_bins_annotated)
all_bins_annotated_df <- all_bins_annotated_df[complete.cases(all_bins_annotated_df),]
# get list of all genes in the ChIP peak
all_genes <- all_bins_annotated_df$annot.symbol
all_genes <- unique(all_genes)
# call biomart function to interconvert all gene names to entrez ids
gene_2_id_universe<-getBM(attributes = c('external_gene_name','entrezgene_id'),
                          filters = 'external_gene_name',
                          values = all_genes, # here all_data_gene_name is list of gene names
                          mart = ensembl)
gene_2_id_universe<-gene_2_id_universe[complete.cases(gene_2_id_universe),]
# extract gene_ids for GO enrichment
all_gene_ids<-gene_2_id_universe$entrezgene_id # here all_genes is list of entrez ids
# call the function to interconvert genes to entrez ids of filtered genes
gene_2_id_filtered<-getBM(attributes = c('external_gene_name','entrezgene_id'),
                          filters = 'external_gene_name',
                          values = enriched_genes, 
                          mart = ensembl)
gene_2_id_filtered<-gene_2_id_filtered[complete.cases(gene_2_id_filtered),]
# extract sample gene ids for GO enrichment
filtered_gene_ids<-gene_2_id_filtered$entrezgene_id
# create enrichGO object
ego_chip <- enrichGO(gene = filtered_gene_ids,
                           universe      = all_gene_ids,
                           OrgDb         = org.Hs.eg.db,
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 0.05,
                           readable      = TRUE)
# view the enrichGO dataframe
dim(ego_chip)
head(ego_chip[,c(1,2)],20)
# plot the data
barplot(ego_chip,showCategory = 16,x="p.adjust",ylab="pathways")
dotplot(ego_chip,showCategory=16,x="p.adjust") + ggtitle("overrepresented pathways in HPV infection")
