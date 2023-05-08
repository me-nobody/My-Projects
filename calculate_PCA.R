library(factoextra)
# for dataset before batch correction
before <- as.data.frame(t(qntl_fpkm_merged))
# pca
before.pca <- prcomp(before)
# screeplot
fviz_eig(before.pca)
# plot the PCA
fviz_pca_ind(before.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )+ggtitle("samples-PCA")
# plot correlation of genes
fviz_pca_var(before.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )+ggtitle("genes")
# for dataset after batch correction
after <- as.data.frame(t(corrected_df))
#pca
after.pca <- prcomp(after)
# screeplot
fviz_eig(after.pca)
# plot the PCA
fviz_pca_ind(after.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)+ggtitle("samples-PCA")
