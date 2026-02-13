#############################################################################
#
#	Clust 0 vs All
#
#############################################################################

library(devtools)
library(monocle3)
library(Seurat)
library(SeuratWrappers)  # For conversion helper
library(scater)
library(TSCAN)
library(slingshot)
library(tradeSeq)
library(scran)
library(ggplot2)
library(RColorBrewer)
library(grDevices)
library(DT)
library(pheatmap)




# DEA prep
#Idents(seurat_obj) <- "clust_0.5"  # Set the identity class

# Perform DEA between two groups
#de_markers <- FindMarkers(seurat_obj, ident.1 = "1", ident.2 = "0", test.use = "MAST")


##### Read in data
countMatrices <- readRDS("~/applications/scstudio_Jacek/dea/tokens/heart_organoid_S1_3_epicardium_EPDCs_res_0.5/countMatrices.rds")

metadata <- readRDS("~/applications/scstudio_Jacek/dea/tokens/heart_organoid_S1_3_epicardium_EPDCs_res_0.5/metadata.rds")

dimred <- readRDS("~/applications/scstudio_Jacek/dea/tokens/heart_organoid_S1_3_epicardium_EPDCs_res_0.5/dimred.rds")

hvgs <- readRDS("~/applications/scstudio_Jacek/dea/tokens/heart_organoid_S1_3_epicardium_EPDCs_res_0.5/hvgs.rds")

de_markers <- readRDS("~/applications/scstudio_Jacek/dea/tokens/heart_organoid_S1_3_epicardium_EPDCs_res_0.5/dea.rds")


names(countMatrices)

dim(countMatrices[[1]])
dim(countMatrices$rawCountMatrix)

counts <- countMatrices$rawCountMatrix
data.normalized <- countMatrices$normalization

names(dimred)

pca <- dimred[[1]]$pca$x
dim(pca)

tsne <- dimred$tsne$tsne
dim(tsne)

umap <- dimred$umap$umap
dim(umap)

names(metadata)

clust_res <- "clust_res_0.5"
clust <- metadata[, clust_res]
length(clust)

percentage_mt_genes <- metadata$percentage_mt_genes
length(percentage_mt_genes)

doublets_score <- metadata$doublets_score
length(doublets_score)

dataset <- metadata$dataset
length(dataset)

total_features <- metadata$total_features
length(total_features)

library_size <- metadata$library_size
length(library_size)

library_size_normalization <- metadata$library_size_normalization
length(library_size_normalization)

names(hvgs)
top_hvgs <- hvgs$top_hvgs
length(top_hvgs)


# 1. Create the Seurat object from counts

seurat_obj <- CreateSeuratObject(
  counts = data.normalized,
  assay = "RNA",
  names.field = 1L,
  names.delim = "_",
  meta.data = NULL,
  project = "S1_3_epicardium_EPDCs_res_0.5",
)

# Make sure dimnames match the Seurat object
all(rownames(counts) == rownames(seurat_obj))
all(colnames(counts) == colnames(seurat_obj))


# Set the "data" slot for the "RNA" assay
seurat_obj[["RNA"]] <- SetAssayData(seurat_obj[["RNA"]], slot = "data", new.data = data.normalized)



# Add Clusters, Dataset etc to Seurat Object

# Add the cluster information as metadata
seurat_obj$dataset <- dataset
seurat_obj$clust <- clust
seurat_obj$percentage_mt_genes <- percentage_mt_genes
seurat_obj$doublets_score <- doublets_score
seurat_obj$total_features <- total_features
seurat_obj$library_size <- library_size
seurat_obj$library_size_normalization <- library_size_normalization



# Change results directory
output_dir <- "~/applications/scstudio_Jacek/dea/tokens/heart_organoid_S1_3_epicardium_EPDCs_res_0.5/dea"

if (file.exists(output_dir)){
    setwd(file.path(output_dir))
} else {
    dir.create(file.path(output_dir))
    setwd(file.path(output_dir))
    
}



# Get top 20 DE genes
top_genes <- unique(c(de_markers$dea_clust_0.5_0_vs_all_MAST$Gene[ 1:20 ], 
de_markers$dea_clust_0.5_1_vs_all_MAST$Gene[ 1:20 ], 
de_markers$dea_clust_0.5_2_vs_all_MAST$Gene[ 1:20 ], 
de_markers$dea_clust_0.5_3_vs_all_MAST$Gene[ 1:20 ], 
de_markers$dea_clust_0.5_4_vs_all_MAST$Gene[ 1:20 ], 
de_markers$dea_clust_0.5_5_vs_all_MAST$Gene[ 1:20 ], 
de_markers$dea_clust_0.5_6_vs_all_MAST$Gene[ 1:20 ]))

length(top_genes)


top_genes_selected <- unique(c(de_markers$dea_clust_0.5_0_vs_all_MAST$Gene[ 1:20 ], 
de_markers$dea_clust_0.5_1_vs_all_MAST$Gene[ 1:20 ], 
de_markers$dea_clust_0.5_2_vs_all_MAST$Gene[ 1:20 ], 
de_markers$dea_clust_0.5_3_vs_all_MAST$Gene[ 1:20 ]))


# Remove MALAT1 from visualization since it appears ubiquitously expressed across cells and dominates the expression on the plots. It may dilute color scales in heatmaps because its expression range is orders of magnitude higher than other genes
top_genes <- setdiff(top_genes, "MALAT1")

length(top_genes)


# Extract scaled expression data
heatmap_data <- FetchData(seurat_obj, vars = top_genes)


# Ensure scaling only happens on genes with non-zero variance
heatmap_data <- heatmap_data[apply(heatmap_data, 1, sd) > 0, ]

# Scale and remove rows with NA/Inf after scaling
heatmap_matrix <- t(scale(t(heatmap_data)))

# Remove rows with NA/Inf
heatmap_matrix <- t(heatmap_matrix[complete.cases(heatmap_matrix), ])



# Annotate cells
annotation <- seurat_obj@meta.data[, c("clust", "percentage_mt_genes", "library_size", "nCount_RNA", "nFeature_RNA", "doublets_score")]

# Make sure that data and annotation are of the same length
heatmap_matrix <- heatmap_matrix[ , colnames(heatmap_matrix) %in% rownames(annotation) , drop = FALSE ]
annotation <- annotation[ rownames(annotation) %in% colnames(heatmap_matrix), , drop = FALSE ]

# Define color palette (blue-white-red)
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)


# Define custom colors for one annotation
ann_colors <- list(
  clust = c("clust_0" = "brown1", "clust_1" = "dodgerblue", "clust_2" = "grey", "clust_3" = "gold", "clust_4" = "palegreen3", "clust_5" = "peru", "clust_6" = "bisque")
)


png("dea_MAST_top_20_heatmap.png", width = 3500, height = 3500, pointsize = 8, units = "px", res = 300)
pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 8,
         color = heatmap_colors,
         scale = "none")

invisible(dev.off())


# Now only include clusters 0,1,2 and 3
annotation_selected <- annotation[ annotation$clust %in% c("clust_0", "clust_1", "clust_2", "clust_3") ,  ]
heatmap_matrix_selected <- heatmap_matrix[ , colnames(heatmap_matrix) %in% rownames(annotation_selected), drop = FALSE  ]

heatmap_matrix_selected <- heatmap_matrix_selected[ rownames(heatmap_matrix_selected) %in% top_genes_selected,   ]


# Define custom colors for one annotation
ann_colors_selected <- list(
  clust = c("clust_0" = "brown1", "clust_1" = "dodgerblue", "clust_2" = "grey", "clust_3" = "gold")
)

png("dea_MAST_top_20_heatmap_seleted_clust.png", width = 3500, height = 3000, pointsize = 8, units = "px", res = 300)
pheatmap(heatmap_matrix_selected,
         annotation_col = annotation_selected,
         annotation_colors = ann_colors_selected,
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 8,
         color = heatmap_colors,
         scale = "none")

invisible(dev.off())



##### Differenatial expression analysis results visualisation and tables
for (i in 1:length(de_markers) ) {
  
  # First, remove entries with P-values = 1
  de_markers_filtered <- de_markers[[i]][ de_markers[[i]]$"Adjusted p-value" < 1 , ]

  #####  P-value histograms
  png(file = paste0(names(de_markers)[i], "_P_hist.png", sep=""), width = 6, height = 4, units = "in", res = 300, pointsize = 10)

  # Remove P-values = 1
  histogram <- hist(de_markers_filtered$"Adjusted p-value", breaks=seq(0,1,by= 0.01), main="", xlab="Adjusted P-value")
  exprected_p.value <- mean(histogram$counts)+(sd(histogram$counts)*1)
  abline(v=0.05,col="red")
  abline(h=exprected_p.value,col="blue")
  invisible(dev.off())
  
  
  #  Volcano plots of log2 fold-changes versus significance (p-values, +label top 10)
  png(file = paste0(names(de_markers)[i], "_volcano_plot.png", sep = ""), width = 6, height = 4, units = "in", res = 300, pointsize = 10)
  plot(de_markers[[i]]$"Log2FC (mean)",-log10(de_markers[[i]]$"Adjusted p-value"),pch=16,cex=0.5,xlab="Average log2 fold-change",ylab="-log10(adjusted P-value)",main="",col="grey")
  
  #  Highlight genes with logFC above specified threshold
  points(de_markers[[i]][abs(de_markers[[i]]$"Log2FC (mean)")>1,"Log2FC (mean)"],-log10(de_markers[[i]][abs(de_markers[[i]]$"Log2FC (mean)")>1,"Adjusted p-value"]),cex=0.5,pch=16)
  abline(h=-log10(0.05),col="red", lty = 2)
  
  #  Label top 10 most significant genes
  ord <- order(-log10(de_markers[[i]]$"Adjusted p-value"),decreasing=TRUE)
  top10 <- ord[1:10]
  text(de_markers[[i]][top10,"Log2FC (mean)"],-log10(de_markers[[i]][top10,"Adjusted p-value"]), labels=de_markers[[i]]$Gene[top10],cex=0.6,col="blue")
  invisible(dev.off())

  # Write DE results into a files
  write.table(de_markers[[i]], file = paste0(names(de_markers)[i], "_DE_stats.txt", sep=""), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA, append = FALSE)
}





#############################################################################
#
#	S1 vs S3 (endothelial_cells)
#
#############################################################################

library(devtools)
library(monocle3)
library(Seurat)
library(SeuratWrappers)  # For conversion helper
library(scater)
library(TSCAN)
library(slingshot)
library(tradeSeq)
library(scran)
library(ggplot2)
library(RColorBrewer)
library(grDevices)
library(DT)
library(pheatmap)




# DEA prep
#Idents(seurat_obj) <- "clust_0.5"  # Set the identity class

# Perform DEA between two groups
#de_markers <- FindMarkers(seurat_obj, ident.1 = "1", ident.2 = "0", test.use = "MAST")


##### Read in data
countMatrices <- readRDS("~/applications/scstudio_Jacek/dea/tokens/heart_organoid_S1_3_endothelial_cells/countMatrices.rds")

metadata <- readRDS("~/applications/scstudio_Jacek/dea/tokens/heart_organoid_S1_3_endothelial_cells/metadata.rds")

dimred <- readRDS("~/applications/scstudio_Jacek/dea/tokens/heart_organoid_S1_3_endothelial_cells/dimred.rds")

hvgs <- readRDS("~/applications/scstudio_Jacek/dea/tokens/heart_organoid_S1_3_endothelial_cells/hvgs.rds")

de_markers <- readRDS("~/applications/scstudio_Jacek/dea/tokens/heart_organoid_S1_3_endothelial_cells/dea.rds")


names(countMatrices)

dim(countMatrices[[1]])
dim(countMatrices$rawCountMatrix)

counts <- countMatrices$rawCountMatrix
data.normalized <- countMatrices$normalization

names(dimred)

pca <- dimred[[1]]$pca$x
dim(pca)

tsne <- dimred$tsne$tsne
dim(tsne)

umap <- dimred$umap$umap
dim(umap)

names(metadata)

clust_res <- "clust_res_0.1"
clust <- metadata[, clust_res]
length(clust)

percentage_mt_genes <- metadata$percentage_mt_genes
length(percentage_mt_genes)

doublets_score <- metadata$doublets_score
length(doublets_score)

dataset <- metadata$dataset
length(dataset)

total_features <- metadata$total_features
length(total_features)

library_size <- metadata$library_size
length(library_size)

library_size_normalization <- metadata$library_size_normalization
length(library_size_normalization)

names(hvgs)
top_hvgs <- hvgs$top_hvgs
length(top_hvgs)


# 1. Create the Seurat object from counts

seurat_obj <- CreateSeuratObject(
  counts = data.normalized,
  assay = "RNA",
  names.field = 1L,
  names.delim = "_",
  meta.data = NULL,
  project = "S1_3_epicardium_EPDCs_res_0.5",
)

# Make sure dimnames match the Seurat object
all(rownames(counts) == rownames(seurat_obj))
all(colnames(counts) == colnames(seurat_obj))


# Set the "data" slot for the "RNA" assay
seurat_obj[["RNA"]] <- SetAssayData(seurat_obj[["RNA"]], slot = "data", new.data = data.normalized)



# Add Clusters, Dataset etc to Seurat Object

# Add the cluster information as metadata
seurat_obj$dataset <- dataset
seurat_obj$clust <- clust
seurat_obj$percentage_mt_genes <- percentage_mt_genes
seurat_obj$doublets_score <- doublets_score
seurat_obj$total_features <- total_features
seurat_obj$library_size <- library_size
seurat_obj$library_size_normalization <- library_size_normalization



# Change results directory
output_dir <- "~/applications/scstudio_Jacek/dea/tokens/heart_organoid_S1_3_endothelial_cells/dea"

if (file.exists(output_dir)){
    setwd(file.path(output_dir))
} else {
    dir.create(file.path(output_dir))
    setwd(file.path(output_dir))
    
}



# Get top 20 DE genes
top_genes <- unique(c(de_markers$dea_S1_vs_S3_MAST$Gene[ 1:20 ]))

length(top_genes)


# Remove MALAT1 from visualization since it appears ubiquitously expressed across cells and dominates the expression on the plots. It may dilute color scales in heatmaps because its expression range is orders of magnitude higher than other genes
top_genes <- setdiff(top_genes, "MALAT1")

length(top_genes)


# Extract scaled expression data
heatmap_data <- FetchData(seurat_obj, vars = top_genes)


# Ensure scaling only happens on genes with non-zero variance
heatmap_data <- heatmap_data[apply(heatmap_data, 1, sd) > 0, ]

# Scale and remove rows with NA/Inf after scaling
heatmap_matrix <- t(scale(t(heatmap_data)))

# Remove rows with NA/Inf
heatmap_matrix <- t(heatmap_matrix[complete.cases(heatmap_matrix), ])



# Annotate cells
annotation <- seurat_obj@meta.data[, c("dataset", "clust", "percentage_mt_genes", "library_size", "nCount_RNA", "nFeature_RNA", "doublets_score")]

# Make sure that data and annotation are of the same length
heatmap_matrix <- heatmap_matrix[ , colnames(heatmap_matrix) %in% rownames(annotation) , drop = FALSE ]
annotation <- annotation[ rownames(annotation) %in% colnames(heatmap_matrix), , drop = FALSE ]

# Define color palette (blue-white-red)
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)


# Define custom colors for one annotation
ann_colors <- list(
  clust = c("0" = "brown1", "1" = "dodgerblue", "2" = "grey"),
  dataset = c("S3" = "brown1", "S1" = "dodgerblue")
)


png("dea_MAST_top_20_heatmap.png", width = 3500, height = 2300, pointsize = 8, units = "px", res = 300)
pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 8,
         color = heatmap_colors,
         scale = "none")

invisible(dev.off())


##### Differenatial expression analysis results visualisation and tables
for (i in 1:length(de_markers) ) {
  
  # First, remove entries with P-values = 1
  de_markers_filtered <- de_markers[[i]][ de_markers[[i]]$"Adjusted p-value" < 1 , ]

  #####  P-value histograms
  png(file = paste0(names(de_markers)[i], "_P_hist.png", sep=""), width = 6, height = 4, units = "in", res = 300, pointsize = 10)

  # Remove P-values = 1
  histogram <- hist(de_markers_filtered$"Adjusted p-value", breaks=seq(0,1,by= 0.01), main="", xlab="Adjusted P-value")
  exprected_p.value <- mean(histogram$counts)+(sd(histogram$counts)*1)
  abline(v=0.05,col="red")
  abline(h=exprected_p.value,col="blue")
  invisible(dev.off())
  
  
  #  Volcano plots of log2 fold-changes versus significance (p-values, +label top 10)
  png(file = paste0(names(de_markers)[i], "_volcano_plot.png", sep = ""), width = 6, height = 4, units = "in", res = 300, pointsize = 10)
  plot(de_markers[[i]]$"Log2FC (mean)",-log10(de_markers[[i]]$"Adjusted p-value"),pch=16,cex=0.5,xlab="Average log2 fold-change",ylab="-log10(adjusted P-value)",main="",col="grey")
  
  #  Highlight genes with logFC above specified threshold
  points(de_markers[[i]][abs(de_markers[[i]]$"Log2FC (mean)")>1,"Log2FC (mean)"],-log10(de_markers[[i]][abs(de_markers[[i]]$"Log2FC (mean)")>1,"Adjusted p-value"]),cex=0.5,pch=16)
  abline(h=-log10(0.05),col="red", lty = 2)
  
  #  Label top 10 most significant genes
  ord <- order(-log10(de_markers[[i]]$"Adjusted p-value"),decreasing=TRUE)
  top10 <- ord[1:10]
  text(de_markers[[i]][top10,"Log2FC (mean)"],-log10(de_markers[[i]][top10,"Adjusted p-value"]), labels=de_markers[[i]]$Gene[top10],cex=0.6,col="blue")
  invisible(dev.off())

  # Write DE results into a files
  write.table(de_markers[[i]], file = paste0(names(de_markers)[i], "_DE_stats.txt", sep=""), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA, append = FALSE)
}






#############################################################################
#
#	S1 vs S3 (fibroblasts)
#
#############################################################################

library(devtools)
library(monocle3)
library(Seurat)
library(SeuratWrappers)  # For conversion helper
library(scater)
library(TSCAN)
library(slingshot)
library(tradeSeq)
library(scran)
library(ggplot2)
library(RColorBrewer)
library(grDevices)
library(DT)
library(pheatmap)




# DEA prep
#Idents(seurat_obj) <- "clust_0.5"  # Set the identity class

# Perform DEA between two groups
#de_markers <- FindMarkers(seurat_obj, ident.1 = "1", ident.2 = "0", test.use = "MAST")


##### Read in data
countMatrices <- readRDS("~/applications/scstudio_Jacek/dea/tokens/heart_organoid_S1_3_fibroblasts/countMatrices.rds")

metadata <- readRDS("~/applications/scstudio_Jacek/dea/tokens/heart_organoid_S1_3_fibroblasts/metadata.rds")

dimred <- readRDS("~/applications/scstudio_Jacek/dea/tokens/heart_organoid_S1_3_fibroblasts/dimred.rds")

hvgs <- readRDS("~/applications/scstudio_Jacek/dea/tokens/heart_organoid_S1_3_fibroblasts/hvgs.rds")

de_markers <- readRDS("~/applications/scstudio_Jacek/dea/tokens/heart_organoid_S1_3_fibroblasts/dea.rds")


names(countMatrices)

dim(countMatrices[[1]])
dim(countMatrices$rawCountMatrix)

counts <- countMatrices$rawCountMatrix
data.normalized <- countMatrices$normalization

names(dimred)

pca <- dimred[[1]]$pca$x
dim(pca)

tsne <- dimred$tsne$tsne
dim(tsne)

umap <- dimred$umap$umap
dim(umap)

names(metadata)

clust_res <- "clust_res_0.2"
clust <- metadata[, clust_res]
length(clust)

percentage_mt_genes <- metadata$percentage_mt_genes
length(percentage_mt_genes)

doublets_score <- metadata$doublets_score
length(doublets_score)

dataset <- metadata$dataset
length(dataset)

total_features <- metadata$total_features
length(total_features)

library_size <- metadata$library_size
length(library_size)

library_size_normalization <- metadata$library_size_normalization
length(library_size_normalization)

names(hvgs)
top_hvgs <- hvgs$top_hvgs
length(top_hvgs)


# 1. Create the Seurat object from counts

seurat_obj <- CreateSeuratObject(
  counts = data.normalized,
  assay = "RNA",
  names.field = 1L,
  names.delim = "_",
  meta.data = NULL,
  project = "S1_3_epicardium_EPDCs_res_0.5",
)

# Make sure dimnames match the Seurat object
all(rownames(counts) == rownames(seurat_obj))
all(colnames(counts) == colnames(seurat_obj))


# Set the "data" slot for the "RNA" assay
seurat_obj[["RNA"]] <- SetAssayData(seurat_obj[["RNA"]], slot = "data", new.data = data.normalized)



# Add Clusters, Dataset etc to Seurat Object

# Add the cluster information as metadata
seurat_obj$dataset <- dataset
seurat_obj$clust <- clust
seurat_obj$percentage_mt_genes <- percentage_mt_genes
seurat_obj$doublets_score <- doublets_score
seurat_obj$total_features <- total_features
seurat_obj$library_size <- library_size
seurat_obj$library_size_normalization <- library_size_normalization



# Change results directory
output_dir <- "~/applications/scstudio_Jacek/dea/tokens/heart_organoid_S1_3_fibroblasts/dea"

if (file.exists(output_dir)){
    setwd(file.path(output_dir))
} else {
    dir.create(file.path(output_dir))
    setwd(file.path(output_dir))
    
}



# Get top 20 DE genes
top_genes <- unique(c(de_markers$dea_S1_vs_S3_MAST$Gene[ 1:20 ]))

length(top_genes)


# Remove MALAT1 from visualization since it appears ubiquitously expressed across cells and dominates the expression on the plots. It may dilute color scales in heatmaps because its expression range is orders of magnitude higher than other genes
top_genes <- setdiff(top_genes, "MALAT1")

length(top_genes)


# Extract scaled expression data
heatmap_data <- FetchData(seurat_obj, vars = top_genes)


# Ensure scaling only happens on genes with non-zero variance
heatmap_data <- heatmap_data[apply(heatmap_data, 1, sd) > 0, ]

# Scale and remove rows with NA/Inf after scaling
heatmap_matrix <- t(scale(t(heatmap_data)))

# Remove rows with NA/Inf
heatmap_matrix <- t(heatmap_matrix[complete.cases(heatmap_matrix), ])



# Annotate cells
annotation <- seurat_obj@meta.data[, c("dataset", "clust", "percentage_mt_genes", "library_size", "nCount_RNA", "nFeature_RNA", "doublets_score")]

# Make sure that data and annotation are of the same length
heatmap_matrix <- heatmap_matrix[ , colnames(heatmap_matrix) %in% rownames(annotation) , drop = FALSE ]
annotation <- annotation[ rownames(annotation) %in% colnames(heatmap_matrix), , drop = FALSE ]

# Define color palette (blue-white-red)
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)


# Define custom colors for one annotation
ann_colors <- list(
  clust = c("0" = "brown1", "1" = "dodgerblue", "2" = "grey", "3" = "gold"),
  dataset = c("S3" = "brown1", "S1" = "dodgerblue")
)


png("dea_MAST_top_20_heatmap.png", width = 3500, height = 2300, pointsize = 8, units = "px", res = 300)
pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 8,
         color = heatmap_colors,
         scale = "none")

invisible(dev.off())


##### Differenatial expression analysis results visualisation and tables
for (i in 1:length(de_markers) ) {
  
  # First, remove entries with P-values = 1
  de_markers_filtered <- de_markers[[i]][ de_markers[[i]]$"Adjusted p-value" < 1 , ]

  #####  P-value histograms
  png(file = paste0(names(de_markers)[i], "_P_hist.png", sep=""), width = 6, height = 4, units = "in", res = 300, pointsize = 10)

  # Remove P-values = 1
  histogram <- hist(de_markers_filtered$"Adjusted p-value", breaks=seq(0,1,by= 0.01), main="", xlab="Adjusted P-value")
  exprected_p.value <- mean(histogram$counts)+(sd(histogram$counts)*1)
  abline(v=0.05,col="red")
  abline(h=exprected_p.value,col="blue")
  invisible(dev.off())
  
  
  #  Volcano plots of log2 fold-changes versus significance (p-values, +label top 10)
  png(file = paste0(names(de_markers)[i], "_volcano_plot.png", sep = ""), width = 6, height = 4, units = "in", res = 300, pointsize = 10)
  plot(de_markers[[i]]$"Log2FC (mean)",-log10(de_markers[[i]]$"Adjusted p-value"),pch=16,cex=0.5,xlab="Average log2 fold-change",ylab="-log10(adjusted P-value)",main="",col="grey")
  
  #  Highlight genes with logFC above specified threshold
  points(de_markers[[i]][abs(de_markers[[i]]$"Log2FC (mean)")>1,"Log2FC (mean)"],-log10(de_markers[[i]][abs(de_markers[[i]]$"Log2FC (mean)")>1,"Adjusted p-value"]),cex=0.5,pch=16)
  abline(h=-log10(0.05),col="red", lty = 2)
  
  #  Label top 10 most significant genes
  ord <- order(-log10(de_markers[[i]]$"Adjusted p-value"),decreasing=TRUE)
  top10 <- ord[1:10]
  text(de_markers[[i]][top10,"Log2FC (mean)"],-log10(de_markers[[i]][top10,"Adjusted p-value"]), labels=de_markers[[i]]$Gene[top10],cex=0.6,col="blue")
  invisible(dev.off())

  # Write DE results into a files
  write.table(de_markers[[i]], file = paste0(names(de_markers)[i], "_DE_stats.txt", sep=""), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA, append = FALSE)
}






#############################################################################
#
#	S1 vs S3 (mural cells)
#
#############################################################################

library(devtools)
library(monocle3)
library(Seurat)
library(SeuratWrappers)  # For conversion helper
library(scater)
library(TSCAN)
library(slingshot)
library(tradeSeq)
library(scran)
library(ggplot2)
library(RColorBrewer)
library(grDevices)
library(DT)
library(pheatmap)




# DEA prep
#Idents(seurat_obj) <- "clust_0.5"  # Set the identity class

# Perform DEA between two groups
#de_markers <- FindMarkers(seurat_obj, ident.1 = "1", ident.2 = "0", test.use = "MAST")


##### Read in data
countMatrices <- readRDS("~/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_mural_cells/countMatrices.rds")

metadata <- readRDS("~/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_mural_cells/metadata.rds")

dimred <- readRDS("~/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_mural_cells/dimred.rds")

hvgs <- readRDS("~/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_mural_cells/hvgs.rds")

de_markers <- readRDS("~/applications/scstudio_Jacek/dea/tokens/heart_organoid_S1_3_mural_cells/dea.rds")


names(countMatrices)

dim(countMatrices[[1]])
dim(countMatrices$rawCountMatrix)

counts <- countMatrices$rawCountMatrix
data.normalized <- countMatrices$normalization

names(dimred)

pca <- dimred[[1]]$pca$x
dim(pca)

tsne <- dimred$tsne$tsne
dim(tsne)

umap <- dimred$umap$umap
dim(umap)

names(metadata)

clust_res <- "clust_res_0.1"
clust <- metadata[, clust_res]
length(clust)

percentage_mt_genes <- metadata$percentage_mt_genes
length(percentage_mt_genes)

doublets_score <- metadata$doublets_score
length(doublets_score)

dataset <- metadata$dataset
length(dataset)

total_features <- metadata$total_features
length(total_features)

library_size <- metadata$library_size
length(library_size)

library_size_normalization <- metadata$library_size_normalization
length(library_size_normalization)

names(hvgs)
top_hvgs <- hvgs$top_hvgs
length(top_hvgs)


# 1. Create the Seurat object from counts

seurat_obj <- CreateSeuratObject(
  counts = data.normalized,
  assay = "RNA",
  names.field = 1L,
  names.delim = "_",
  meta.data = NULL,
  project = "S1_3_epicardium_EPDCs_res_0.5",
)

# Make sure dimnames match the Seurat object
all(rownames(counts) == rownames(seurat_obj))
all(colnames(counts) == colnames(seurat_obj))


# Set the "data" slot for the "RNA" assay
seurat_obj[["RNA"]] <- SetAssayData(seurat_obj[["RNA"]], slot = "data", new.data = data.normalized)



# Add Clusters, Dataset etc to Seurat Object

# Add the cluster information as metadata
seurat_obj$dataset <- dataset
seurat_obj$clust <- clust
seurat_obj$percentage_mt_genes <- percentage_mt_genes
seurat_obj$doublets_score <- doublets_score
seurat_obj$total_features <- total_features
seurat_obj$library_size <- library_size
seurat_obj$library_size_normalization <- library_size_normalization



# Change results directory
output_dir <- "~/applications/scstudio_Jacek/dea/tokens/heart_organoid_S1_3_mural_cells/dea"

if (file.exists(output_dir)){
    setwd(file.path(output_dir))
} else {
    dir.create(file.path(output_dir))
    setwd(file.path(output_dir))
    
}



# Get top 20 DE genes
top_genes <- unique(c(de_markers$dea_S1_vs_S3_MAST$Gene[ 1:20 ]))

length(top_genes)


# Remove MALAT1 from visualization since it appears ubiquitously expressed across cells and dominates the expression on the plots. It may dilute color scales in heatmaps because its expression range is orders of magnitude higher than other genes
top_genes <- setdiff(top_genes, "MALAT1")

length(top_genes)


# Extract scaled expression data
heatmap_data <- FetchData(seurat_obj, vars = top_genes)


# Ensure scaling only happens on genes with non-zero variance
heatmap_data <- heatmap_data[apply(heatmap_data, 1, sd) > 0, ]

# Scale and remove rows with NA/Inf after scaling
heatmap_matrix <- t(scale(t(heatmap_data)))

# Remove rows with NA/Inf
heatmap_matrix <- t(heatmap_matrix[complete.cases(heatmap_matrix), ])



# Annotate cells
annotation <- seurat_obj@meta.data[, c("dataset", "clust", "percentage_mt_genes", "library_size", "nCount_RNA", "nFeature_RNA", "doublets_score")]

# Make sure that data and annotation are of the same length
heatmap_matrix <- heatmap_matrix[ , colnames(heatmap_matrix) %in% rownames(annotation) , drop = FALSE ]
annotation <- annotation[ rownames(annotation) %in% colnames(heatmap_matrix), , drop = FALSE ]

# Define color palette (blue-white-red)
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)


# Define custom colors for one annotation
ann_colors <- list(
  clust = c("0" = "brown1", "1" = "dodgerblue", "2" = "grey"),
  dataset = c("S3" = "brown1", "S1" = "dodgerblue")
)


png("dea_MAST_top_20_heatmap.png", width = 3500, height = 2300, pointsize = 8, units = "px", res = 300)
pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 8,
         color = heatmap_colors,
         scale = "none")

invisible(dev.off())


##### Differenatial expression analysis results visualisation and tables
for (i in 1:length(de_markers) ) {
  
  # First, remove entries with P-values = 1
  de_markers_filtered <- de_markers[[i]][ de_markers[[i]]$"Adjusted p-value" < 1 , ]

  #####  P-value histograms
  png(file = paste0(names(de_markers)[i], "_P_hist.png", sep=""), width = 6, height = 4, units = "in", res = 300, pointsize = 10)

  # Remove P-values = 1
  histogram <- hist(de_markers_filtered$"Adjusted p-value", breaks=seq(0,1,by= 0.01), main="", xlab="Adjusted P-value")
  exprected_p.value <- mean(histogram$counts)+(sd(histogram$counts)*1)
  abline(v=0.05,col="red")
  abline(h=exprected_p.value,col="blue")
  invisible(dev.off())
  
  
  #  Volcano plots of log2 fold-changes versus significance (p-values, +label top 10)
  png(file = paste0(names(de_markers)[i], "_volcano_plot.png", sep = ""), width = 6, height = 4, units = "in", res = 300, pointsize = 10)
  plot(de_markers[[i]]$"Log2FC (mean)",-log10(de_markers[[i]]$"Adjusted p-value"),pch=16,cex=0.5,xlab="Average log2 fold-change",ylab="-log10(adjusted P-value)",main="",col="grey")
  
  #  Highlight genes with logFC above specified threshold
  points(de_markers[[i]][abs(de_markers[[i]]$"Log2FC (mean)")>1,"Log2FC (mean)"],-log10(de_markers[[i]][abs(de_markers[[i]]$"Log2FC (mean)")>1,"Adjusted p-value"]),cex=0.5,pch=16)
  abline(h=-log10(0.05),col="red", lty = 2)
  
  #  Label top 10 most significant genes
  ord <- order(-log10(de_markers[[i]]$"Adjusted p-value"),decreasing=TRUE)
  top10 <- ord[1:10]
  text(de_markers[[i]][top10,"Log2FC (mean)"],-log10(de_markers[[i]][top10,"Adjusted p-value"]), labels=de_markers[[i]]$Gene[top10],cex=0.6,col="blue")
  invisible(dev.off())

  # Write DE results into a files
  write.table(de_markers[[i]], file = paste0(names(de_markers)[i], "_DE_stats.txt", sep=""), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA, append = FALSE)
}





#############################################################################
#
#	S1 vs S3 (epicardium / EPDCs)
#
#############################################################################

library(devtools)
library(monocle3)
library(Seurat)
library(SeuratWrappers)  # For conversion helper
library(scater)
library(TSCAN)
library(slingshot)
library(tradeSeq)
library(scran)
library(ggplot2)
library(RColorBrewer)
library(grDevices)
library(DT)
library(pheatmap)




# DEA prep
#Idents(seurat_obj) <- "clust_0.5"  # Set the identity class

# Perform DEA between two groups
#de_markers <- FindMarkers(seurat_obj, ident.1 = "1", ident.2 = "0", test.use = "MAST")


##### Read in data
countMatrices <- readRDS("~/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_epicardium_EPDCs/countMatrices.rds")

metadata <- readRDS("~/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_epicardium_EPDCs/metadata.rds")

dimred <- readRDS("~/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_epicardium_EPDCs/dimred.rds")

hvgs <- readRDS("~/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_epicardium_EPDCs/hvgs.rds")

de_markers <- readRDS("~/applications/scstudio_Jacek/dea/tokens/heart_organoid_S1_3_epicardium_EPDCs/dea.rds")


names(countMatrices)

dim(countMatrices[[1]])
dim(countMatrices$rawCountMatrix)

counts <- countMatrices$rawCountMatrix
data.normalized <- countMatrices$normalization

names(dimred)

pca <- dimred[[1]]$pca$x
dim(pca)

tsne <- dimred$tsne$tsne
dim(tsne)

umap <- dimred$umap$umap
dim(umap)

names(metadata)

clust_res <- "clust_res_0.5"
clust <- metadata[, clust_res]
length(clust)

percentage_mt_genes <- metadata$percentage_mt_genes
length(percentage_mt_genes)

doublets_score <- metadata$doublets_score
length(doublets_score)

dataset <- metadata$dataset
length(dataset)

total_features <- metadata$total_features
length(total_features)

library_size <- metadata$library_size
length(library_size)

library_size_normalization <- metadata$library_size_normalization
length(library_size_normalization)

names(hvgs)
top_hvgs <- hvgs$top_hvgs
length(top_hvgs)


# 1. Create the Seurat object from counts

seurat_obj <- CreateSeuratObject(
  counts = data.normalized,
  assay = "RNA",
  names.field = 1L,
  names.delim = "_",
  meta.data = NULL,
  project = "S1_3_epicardium_EPDCs_res_0.5",
)

# Make sure dimnames match the Seurat object
all(rownames(counts) == rownames(seurat_obj))
all(colnames(counts) == colnames(seurat_obj))


# Set the "data" slot for the "RNA" assay
seurat_obj[["RNA"]] <- SetAssayData(seurat_obj[["RNA"]], slot = "data", new.data = data.normalized)



# Add Clusters, Dataset etc to Seurat Object

# Add the cluster information as metadata
seurat_obj$dataset <- dataset
seurat_obj$clust <- clust
seurat_obj$percentage_mt_genes <- percentage_mt_genes
seurat_obj$doublets_score <- doublets_score
seurat_obj$total_features <- total_features
seurat_obj$library_size <- library_size
seurat_obj$library_size_normalization <- library_size_normalization



# Change results directory
output_dir <- "~/applications/scstudio_Jacek/dea/tokens/heart_organoid_S1_3_epicardium_EPDCs/dea"

if (file.exists(output_dir)){
    setwd(file.path(output_dir))
} else {
    dir.create(file.path(output_dir))
    setwd(file.path(output_dir))
    
}



# Get top 20 DE genes
top_genes <- unique(c(de_markers$dea_S1_vs_S3_MAST$Gene[ 1:20 ]))

length(top_genes)


# Remove MALAT1 from visualization since it appears ubiquitously expressed across cells and dominates the expression on the plots. It may dilute color scales in heatmaps because its expression range is orders of magnitude higher than other genes
top_genes <- setdiff(top_genes, "MALAT1")

length(top_genes)


# Extract scaled expression data
heatmap_data <- FetchData(seurat_obj, vars = top_genes)


# Ensure scaling only happens on genes with non-zero variance
heatmap_data <- heatmap_data[apply(heatmap_data, 1, sd) > 0, ]

# Scale and remove rows with NA/Inf after scaling
heatmap_matrix <- t(scale(t(heatmap_data)))

# Remove rows with NA/Inf
heatmap_matrix <- t(heatmap_matrix[complete.cases(heatmap_matrix), ])



# Annotate cells
annotation <- seurat_obj@meta.data[, c("dataset", "clust", "percentage_mt_genes", "library_size", "nCount_RNA", "nFeature_RNA", "doublets_score")]

# Make sure that data and annotation are of the same length
heatmap_matrix <- heatmap_matrix[ , colnames(heatmap_matrix) %in% rownames(annotation) , drop = FALSE ]
annotation <- annotation[ rownames(annotation) %in% colnames(heatmap_matrix), , drop = FALSE ]

# Define color palette (blue-white-red)
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)


# Define custom colors for one annotation
ann_colors <- list(
  clust = c("0" = "brown1", "1" = "dodgerblue", "2" = "grey", "3" = "gold", "4" = "palegreen3", "5" = "peru", "6" = "bisque"),
  dataset = c("S3" = "brown1", "S1" = "dodgerblue")
)


png("dea_MAST_top_20_heatmap.png", width = 3500, height = 2300, pointsize = 8, units = "px", res = 300)
pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 8,
         color = heatmap_colors,
         scale = "none")

invisible(dev.off())


##### Differenatial expression analysis results visualisation and tables
for (i in 1:length(de_markers) ) {
  
  # First, remove entries with P-values = 1
  de_markers_filtered <- de_markers[[i]][ de_markers[[i]]$"Adjusted p-value" < 1 , ]

  #####  P-value histograms
  png(file = paste0(names(de_markers)[i], "_P_hist.png", sep=""), width = 6, height = 4, units = "in", res = 300, pointsize = 10)

  # Remove P-values = 1
  histogram <- hist(de_markers_filtered$"Adjusted p-value", breaks=seq(0,1,by= 0.01), main="", xlab="Adjusted P-value")
  exprected_p.value <- mean(histogram$counts)+(sd(histogram$counts)*1)
  abline(v=0.05,col="red")
  abline(h=exprected_p.value,col="blue")
  invisible(dev.off())
  
  
  #  Volcano plots of log2 fold-changes versus significance (p-values, +label top 10)
  png(file = paste0(names(de_markers)[i], "_volcano_plot.png", sep = ""), width = 6, height = 4, units = "in", res = 300, pointsize = 10)
  plot(de_markers[[i]]$"Log2FC (mean)",-log10(de_markers[[i]]$"Adjusted p-value"),pch=16,cex=0.5,xlab="Average log2 fold-change",ylab="-log10(adjusted P-value)",main="",col="grey")
  
  #  Highlight genes with logFC above specified threshold
  points(de_markers[[i]][abs(de_markers[[i]]$"Log2FC (mean)")>1,"Log2FC (mean)"],-log10(de_markers[[i]][abs(de_markers[[i]]$"Log2FC (mean)")>1,"Adjusted p-value"]),cex=0.5,pch=16)
  abline(h=-log10(0.05),col="red", lty = 2)
  
  #  Label top 10 most significant genes
  ord <- order(-log10(de_markers[[i]]$"Adjusted p-value"),decreasing=TRUE)
  top10 <- ord[1:10]
  text(de_markers[[i]][top10,"Log2FC (mean)"],-log10(de_markers[[i]][top10,"Adjusted p-value"]), labels=de_markers[[i]]$Gene[top10],cex=0.6,col="blue")
  invisible(dev.off())

  # Write DE results into a files
  write.table(de_markers[[i]], file = paste0(names(de_markers)[i], "_DE_stats.txt", sep=""), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA, append = FALSE)
}





