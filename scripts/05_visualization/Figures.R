
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
library(reshape2)
library(dplyr)
library(tidyr)
library(ggridges)
library(ggrepel)



##### Set working dir to MS folder

setwd("~/data/heart_organoids/scRNA-seq/docs/MS/Figures")


#############################################################################
#
#	Filtered data (S1 + 3, all cells)
#
#############################################################################


##### Set working dir to MS folder

setwd("~/data/heart_organoids/scRNA-seq/docs/MS/Figures")


##### Read in data

countMatrices <- readRDS("/mnt/scratch/home/jacek/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_annot/countMatrices.rds")

metadata <- readRDS("/mnt/scratch/home/jacek/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_annot/metadata.rds")

dimred <- readRDS("/mnt/scratch/home/jacek/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_annot/dimred.rds")

hvgs <- readRDS("/mnt/scratch/home/jacek/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_annot/hvgs.rds")



names(countMatrices)

cell_type <- metadata$cell_type
length(cell_type)

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

clust <- metadata$clust_res_0.1
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
  counts = counts,
  assay = "RNA",
  names.field = 1L,
  names.delim = "_",
  meta.data = NULL,
  project = "S1_3",
)


# Make sure dimnames match the Seurat object
all(rownames(counts) == rownames(seurat_obj))
all(colnames(counts) == colnames(seurat_obj))



# Set the "data" slot for the "RNA" assay
seurat_obj[["RNA"]] <- SetAssayData(seurat_obj[["RNA"]], slot = "data", new.data = data.normalized)



# Download barcode whitelist

whitelist_1 = read.table(gzfile("/mnt/scratch/home/jacek/data/heart_organoids/scRNA-seq/Sample1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"),sep="\t")
whitelist_2 = read.table(gzfile("/mnt/scratch/home/jacek/data/heart_organoids/scRNA-seq/Sample2/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"),sep="\t")
whitelist_3 = read.table(gzfile("/mnt/scratch/home/jacek/data/heart_organoids/scRNA-seq/Sample3/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"),sep="\t")

# Combine all barcodes lists
whitelist <- rbind(whitelist_1, whitelist_2, whitelist_3)

# Number of cells
n_cells <- ncol(seurat_obj)

# Assign to cells in your matrix
valid_barcodes <- unique(unlist(whitelist))[1:n_cells]
colnames(seurat_obj) <- valid_barcodes


# 2. Add variable features
VariableFeatures(seurat_obj) <- top_hvgs


# 3. Add dimensionality reductions

# Add PCA
# Make sure it is a numeric matrix with rownames = cell names

pca <- as.matrix(pca)
rownames(pca) <- colnames(seurat_obj)  # Ensure rows are cells

seurat_obj[["pca"]] <- Seurat::CreateDimReducObject(
  embeddings = pca,
  key = "PC_",
  assay = DefaultAssay(seurat_obj)
)


# Add UMAP

umap <- as.matrix(umap)

# Set rownames of umap to match Seurat cell names (i.e., colnames of the count matrix)
rownames(umap) <- colnames(seurat_obj)

seurat_obj[["umap"]] <- Seurat::CreateDimReducObject(
  embeddings = umap,
  key = "UMAP_",
  assay = DefaultAssay(seurat_obj)
)


# Add t-SNE

tsne <- as.matrix(tsne)

# Set rownames of tsne to match Seurat cell names (i.e., colnames of the count matrix)
rownames(tsne) <- colnames(seurat_obj)

# Now you can safely create and add the tsne slot
seurat_obj[["tsne"]] <- Seurat::CreateDimReducObject(
  embeddings = tsne,
  key = "tSNE_",
  assay = DefaultAssay(seurat_obj)
)


# Add PCA
#seurat_obj[["pca"]] <- CreateDimReducObject(embeddings = as.matrix(pca), key = "PC_", assay = DefaultAssay(seurat_obj))

# Add t-SNE
#seurat_obj[["tsne"]] <- CreateDimReducObject(embeddings = as.matrix(tsne), key = "tSNE_", assay = DefaultAssay(seurat_obj))

# Add UMAP
#seurat_obj[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(umap), key = "UMAP_", assay = DefaultAssay(seurat_obj))



# 4. Add Clusters, Dataset etc to Seurat Object

# Add the cluster information as metadata
seurat_obj$dataset <- dataset
seurat_obj$clust <- clust
seurat_obj$cell_type <- cell_type
seurat_obj$percentage_mt_genes <- percentage_mt_genes
seurat_obj$doublets_score <- doublets_score
seurat_obj$total_features <- total_features
seurat_obj$library_size <- library_size
seurat_obj$library_size_normalization <- library_size_normalization


head(seurat_obj)



# Update cell type names

annotation <- seurat_obj@meta.data

annotation[ annotation$cell_type == "Mural_cells", c("cell_type") ] <- "Mural cells"
annotation[ annotation$cell_type == "Epicardium_EPDCs", c("cell_type") ] <- "Epicardium/EPDCs"
annotation[ annotation$cell_type == "Proliferative_fibroblasts", c("cell_type") ] <- "Proliferative fibroblasts"
annotation[ annotation$cell_type == "Endothelial_cells", c("cell_type") ] <- "Endothelial cells"
annotation[ annotation$cell_type == "Endoderm-derived_epithelial_cells", c("cell_type") ] <- "Endoderm-derived epithelial cells"
annotation[ annotation$cell_type == "Hepatic-derived_cells", c("cell_type") ] <- "Hepatic-derived cells"

annotation[ annotation$dataset == "S1", c("dataset") ] <- "VEGF"
annotation[ annotation$dataset == "S3", c("dataset") ] <- "PDGFBB"

seurat_obj <- AddMetaData(seurat_obj, annotation$cell_type, col.name ="cell_type")
seurat_obj <- AddMetaData(seurat_obj, annotation$dataset, col.name ="condition")


# Define desired order of cell types
desired_order <- c("Epicardium/EPDCs" , "Fibroblasts", "Proliferative fibroblasts", "Mural cells", "Endothelial cells", "Cardiomyocytes", "Endoderm-derived epithelial cells", "Hepatic-derived cells")


# Get cell (column) names ordered by cell_type
cell_order <- annotation %>%
  dplyr::mutate(cell = rownames(.)) %>%
  dplyr::arrange(factor(cell_type, levels = desired_order)) %>%
  pull(cell)


# Convert cell_type to a factor with specified order
seurat_obj$cell_type <- factor(seurat_obj$cell_type, levels = desired_order)


#############################################################################
#	t-SNE plot (all cells)
#############################################################################

# Basic DimPlot using existing colors
p <- DimPlot(seurat_obj,
             reduction = "tsne",
             group.by = "cell_type",
             label = FALSE,                 # No labels on plot
             pt.size = 0.5) +               # Adjust point size
  theme_classic() +
  labs(x = "t-SNE 1", y = "t-SNE 2") +     # <-- Add axis titles
  theme(legend.position = "right",         # Legend on the side
        #axis.line = element_blank(),
        #axis.ticks = element_blank(),
        #axis.text = element_blank(),
        axis.title = element_blank())      # Clean up axis


# Preserve Current Colors
current_colors <- Seurat::DiscretePalette(length(unique(seurat_obj$cell_type)))
names(current_colors) <- levels(seurat_obj$cell_type)

cell_type_colors <- c("yellow2", "mediumpurple2", "darkgoldenrod3", "darkolivegreen3", "lightpink2", "coral4", "burlywood4", "midnightblue")
names(cell_type_colors) <- c("Epicardium/EPDCs", "Fibroblasts", "Proliferative fibroblasts", "Mural cells", "Endothelial cells", "Cardiomyocytes", "Endoderm-derived epithelial cells", "Hepatic-derived cells")


p <- p + scale_color_manual(values = cell_type_colors)

# Legend Improvements
p <- p + guides(color = guide_legend(
                  title = "Cell Type",      # Legend title
                  override.aes = list(size = 3), # Legend dot size
                  ncol = 1))                # One column (vertical)

p <- p + coord_fixed()   # Ensures equal scaling of x and y axes

p <- p + ggtitle("t-SNE of all cell types")

p <- p + labs(x = "t-SNE 1", y = "t-SNE 2")


# Save the plot at 300–600 dpi (journal-ready):
ggsave("tSNE_all_cell_types.png", p, width = 7, height = 5, dpi = 600)  # PNG (raster)
ggsave("tSNE_all_cell_types.pdf", p, width = 7, height = 5)      




#############################################################################
#	Heatmap for selected genes (all cells)
#############################################################################

genes <- c("TNNT2", "PECAM1", "TTR", "WT1", "TBX18", "SEMA3D", "POSTN", "PTN", "DCN", "TCF21", "TOP2A", "PDGFRB", "RGS5", "MGP", "ACTA2", "ACTN1", "CDH5")


# Extract scaled expression data
heatmap_data <- FetchData(seurat_obj, vars = genes)

# Ensure scaling only happens on genes with non-zero variance
heatmap_data <- heatmap_data[apply(heatmap_data, 1, sd) > 0, ]

# Scale and remove rows with NA/Inf after scaling
heatmap_matrix <- t(scale(t(heatmap_data)))

# Remove rows with NA/Inf
heatmap_matrix <- t(heatmap_matrix[complete.cases(heatmap_matrix), ])



# Annotate cells
annotation <- seurat_obj@meta.data[, c("cell_type", "condition", "total_features", "library_size", "percentage_mt_genes")]

colnames(annotation) <- c("cell_type", "condition", "total_features", "library_size", "percentage_mt_genes")

# Make sure that data and annotation are of the same length
heatmap_matrix <- heatmap_matrix[ , colnames(heatmap_matrix) %in% rownames(annotation) , drop = FALSE ]
annotation <- annotation[ rownames(annotation) %in% colnames(heatmap_matrix), , drop = FALSE ]


# Reorder heatmap matrix columns and annotation rows
heatmap_matrix <- heatmap_matrix[, cell_order]
annotation <- annotation[cell_order, ]

# Define color palette
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)


# Define custom colors for one annotation
ann_colors <- list(
  cell_type = cell_type_colors,
  condition = c("VEGF" = "gold", "PDGFBB" = "firebrick")
)



p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "none")

png("heatmap_all_cell_types_blue_red_annot.png", width = 3500, height = 2000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())


# Annotate cells
annotation <- seurat_obj@meta.data[, c("cell_type", "condition")]
colnames(annotation) <- c("cell_type", "condition")

p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "none")

png("heatmap_all_cell_types_blue_red.png", width = 3500, height = 1500, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())



# Define color palette
heatmap_colors <- viridis::viridis(100)

p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "none")

png("heatmap_all_cell_types_blue_yellow.png", width = 3500, height = 1500, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())




# Define color palette
heatmap_colors <- colorRampPalette(c("royalblue4", "white", "red4"))(100)

p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "none")

png("heatmap_all_cell_types_blue_red_2.png", width = 3500, height = 1500, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())



#############################################################################
#	Dotplot for selected genes (all cells)
############################################################################

# Reverse the y-axis order (top-to-bottom)
seurat_obj$cell_type <- factor(seurat_obj$cell_type, levels = rev(desired_order))

# Generate dot plot grouped by cell types
p <- DotPlot(seurat_obj,
        features = genes,
        group.by = "cell_type",
        scale = FALSE) +
  scale_color_gradient(low = "lightgrey", high = "blue") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_all_cell_types.png", p, width = 10, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_all_cell_types.pdf", p, width = 10, height = 4)  



# Generate dot plot grouped by cell types
p <- DotPlot(seurat_obj,
        features = genes,
        group.by = "cell_type",
        scale = FALSE) +
  scale_color_gradient(low = "lightgrey", high = "red") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_all_cell_types_grey_red.png", p, width = 10, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_all_cell_types_grey_red.pdf", p, width = 10, height = 4)  


# Generate dot plot grouped by cell types
p <- DotPlot(seurat_obj,
        features = genes,
        group.by = "cell_type",
        scale = FALSE) +
  scale_color_gradient(low = "blue", high = "red") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_all_cell_types_blue_red.png", p, width = 10, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_all_cell_types_blue_red.pdf", p, width = 10, height = 4)  



# Generate dot plot grouped by cell types
p <- DotPlot(seurat_obj,
        features = genes,
        group.by = "cell_type",
        scale = FALSE) +
  scale_color_gradient(low = "royalblue", high = "red3") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_all_cell_types_blue_red_2.png", p, width = 10, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_all_cell_types_blue_red_2.pdf", p, width = 10, height = 4)  



# Generate dot plot grouped by cell types
dotplot_colors <- viridis::viridis(2)
p <- DotPlot(seurat_obj,
        features = genes,
        group.by = "cell_type",
        scale = FALSE) +
  scale_color_gradient(low = dotplot_colors[1], high = dotplot_colors[2]) +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_all_cell_types_blue_yellow.png", p, width = 10, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_all_cell_types_blue_yellow.pdf", p, width = 10, height = 4)  



# Generate dot plot grouped by cell types
p <- DotPlot(seurat_obj,
        features = genes,
        group.by = "cell_type",
        scale = FALSE) +
  scale_color_gradient(low = "lightgrey", high = "darkgoldenrod2") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_all_cell_types_darkgoldenrod.png", p, width = 10, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_all_cell_types_darkgoldenrod.pdf", p, width = 10, height = 4)  


# Generate dot plot grouped by cell types
p <- DotPlot(seurat_obj,
        features = genes,
        group.by = "cell_type",
        scale = FALSE) +
  scale_color_gradient(low = "khaki", high = "darkred") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_all_cell_types_yellow_darkred.png", p, width = 10, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_all_cell_types_yellow_darkred.pdf", p, width = 10, height = 4) 



# Reverse the y-axis order again
seurat_obj$cell_type <- factor(seurat_obj$cell_type, levels = desired_order)




#############################################################################
#
#	Filtered data (S1 + 3, fibroblasts)
#
#############################################################################


##### Set working dir to MS folder

setwd("~/data/heart_organoids/scRNA-seq/docs/MS/Figures")


##### Read in data

countMatrices <- readRDS("/mnt/scratch/home/jacek/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_fibroblasts/countMatrices.rds")

metadata <- readRDS("/mnt/scratch/home/jacek/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_fibroblasts/metadata.rds")

dimred <- readRDS("/mnt/scratch/home/jacek/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_fibroblasts/dimred.rds")

hvgs <- readRDS("/mnt/scratch/home/jacek/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_fibroblasts/hvgs.rds")


names(countMatrices)

cell_type <- metadata$cell_type
length(cell_type)

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

clust <- metadata$clust_res_0.2
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
  counts = counts,
  assay = "RNA",
  names.field = 1L,
  names.delim = "_",
  meta.data = NULL,
  project = "S1_3",
)


# Make sure dimnames match the Seurat object
all(rownames(counts) == rownames(seurat_obj))
all(colnames(counts) == colnames(seurat_obj))



# Set the "data" slot for the "RNA" assay
seurat_obj[["RNA"]] <- SetAssayData(seurat_obj[["RNA"]], slot = "data", new.data = data.normalized)



# Download barcode whitelist

whitelist_1 = read.table(gzfile("/mnt/scratch/home/jacek/data/heart_organoids/scRNA-seq/Sample1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"),sep="\t")
whitelist_2 = read.table(gzfile("/mnt/scratch/home/jacek/data/heart_organoids/scRNA-seq/Sample2/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"),sep="\t")
whitelist_3 = read.table(gzfile("/mnt/scratch/home/jacek/data/heart_organoids/scRNA-seq/Sample3/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"),sep="\t")

# Combine all barcodes lists
whitelist <- rbind(whitelist_1, whitelist_2, whitelist_3)

# Number of cells
n_cells <- ncol(seurat_obj)

# Assign to cells in your matrix
valid_barcodes <- unique(unlist(whitelist))[1:n_cells]
colnames(seurat_obj) <- valid_barcodes


# 2. Add variable features
VariableFeatures(seurat_obj) <- top_hvgs


# 3. Add dimensionality reductions

# Add PCA
# Make sure it is a numeric matrix with rownames = cell names

pca <- as.matrix(pca)
rownames(pca) <- colnames(seurat_obj)  # Ensure rows are cells

seurat_obj[["pca"]] <- Seurat::CreateDimReducObject(
  embeddings = pca,
  key = "PC_",
  assay = DefaultAssay(seurat_obj)
)


# Add UMAP

umap <- as.matrix(umap)

# Set rownames of umap to match Seurat cell names (i.e., colnames of the count matrix)
rownames(umap) <- colnames(seurat_obj)

seurat_obj[["umap"]] <- Seurat::CreateDimReducObject(
  embeddings = umap,
  key = "UMAP_",
  assay = DefaultAssay(seurat_obj)
)


# Add t-SNE

tsne <- as.matrix(tsne)

# Set rownames of tsne to match Seurat cell names (i.e., colnames of the count matrix)
rownames(tsne) <- colnames(seurat_obj)

# Now you can safely create and add the tsne slot
seurat_obj[["tsne"]] <- Seurat::CreateDimReducObject(
  embeddings = tsne,
  key = "tSNE_",
  assay = DefaultAssay(seurat_obj)
)


# Add PCA
#seurat_obj[["pca"]] <- CreateDimReducObject(embeddings = as.matrix(pca), key = "PC_", assay = DefaultAssay(seurat_obj))

# Add t-SNE
#seurat_obj[["tsne"]] <- CreateDimReducObject(embeddings = as.matrix(tsne), key = "tSNE_", assay = DefaultAssay(seurat_obj))

# Add UMAP
#seurat_obj[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(umap), key = "UMAP_", assay = DefaultAssay(seurat_obj))



# 4. Add Clusters, Dataset etc to Seurat Object

# Add the cluster information as metadata
seurat_obj$dataset <- dataset
seurat_obj$clust <- clust
seurat_obj$cell_type <- cell_type
seurat_obj$percentage_mt_genes <- percentage_mt_genes
seurat_obj$doublets_score <- doublets_score
seurat_obj$total_features <- total_features
seurat_obj$library_size <- library_size
seurat_obj$library_size_normalization <- library_size_normalization



head(seurat_obj)


# Update cell type names

annotation <- seurat_obj@meta.data

annotation[ annotation$dataset == "S1", c("dataset") ] <- "VEGF"
annotation[ annotation$dataset == "S3", c("dataset") ] <- "PDGFBB"

seurat_obj <- AddMetaData(seurat_obj, annotation$dataset, col.name ="condition")
seurat_obj$condition <- factor(seurat_obj$condition, levels = c("VEGF", "PDGFBB"))


# Define desired order of cell types
desired_order <- c(0,1,2,3)


# Get cell (column) names ordered by cluster
cell_order <- annotation %>%
  dplyr::mutate(cell = rownames(.)) %>%
  dplyr::arrange(factor(clust, levels = desired_order)) %>%
  pull(cell)


# Convert clust to a factor with specified order
seurat_obj$clust <- factor(seurat_obj$clust, levels = desired_order)


#############################################################################
#	t-SNE plot (fibroblasts)
#############################################################################

# Basic DimPlot using existing colors
p <- DimPlot(seurat_obj,
             reduction = "tsne",
             group.by = "clust",
             label = FALSE,                 # No labels on plot
             pt.size = 0.5) +               # Adjust point size
  theme_classic() +
  labs(x = "t-SNE 1", y = "t-SNE 2") +     # <-- Add axis titles
  theme(legend.position = "right",         # Legend on the side
        #axis.line = element_blank(),
        #axis.ticks = element_blank(),
        #axis.text = element_blank(),
        axis.title = element_blank())      # Clean up axis


# Preserve Current Colors
current_colors <- Seurat::DiscretePalette(length(unique(seurat_obj$clust)))
names(current_colors) <- levels(seurat_obj$clust)

cell_type_colors <- c("darkolivegreen3", "cadetblue", "darkorchid", "saddlebrown")
names(cell_type_colors) <- c("0", "1", "2", "3")


p <- p + scale_color_manual(values = cell_type_colors)

# Legend Improvements
p <- p + guides(color = guide_legend(
                  title = "Cluster",      # Legend title
                  override.aes = list(size = 3), # Legend dot size
                  ncol = 1))                # One column (vertical)

p <- p + coord_fixed()   # Ensures equal scaling of x and y axes

p <- p + ggtitle("t-SNE of fibroblasts")

p <- p + labs(x = "t-SNE 1", y = "t-SNE 2")


# Save the plot at 300–600 dpi (journal-ready):
ggsave("tSNE_fibroblasts.png", p, width = 7, height = 5, dpi = 600)  # PNG (raster)
ggsave("tSNE_fibroblasts.pdf", p, width = 7, height = 5)      


# Basic DimPlot using existing colors (split samples)
p <- DimPlot(seurat_obj,
             reduction = "tsne",
             group.by = "clust",
             split.by = "condition",   # Splits plot by sample,
             label = FALSE,                 # No labels on plot
             pt.size = 0.5) +               # Adjust point size
  theme_classic() +
  labs(x = "t-SNE 1", y = "t-SNE 2") +     # <-- Add axis titles
  theme(legend.position = "right",         # Legend on the side
        #axis.line = element_blank(),
        #axis.ticks = element_blank(),
        #axis.text = element_blank(),
        axis.title = element_blank())      # Clean up axis


# Preserve Current Colors
current_colors <- Seurat::DiscretePalette(length(unique(seurat_obj$clust)))
names(current_colors) <- levels(seurat_obj$clust)


p <- p + scale_color_manual(values = cell_type_colors)

# Legend Improvements
p <- p + guides(color = guide_legend(
                  title = "Cluster",      # Legend title
                  override.aes = list(size = 3), # Legend dot size
                  ncol = 1))                # One column (vertical)

p <- p + coord_fixed()   # Ensures equal scaling of x and y axes

p <- p + ggtitle("t-SNE of fibroblasts")

p <- p + labs(x = "t-SNE 1", y = "t-SNE 2")


# Save the plot at 300–600 dpi (journal-ready):
ggsave("tSNE_fibroblasts_by_sample.png", p, width = 7, height = 5, dpi = 600)  # PNG (raster)
ggsave("tSNE_fibroblasts_by_sample.pdf", p, width = 7, height = 5)      



#############################################################################
#	Heatmap for selected genes (fibroblasts)
#############################################################################

genes <- c("VEGFA", "EFNB2", "ACKR3", "EGR1", "PLXNA4", "ACTA2", "IGFBP7", "LUM", "TGFBI", "BGN", "UACA", "ITGA1", "TPM2", "SEMA3D")


# Extract scaled expression data
heatmap_data <- FetchData(seurat_obj, vars = genes)

# Reorder heatmap matrix columns and annotation rows
heatmap_data <- heatmap_data[ cell_order, ]

# Ensure scaling only happens on genes with non-zero variance
heatmap_data <- heatmap_data[apply(heatmap_data, 1, sd) > 0, ]

# Scale and remove rows with NA/Inf after scaling
heatmap_matrix <- t(scale(t(heatmap_data)))

# Remove rows with NA/Inf
heatmap_matrix <- t(heatmap_matrix[complete.cases(heatmap_matrix), ])



# Annotate cells
annotation <- seurat_obj@meta.data[, c("clust", "condition", "total_features", "library_size", "percentage_mt_genes")]
colnames(annotation) <- c("cluster", "condition", "total_features", "library_size", "percentage_mt_genes")

# Make sure that data and annotation are of the same length
heatmap_matrix <- heatmap_matrix[ , colnames(heatmap_matrix) %in% rownames(annotation) , drop = FALSE ]
annotation <- annotation[ rownames(annotation) %in% colnames(heatmap_matrix), , drop = FALSE ]

# Define color palette
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)


# Define custom colors for one annotation
ann_colors <- list(
  cluster = cell_type_colors,
  condition = c("VEGF" = "gold", "PDGFBB" = "firebrick")
)



p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "none")

png("heatmap_fibroblasts_annot.png", width = 3500, height = 2000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())


# Annotate cells
annotation <- seurat_obj@meta.data[, c("clust", "condition")]
colnames(annotation) <- c("cluster", "condition")

p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "none")

png("heatmap_fibroblasts.png", width = 3500, height = 1500, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())



# Define color palette
heatmap_colors <- viridis::viridis(100)

p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "none")

png("heatmap_fibroblasts_blue_yellow.png", width = 3500, height = 1500, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())




# Define color palette
heatmap_colors <- colorRampPalette(c("royalblue4", "white", "red4"))(100)

p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "none")

png("heatmap_fibroblasts_blue_red_2.png", width = 3500, height = 1500, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())





#############################################################################
#	Violin plots for selected genes (fibroblasts)
#############################################################################


# Violin plot grouped by clusters
p <- VlnPlot(seurat_obj,
             features = genes,
             group.by = "clust",
             pt.size = 0,                    # no points
             cols = cell_type_colors,
             combine = TRUE)                 # one figure with facets


# Adjust facet strip title size
p <- p + theme(strip.text.x = element_text(size = 6))

ggsave("violin_fibroblasts_clusters.png", p, width = 7, height = 10, dpi = 600)   # PNG
ggsave("violin_fibroblasts_clusters.pdf", p, width = 7, height = 10)              # PDF (vector)



# Violin plot grouped by clusters (split by sample)
p <- VlnPlot(seurat_obj,
             features = genes,
             group.by = "clust",
             split.by = "condition",   # Splits plot by sample,
             pt.size = 0,                    # no points
             cols = ann_colors$condition,
             combine = TRUE)                 # one figure with facets


# Adjust facet strip title size
p <- p + theme(strip.text.x = element_text(size = 6))

p <- p + 
  scale_fill_manual(
    values = ann_colors$condition,  # pick colors for samples
    name = "Condition"                                  # legend title
  ) +
  theme(legend.position = "right")   


ggsave("violin_fibroblasts_samples.png", p, width = 10, height = 10, dpi = 600)   # PNG
ggsave("violin_fibroblasts_samples.pdf", p, width = 10, height = 10)              # PDF (vector)





#############################################################################
#
#	Filtered data (S1 + 3, endothelial cells)
#
#############################################################################


##### Set working dir to MS folder

setwd("~/data/heart_organoids/scRNA-seq/docs/MS/Figures")


##### Read in data

countMatrices <- readRDS("/mnt/scratch/home/jacek/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_endothelial_cells/countMatrices.rds")

metadata <- readRDS("/mnt/scratch/home/jacek/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_endothelial_cells/metadata.rds")

dimred <- readRDS("/mnt/scratch/home/jacek/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_endothelial_cells/dimred.rds")

hvgs <- readRDS("/mnt/scratch/home/jacek/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_endothelial_cells/hvgs.rds")


names(countMatrices)

cell_type <- metadata$cell_type
length(cell_type)

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

clust <- metadata$clust_res_0.1
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
  counts = counts,
  assay = "RNA",
  names.field = 1L,
  names.delim = "_",
  meta.data = NULL,
  project = "S1_3",
)


# Make sure dimnames match the Seurat object
all(rownames(counts) == rownames(seurat_obj))
all(colnames(counts) == colnames(seurat_obj))



# Set the "data" slot for the "RNA" assay
seurat_obj[["RNA"]] <- SetAssayData(seurat_obj[["RNA"]], slot = "data", new.data = data.normalized)



# Download barcode whitelist

whitelist_1 = read.table(gzfile("/mnt/scratch/home/jacek/data/heart_organoids/scRNA-seq/Sample1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"),sep="\t")
whitelist_2 = read.table(gzfile("/mnt/scratch/home/jacek/data/heart_organoids/scRNA-seq/Sample2/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"),sep="\t")
whitelist_3 = read.table(gzfile("/mnt/scratch/home/jacek/data/heart_organoids/scRNA-seq/Sample3/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"),sep="\t")

# Combine all barcodes lists
whitelist <- rbind(whitelist_1, whitelist_2, whitelist_3)

# Number of cells
n_cells <- ncol(seurat_obj)

# Assign to cells in your matrix
valid_barcodes <- unique(unlist(whitelist))[1:n_cells]
colnames(seurat_obj) <- valid_barcodes


# 2. Add variable features
VariableFeatures(seurat_obj) <- top_hvgs


# 3. Add dimensionality reductions

# Add PCA
# Make sure it is a numeric matrix with rownames = cell names

pca <- as.matrix(pca)
rownames(pca) <- colnames(seurat_obj)  # Ensure rows are cells

seurat_obj[["pca"]] <- Seurat::CreateDimReducObject(
  embeddings = pca,
  key = "PC_",
  assay = DefaultAssay(seurat_obj)
)


# Add UMAP

umap <- as.matrix(umap)

# Set rownames of umap to match Seurat cell names (i.e., colnames of the count matrix)
rownames(umap) <- colnames(seurat_obj)

seurat_obj[["umap"]] <- Seurat::CreateDimReducObject(
  embeddings = umap,
  key = "UMAP_",
  assay = DefaultAssay(seurat_obj)
)


# Add t-SNE

tsne <- as.matrix(tsne)

# Set rownames of tsne to match Seurat cell names (i.e., colnames of the count matrix)
rownames(tsne) <- colnames(seurat_obj)

# Now you can safely create and add the tsne slot
seurat_obj[["tsne"]] <- Seurat::CreateDimReducObject(
  embeddings = tsne,
  key = "tSNE_",
  assay = DefaultAssay(seurat_obj)
)


# Add PCA
#seurat_obj[["pca"]] <- CreateDimReducObject(embeddings = as.matrix(pca), key = "PC_", assay = DefaultAssay(seurat_obj))

# Add t-SNE
#seurat_obj[["tsne"]] <- CreateDimReducObject(embeddings = as.matrix(tsne), key = "tSNE_", assay = DefaultAssay(seurat_obj))

# Add UMAP
#seurat_obj[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(umap), key = "UMAP_", assay = DefaultAssay(seurat_obj))



# 4. Add Clusters, Dataset etc to Seurat Object

# Add the cluster information as metadata
seurat_obj$dataset <- dataset
seurat_obj$clust <- clust
seurat_obj$cell_type <- cell_type
seurat_obj$percentage_mt_genes <- percentage_mt_genes
seurat_obj$doublets_score <- doublets_score
seurat_obj$total_features <- total_features
seurat_obj$library_size <- library_size
seurat_obj$library_size_normalization <- library_size_normalization



head(seurat_obj)


# Update cell type names

annotation <- seurat_obj@meta.data

annotation[ annotation$dataset == "S1", c("dataset") ] <- "VEGF"
annotation[ annotation$dataset == "S3", c("dataset") ] <- "PDGFBB"

seurat_obj <- AddMetaData(seurat_obj, annotation$dataset, col.name ="condition")
seurat_obj$condition <- factor(seurat_obj$condition, levels = c("VEGF", "PDGFBB"))


# Define desired order of cell types
desired_order <- c(0,1,2)


# Get cell (column) names ordered by cluster
cell_order <- annotation %>%
  dplyr::mutate(cell = rownames(.)) %>%
  dplyr::arrange(factor(clust, levels = desired_order)) %>%
  pull(cell)


# Convert clust to a factor with specified order
seurat_obj$clust <- factor(seurat_obj$clust, levels = desired_order)


#############################################################################
#	t-SNE plot (endothelial cells)
#############################################################################

# Basic DimPlot using existing colors
p <- DimPlot(seurat_obj,
             reduction = "tsne",
             group.by = "clust",
             label = FALSE,                 # No labels on plot
             pt.size = 0.5) +               # Adjust point size
  theme_classic() +
  labs(x = "t-SNE 1", y = "t-SNE 2") +     # <-- Add axis titles
  theme(legend.position = "right",         # Legend on the side
        #axis.line = element_blank(),
        #axis.ticks = element_blank(),
        #axis.text = element_blank(),
        axis.title = element_blank())      # Clean up axis


# Preserve Current Colors
current_colors <- Seurat::DiscretePalette(length(unique(seurat_obj$clust)))
names(current_colors) <- levels(seurat_obj$clust)

cell_type_colors <- c("cornflowerblue", "darkgoldenrod2", "darkolivegreen2")
names(cell_type_colors) <- c("0", "1", "2")


p <- p + scale_color_manual(values = cell_type_colors)

# Legend Improvements
p <- p + guides(color = guide_legend(
                  title = "Cluster",      # Legend title
                  override.aes = list(size = 3), # Legend dot size
                  ncol = 1))                # One column (vertical)

p <- p + coord_fixed()   # Ensures equal scaling of x and y axes

p <- p + ggtitle("t-SNE of endothelial cells")

p <- p + labs(x = "t-SNE 1", y = "t-SNE 2")


# Save the plot at 300–600 dpi (journal-ready):
ggsave("tSNE_endothelial_cells.png", p, width = 7, height = 5, dpi = 600)  # PNG (raster)
ggsave("tSNE_endothelial_cells.pdf", p, width = 7, height = 5)      


# Basic DimPlot using existing colors (split samples)
p <- DimPlot(seurat_obj,
             reduction = "tsne",
             group.by = "clust",
             split.by = "condition",   # Splits plot by sample,
             label = FALSE,                 # No labels on plot
             pt.size = 0.5) +               # Adjust point size
  theme_classic() +
  labs(x = "t-SNE 1", y = "t-SNE 2") +     # <-- Add axis titles
  theme(legend.position = "right",         # Legend on the side
        #axis.line = element_blank(),
        #axis.ticks = element_blank(),
        #axis.text = element_blank(),
        axis.title = element_blank())      # Clean up axis


# Preserve Current Colors
current_colors <- Seurat::DiscretePalette(length(unique(seurat_obj$clust)))
names(current_colors) <- levels(seurat_obj$clust)


p <- p + scale_color_manual(values = cell_type_colors)

# Legend Improvements
p <- p + guides(color = guide_legend(
                  title = "Cluster",      # Legend title
                  override.aes = list(size = 3), # Legend dot size
                  ncol = 1))                # One column (vertical)

p <- p + coord_fixed()   # Ensures equal scaling of x and y axes

p <- p + ggtitle("t-SNE of endothelial cells")

p <- p + labs(x = "t-SNE 1", y = "t-SNE 2")


# Save the plot at 300–600 dpi (journal-ready):
ggsave("tSNE_endothelial_cells_by_sample.png", p, width = 7, height = 5, dpi = 600)  # PNG (raster)
ggsave("tSNE_endothelial_cells_by_sample.pdf", p, width = 7, height = 5)      




#############################################################################
#	t-SNE plots for selected genes (endothelial cells)
#############################################################################


####################
genes <- c("A2M", "FABP4", "CD36", "COL15A1", "IGFBP7", "FABP5", "TNNT3", "CLDN5", "TCF4", "TM4SF18", "PDGFB", "C8orf4", "CAV1", "APLNR", "COL4A1", "RGCC", "CXorf36", "VWF", "JAM2", "UACA", "COL4A2", "UTRN", "HLA-B", "DST", "ADGRL4", "C1QTNF9", "CA4", "GIMAP7", "ROBO4", "CD74", "MEOX1", "MEF2C", "MTSS1", "SDPR", "GIMAP4", "APOD", "FLT4", "GMFG", "RBP7", "ARHGAP18", "ITGA6")
####################

# Generate one plot per gene (faceted automatically if combine=TRUE)
p <- FeaturePlot(seurat_obj,
            features = genes,
            min.cutoff = "q05",
            max.cutoff = "q90", # Uses max.cutoff to avoid oversaturation from extreme outliers
            reduction = "tsne",     # use t-SNE coordinates
            #cols = c("gray90", "lightgray", "darkred"),  # low → high expression
            cols = c("khaki", "khaki", "darkred"),  # low → high expression
            combine = TRUE,         # one panel per gene
            pt.size = 0.01,
            order = TRUE) &       # <-- ensures high-expression plotted last
  theme(
    plot.title = element_text(size=12),
    axis.title = element_text(size=10),
    axis.text  = element_text(size=8),
    legend.text = element_text(size=8),
    legend.title = element_text(size=9)
  )


# Save the plot at 300–600 dpi (journal-ready):
ggsave("t-SNE_endothelial_cells.png", p, width = 14, height = 24, dpi = 600)  # PNG (raster)
ggsave("t-SNE_endothelial_cells.pdf", p, width = 14, height = 24)  



#############################################################################
#	Violin plots for selected genes (endothelial cells)
#############################################################################

############
genes <- c("DLL4", "JAG1", "CXCR4", "EFNB2", "VEGFC", "PLVAP", "FLRT2", "NRG1", "EPB41", "BST2", "SELENOP", "APLNR", "NR2F2", "FLRT3", "GATA4", "CD36",  "NPR3", "NRP2", "UNC5B", "NOTCH4")
############


# Violin plot grouped by clusters
p <- VlnPlot(seurat_obj,
             features = genes,
             group.by = "clust",
             pt.size = 0,                    # no points
             cols = cell_type_colors,
             combine = TRUE,                 # one figure with facets
             ncol = 5)                       # force 5 columns


# Adjust facet strip title size
p <- p + theme(strip.text.x = element_text(size = 6))

ggsave("violin_endothelial_cells_clusters.png", p, width = 9, height = 9, dpi = 600)   # PNG
ggsave("violin_endothelial_cells_clusters.pdf", p, width = 9, height = 9)              # PDF (vector)



# Violin plot grouped by clusters (split by sample)

# Define custom colors for one annotation
ann_colors <- list(
  cluster = cell_type_colors,
  condition = c("VEGF" = "gold", "PDGFBB" = "firebrick")
)


p <- VlnPlot(seurat_obj,
             features = genes,
             group.by = "clust",
             split.by = "condition",   # Splits plot by sample,
             pt.size = 0,                    # no points
             cols = ann_colors$condition,
             combine = TRUE,                 # one figure with facets
             ncol = 5)                       # force 5 columns

# Adjust facet strip title size
p <- p + theme(strip.text.x = element_text(size = 6))

p <- p + 
  scale_fill_manual(
    values = ann_colors$condition,  # pick colors for samples
    name = "Condition"                                  # legend title
  ) +
  theme(legend.position = "right")   


ggsave("violin_endothelial_cells_samples.png", p, width = 11, height = 9, dpi = 600)   # PNG
ggsave("violin_endothelial_cells_samples.pdf", p, width = 11, height = 9)              # PDF (vector)



############
genes <- c("A2M", "FABP4", "CD36", "COL15A1", "IGFBP7", "FABP5", "TNNT3", "CLDN5", "TCF4", "TM4SF18", "PDGFB", "C8orf4", "CAV1", "APLNR", "COL4A1", "RGCC", "CXorf36", "VWF", "JAM2", "UACA", "COL4A2", "UTRN", "HLA-B", "DST", "ADGRL4", "C1QTNF9", "CA4", "GIMAP7", "ROBO4", "CD74", "MEOX1", "MEF2C", "MTSS1", "SDPR", "GIMAP4", "APOD", "FLT4", "GMFG", "RBP7", "ARHGAP18", "ITGA6")
############


# Violin plot grouped by clusters
p <- VlnPlot(seurat_obj,
             features = genes,
             group.by = "clust",
             pt.size = 0,                    # no points
             cols = cell_type_colors,
             combine = TRUE,                 # one figure with facets
             ncol = 5)                       # force 5 columns


# Adjust facet strip title size
p <- p + theme(strip.text.x = element_text(size = 6))

ggsave("violin_endothelial_cells_geneset_2_clusters.png", p, width = 11, height = 18, dpi = 600)   # PNG
ggsave("violin_endothelial_cells_geneset_2_clusters.pdf", p, width = 11, height = 18)              # PDF (vector)



# Violin plot grouped by clusters (split by sample)

# Define custom colors for one annotation
ann_colors <- list(
  cluster = cell_type_colors,
  condition = c("VEGF" = "gold", "PDGFBB" = "firebrick")
)


p <- VlnPlot(seurat_obj,
             features = genes,
             group.by = "clust",
             split.by = "condition",   # Splits plot by sample,
             pt.size = 0,                    # no points
             cols = ann_colors$condition,
             combine = TRUE,                 # one figure with facets
             ncol = 5)                       # force 5 columns

# Adjust facet strip title size
p <- p + theme(strip.text.x = element_text(size = 6))

p <- p + 
  scale_fill_manual(
    values = ann_colors$condition,  # pick colors for samples
    name = "Condition"                                  # legend title
  ) +
  theme(legend.position = "bottom")   


ggsave("violin_endothelial_cells_geneset_2_samples.png", p, width = 15, height = 18, dpi = 600)   # PNG
ggsave("violin_endothelial_cells_geneset_2_samples.pdf", p, width = 15, height = 18)              # PDF (vector)



############
genes <- c("NRP1", "NRP2", "GJA5", "BMX", "NOTCH1", "NOTCH4", "DLL4", "JAG1", "JAG2", "HEY2", "KDR", "TBX20", "ALK", "EPAS1", "DEPP1", "EPHB4", "LEFTY1", "LEFTY2", "FLT4", "TEK", "NR2F2", "EMCN")
############


# Violin plot grouped by clusters
p <- VlnPlot(seurat_obj,
             features = genes,
             group.by = "clust",
             pt.size = 0,                    # no points
             cols = cell_type_colors,
             combine = TRUE,                 # one figure with facets
             ncol = 5)                       # force 5 columns


# Adjust facet strip title size
p <- p + theme(strip.text.x = element_text(size = 6))

ggsave("violin_endothelial_cells_geneset_3_clusters.png", p, width = 11, height = 12, dpi = 600)   # PNG
ggsave("violin_endothelial_cells_geneset_3_clusters.pdf", p, width = 11, height = 12)              # PDF (vector)



# Violin plot grouped by clusters (split by sample)

# Define custom colors for one annotation
ann_colors <- list(
  cluster = cell_type_colors,
  condition = c("VEGF" = "gold", "PDGFBB" = "firebrick")
)


p <- VlnPlot(seurat_obj,
             features = genes,
             group.by = "clust",
             split.by = "condition",   # Splits plot by sample,
             pt.size = 0,                    # no points
             cols = ann_colors$condition,
             combine = TRUE,                 # one figure with facets
             ncol = 5)                       # force 5 columns

# Adjust facet strip title size
p <- p + theme(strip.text.x = element_text(size = 6))

p <- p + 
  scale_fill_manual(
    values = ann_colors$condition,  # pick colors for samples
    name = "Condition"                                  # legend title
  ) +
  theme(legend.position = "bottom")   


ggsave("violin_endothelial_cells_geneset_3_samples.png", p, width = 15, height = 12, dpi = 600)   # PNG
ggsave("violin_endothelial_cells_geneset_3_samples.pdf", p, width = 15, height = 12)              # PDF (vector)




############
genes <- c("CXCR4", "VEGFC", "EFNB2", "FLRT2", "GATA4", "PDGFB")
############


# Violin plot grouped by clusters
p <- VlnPlot(seurat_obj,
             features = genes,
             group.by = "clust",
             pt.size = 0,                    # no points
             cols = cell_type_colors,
             combine = TRUE,                 # one figure with facets
             ncol = 3)                       # force 3 columns


# Adjust facet strip title size
p <- p + theme(strip.text.x = element_text(size = 6))

ggsave("violin_endothelial_cells_geneset_4_clusters.png", p, width = 7, height = 5, dpi = 600)   # PNG
ggsave("violin_endothelial_cells_geneset_4_clusters.pdf", p, width = 7, height = 5)              # PDF (vector)



# Violin plot grouped by clusters (split by sample)

# Define custom colors for one annotation
ann_colors <- list(
  cluster = cell_type_colors,
  condition = c("VEGF" = "gold", "PDGFBB" = "firebrick")
)


p <- VlnPlot(seurat_obj,
             features = genes,
             group.by = "clust",
             split.by = "condition",   # Splits plot by sample,
             pt.size = 0,                    # no points
             cols = ann_colors$condition,
             combine = TRUE,                 # one figure with facets
             ncol = 3)                       # force 3 columns

# Adjust facet strip title size
p <- p + theme(strip.text.x = element_text(size = 6))

p <- p + 
  scale_fill_manual(
    values = ann_colors$condition,  # pick colors for samples
    name = "Condition"                                  # legend title
  ) +
  theme(legend.position = "bottom")   


ggsave("violin_endothelial_cells_geneset_4_samples.png", p, width = 7, height = 5, dpi = 600)   # PNG
ggsave("violin_endothelial_cells_geneset_4_samples.pdf", p, width = 7, height = 5)              # PDF (vector)




#############################################################################
#	Volcano plots with highlighted selected genes (endothelial cells)
#############################################################################

de_markers <- readRDS("~/applications/scstudio_Jacek/dea/tokens/heart_organoid_S1_3_endothelial_cells/dea.rds")


# List of genes to label
highlight_genes <- c("DLL4", "CXCR4", "EFNB2", "JAG1", "GATA4", "NR2F2", "FLRT2", "VEGFC")


##### Differenatial expression analysis results visualisation

  df <- de_markers[["dea_S1_vs_S3_MAST"]]
  colnames(df) <- make.names(colnames(df))  # make colnames safe
  
  # Make sure Gene column exists
  if (!"Gene" %in% colnames(df)) {
    df$Gene <- rownames(df)
  }
  
  # Classification
  df$Regulation <- "Not significant"
  df$Regulation[df$Log2FC..mean. > 1 & df$Adjusted.p.value < 0.05] <- "Up-regulated"
  df$Regulation[df$Log2FC..mean. < -1 & df$Adjusted.p.value < 0.05] <- "Down-regulated"

  # Highlight genes
  df$Label <- ifelse(df$Gene %in% highlight_genes, df$Gene, NA)
  
  
  # Identify highlighted genes that are also up/downregulated
  df$HighlightColor <- NA
  df$HighlightColor[df$Label %in% df$Gene & df$Regulation == "Up-regulated"] <- "gold"
  df$HighlightColor[df$Label %in% df$Gene & df$Regulation == "Down-regulated"] <- "firebrick"

  # Plot
  p <- ggplot(df, aes(x = Log2FC..mean., 
                      y = -log10(Adjusted.p.value), 
                      color = Regulation)) +
  # First draw all points
  geom_point(alpha = 0.6, size = 1) +
  
  # Highlight points for labelled genes (larger, black border)
  geom_point(
    data = subset(df, !is.na(Label)),
    shape = 21, fill = "gold", color = "black", size = 2.5, stroke = 0.3
  ) +

  # Highlighted points: upregulated
  geom_point(
    data = subset(df, HighlightColor == "gold"),
    aes(fill = HighlightColor),
    shape = 21, color = "black", size = 2.5, stroke = 0.3
  ) +
  
  # Highlighted points: downregulated
  geom_point(
    data = subset(df, HighlightColor == "firebrick"),
    aes(fill = HighlightColor),
    shape = 21, color = "black", size = 2.5, stroke = 0.3
  ) +

    scale_color_manual(values = c("Up-regulated" = "red",
                                  "Down-regulated" = "blue",
                                  "Not significant" = "grey")) +

    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +

  # Labels for highlighted genes
  geom_text_repel(
    data = subset(df, !is.na(Label)),
    aes(label = Label),
    size = 3.5,
    color = "black",
    max.overlaps = 50,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "black",
    bg.color = "white",   # works with ggrepel >= 0.9.5
    bg.r = 0.15
  ) +

  
  # Fill scale for highlighted points
  scale_fill_identity() +

    coord_cartesian(xlim = c(-5, 5)) +
    theme_classic() +
    labs(
      x = "Average log2 fold-change",
      y = "-log10(adjusted p-value)",
      title = "DE genes: endothelial cells",
      color = "Regulation"
    )
  
  # Save
  ggsave(
    filename = "dea_S1_vs_S3_MAST_endothelial_cells_volcano_plot.png",
    plot = p, width = 6, height = 4, dpi = 300
  )




#############################################################################
#
#	Filtered data (S1 + 3, mural cells)
#
#############################################################################


##### Set working dir to MS folder

setwd("~/data/heart_organoids/scRNA-seq/docs/MS/Figures")


##### Read in data

countMatrices <- readRDS("/mnt/scratch/home/jacek/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_mural_cells/countMatrices.rds")

metadata <- readRDS("/mnt/scratch/home/jacek/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_mural_cells/metadata.rds")

dimred <- readRDS("/mnt/scratch/home/jacek/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_mural_cells/dimred.rds")

hvgs <- readRDS("/mnt/scratch/home/jacek/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_mural_cells/hvgs.rds")


names(countMatrices)

cell_type <- metadata$cell_type
length(cell_type)

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

clust <- metadata$clust_res_0.1
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
  counts = counts,
  assay = "RNA",
  names.field = 1L,
  names.delim = "_",
  meta.data = NULL,
  project = "S1_3",
)


# Make sure dimnames match the Seurat object
all(rownames(counts) == rownames(seurat_obj))
all(colnames(counts) == colnames(seurat_obj))



# Set the "data" slot for the "RNA" assay
seurat_obj[["RNA"]] <- SetAssayData(seurat_obj[["RNA"]], slot = "data", new.data = data.normalized)



# Download barcode whitelist

whitelist_1 = read.table(gzfile("/mnt/scratch/home/jacek/data/heart_organoids/scRNA-seq/Sample1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"),sep="\t")
whitelist_2 = read.table(gzfile("/mnt/scratch/home/jacek/data/heart_organoids/scRNA-seq/Sample2/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"),sep="\t")
whitelist_3 = read.table(gzfile("/mnt/scratch/home/jacek/data/heart_organoids/scRNA-seq/Sample3/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"),sep="\t")

# Combine all barcodes lists
whitelist <- rbind(whitelist_1, whitelist_2, whitelist_3)

# Number of cells
n_cells <- ncol(seurat_obj)

# Assign to cells in your matrix
valid_barcodes <- unique(unlist(whitelist))[1:n_cells]
colnames(seurat_obj) <- valid_barcodes


# 2. Add variable features
VariableFeatures(seurat_obj) <- top_hvgs


# 3. Add dimensionality reductions

# Add PCA
# Make sure it is a numeric matrix with rownames = cell names

pca <- as.matrix(pca)
rownames(pca) <- colnames(seurat_obj)  # Ensure rows are cells

seurat_obj[["pca"]] <- Seurat::CreateDimReducObject(
  embeddings = pca,
  key = "PC_",
  assay = DefaultAssay(seurat_obj)
)


# Add UMAP

umap <- as.matrix(umap)

# Set rownames of umap to match Seurat cell names (i.e., colnames of the count matrix)
rownames(umap) <- colnames(seurat_obj)

seurat_obj[["umap"]] <- Seurat::CreateDimReducObject(
  embeddings = umap,
  key = "UMAP_",
  assay = DefaultAssay(seurat_obj)
)


# Add t-SNE

tsne <- as.matrix(tsne)

# Set rownames of tsne to match Seurat cell names (i.e., colnames of the count matrix)
rownames(tsne) <- colnames(seurat_obj)

# Now you can safely create and add the tsne slot
seurat_obj[["tsne"]] <- Seurat::CreateDimReducObject(
  embeddings = tsne,
  key = "tSNE_",
  assay = DefaultAssay(seurat_obj)
)


# Add PCA
#seurat_obj[["pca"]] <- CreateDimReducObject(embeddings = as.matrix(pca), key = "PC_", assay = DefaultAssay(seurat_obj))

# Add t-SNE
#seurat_obj[["tsne"]] <- CreateDimReducObject(embeddings = as.matrix(tsne), key = "tSNE_", assay = DefaultAssay(seurat_obj))

# Add UMAP
#seurat_obj[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(umap), key = "UMAP_", assay = DefaultAssay(seurat_obj))



# 4. Add Clusters, Dataset etc to Seurat Object

# Add the cluster information as metadata
seurat_obj$dataset <- dataset
seurat_obj$clust <- clust
seurat_obj$cell_type <- cell_type
seurat_obj$percentage_mt_genes <- percentage_mt_genes
seurat_obj$doublets_score <- doublets_score
seurat_obj$total_features <- total_features
seurat_obj$library_size <- library_size
seurat_obj$library_size_normalization <- library_size_normalization



head(seurat_obj)


# Update cell type names

annotation <- seurat_obj@meta.data

annotation[ annotation$dataset == "S1", c("dataset") ] <- "VEGF"
annotation[ annotation$dataset == "S3", c("dataset") ] <- "PDGFBB"

seurat_obj <- AddMetaData(seurat_obj, annotation$dataset, col.name ="condition")
seurat_obj$condition <- factor(seurat_obj$condition, levels = c("VEGF", "PDGFBB"))


# Define desired order of cell types
desired_order <- c(0,1,2)


# Get cell (column) names ordered by cluster
cell_order <- annotation %>%
  dplyr::mutate(cell = rownames(.)) %>%
  dplyr::arrange(factor(clust, levels = desired_order)) %>%
  pull(cell)


# Convert clust to a factor with specified order
seurat_obj$clust <- factor(seurat_obj$clust, levels = desired_order)



#############################################################################
#	Volcano plots with highlighted selected genes (mural cells)
#############################################################################

de_markers <- readRDS("~/applications/scstudio_Jacek/dea/tokens/heart_organoid_S1_3_mural_cells/dea.rds")


# List of genes to label
#highlight_genes <- c("COL6A6", "MGP", "OGN", "ELN", "ITM2A", "SERPINF1", "CCDC80", "NTRK2", "CTSK", "COLEC12")
highlight_genes <- c("CXCL12", "NPR3", "NGFR", "SMTN", "RGS5", "CSPG4", "SIRPR3", "MGP", "OGN", "PDGFD", "LTBP1", "CSRP2", "MFAP4", "TLN2", "ELN", "COLEC12")

##### Differenatial expression analysis results visualisation

  df <- de_markers[["dea_S1_vs_S3_MAST"]]
  colnames(df) <- make.names(colnames(df))  # make colnames safe
  
  # Make sure Gene column exists
  if (!"Gene" %in% colnames(df)) {
    df$Gene <- rownames(df)
  }
  
  # Classification
  df$Regulation <- "Not significant"
  df$Regulation[df$Log2FC..mean. > 1 & df$Adjusted.p.value < 0.05] <- "Up-regulated"
  df$Regulation[df$Log2FC..mean. < -1 & df$Adjusted.p.value < 0.05] <- "Down-regulated"

  # Highlight genes
  df$Label <- ifelse(df$Gene %in% highlight_genes, df$Gene, NA)
  
  
  # Identify highlighted genes that are also up/downregulated
  df$HighlightColor <- NA
  df$HighlightColor[df$Label %in% df$Gene & df$Regulation == "Up-regulated"] <- "gold"
  df$HighlightColor[df$Label %in% df$Gene & df$Regulation == "Down-regulated"] <- "firebrick"
  df$HighlightColor[df$Label %in% df$Gene & df$Regulation == "Not significant"] <- "dimgray"

  # Plot
  p <- ggplot(df, aes(x = Log2FC..mean., 
                      y = -log10(Adjusted.p.value), 
                      color = Regulation)) +
  # First draw all points
  geom_point(alpha = 0.6, size = 1) +
  
  # Highlight points for labelled genes (larger, black border)
  geom_point(
    data = subset(df, !is.na(Label)),
    shape = 21, fill = "firebrick", color = "black", size = 2.5, stroke = 0.3
  ) +

  # Highlighted points: upregulated
  geom_point(
    data = subset(df, HighlightColor == "gold"),
    aes(fill = HighlightColor),
    shape = 21, color = "black", size = 2.5, stroke = 0.3
  ) +
  
  # Highlighted points: downregulated
  geom_point(
    data = subset(df, HighlightColor == "firebrick"),
    aes(fill = HighlightColor),
    shape = 21, color = "black", size = 2.5, stroke = 0.3
  ) +
  
  # Highlighted points: not significant
  geom_point(
    data = subset(df, HighlightColor == "dimgray"),
    aes(fill = HighlightColor),
    shape = 21, color = "black", size = 2.5, stroke = 0.3
  ) +
    scale_color_manual(values = c("Up-regulated" = "red",
                                  "Down-regulated" = "blue",
                                  "Not significant" = "grey")) +

    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +

  # Labels for highlighted genes
  geom_text_repel(
    data = subset(df, !is.na(Label)),
    aes(label = Label),
    size = 3.5,
    color = "black",
    max.overlaps = 50,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "black",
    bg.color = "white",   # works with ggrepel >= 0.9.5
    bg.r = 0.15
  ) +

  
  # Fill scale for highlighted points
  scale_fill_identity() +

    coord_cartesian(xlim = c(-2.5, 2.5)) +
    theme_classic() +
    labs(
      x = "Average log2 fold-change",
      y = "-log10(adjusted p-value)",
      title = "DE genes: mural cells",
      color = "Regulation"
    )
  
  # Save
  ggsave(
    filename = "dea_S1_vs_S3_MAST_mural_cells_volcano_plot.png",
    plot = p, width = 6, height = 4, dpi = 300
  )




#############################################################################
#
#	Filtered data (S1 + 3, cardiomyocytes)
#
#############################################################################


##### Set working dir to MS folder

setwd("~/data/heart_organoids/scRNA-seq/docs/MS/Figures")


##### Read in data

countMatrices <- readRDS("/mnt/scratch/home/jacek/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_cardiomyocytes/countMatrices.rds")

metadata <- readRDS("/mnt/scratch/home/jacek/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_cardiomyocytes/metadata.rds")

dimred <- readRDS("/mnt/scratch/home/jacek/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_cardiomyocytes/dimred.rds")

hvgs <- readRDS("/mnt/scratch/home/jacek/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_cardiomyocytes/hvgs.rds")


names(countMatrices)

cell_type <- metadata$cell_type
length(cell_type)

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

clust <- metadata$clust_res_0.1
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
  counts = counts,
  assay = "RNA",
  names.field = 1L,
  names.delim = "_",
  meta.data = NULL,
  project = "S1_3",
)


# Make sure dimnames match the Seurat object
all(rownames(counts) == rownames(seurat_obj))
all(colnames(counts) == colnames(seurat_obj))



# Set the "data" slot for the "RNA" assay
seurat_obj[["RNA"]] <- SetAssayData(seurat_obj[["RNA"]], slot = "data", new.data = data.normalized)



# Download barcode whitelist

whitelist_1 = read.table(gzfile("/mnt/scratch/home/jacek/data/heart_organoids/scRNA-seq/Sample1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"),sep="\t")
whitelist_2 = read.table(gzfile("/mnt/scratch/home/jacek/data/heart_organoids/scRNA-seq/Sample2/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"),sep="\t")
whitelist_3 = read.table(gzfile("/mnt/scratch/home/jacek/data/heart_organoids/scRNA-seq/Sample3/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"),sep="\t")

# Combine all barcodes lists
whitelist <- rbind(whitelist_1, whitelist_2, whitelist_3)

# Number of cells
n_cells <- ncol(seurat_obj)

# Assign to cells in your matrix
valid_barcodes <- unique(unlist(whitelist))[1:n_cells]
colnames(seurat_obj) <- valid_barcodes


# 2. Add variable features
VariableFeatures(seurat_obj) <- top_hvgs


# 3. Add dimensionality reductions

# Add PCA
# Make sure it is a numeric matrix with rownames = cell names

pca <- as.matrix(pca)
rownames(pca) <- colnames(seurat_obj)  # Ensure rows are cells

seurat_obj[["pca"]] <- Seurat::CreateDimReducObject(
  embeddings = pca,
  key = "PC_",
  assay = DefaultAssay(seurat_obj)
)


# Add UMAP

umap <- as.matrix(umap)

# Set rownames of umap to match Seurat cell names (i.e., colnames of the count matrix)
rownames(umap) <- colnames(seurat_obj)

seurat_obj[["umap"]] <- Seurat::CreateDimReducObject(
  embeddings = umap,
  key = "UMAP_",
  assay = DefaultAssay(seurat_obj)
)


# Add t-SNE

tsne <- as.matrix(tsne)

# Set rownames of tsne to match Seurat cell names (i.e., colnames of the count matrix)
rownames(tsne) <- colnames(seurat_obj)

# Now you can safely create and add the tsne slot
seurat_obj[["tsne"]] <- Seurat::CreateDimReducObject(
  embeddings = tsne,
  key = "tSNE_",
  assay = DefaultAssay(seurat_obj)
)


# Add PCA
#seurat_obj[["pca"]] <- CreateDimReducObject(embeddings = as.matrix(pca), key = "PC_", assay = DefaultAssay(seurat_obj))

# Add t-SNE
#seurat_obj[["tsne"]] <- CreateDimReducObject(embeddings = as.matrix(tsne), key = "tSNE_", assay = DefaultAssay(seurat_obj))

# Add UMAP
#seurat_obj[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(umap), key = "UMAP_", assay = DefaultAssay(seurat_obj))



# 4. Add Clusters, Dataset etc to Seurat Object

# Add the cluster information as metadata
seurat_obj$dataset <- dataset
seurat_obj$clust <- clust
seurat_obj$cell_type <- cell_type
seurat_obj$percentage_mt_genes <- percentage_mt_genes
seurat_obj$doublets_score <- doublets_score
seurat_obj$total_features <- total_features
seurat_obj$library_size <- library_size
seurat_obj$library_size_normalization <- library_size_normalization



head(seurat_obj)


# Update cell type names

annotation <- seurat_obj@meta.data

annotation[ annotation$dataset == "S1", c("dataset") ] <- "VEGF"
annotation[ annotation$dataset == "S3", c("dataset") ] <- "PDGFBB"

seurat_obj <- AddMetaData(seurat_obj, annotation$dataset, col.name ="condition")
seurat_obj$condition <- factor(seurat_obj$condition, levels = c("VEGF", "PDGFBB"))


# Define desired order of cell types
desired_order <- c(0,1,2)


# Get cell (column) names ordered by cluster
cell_order <- annotation %>%
  dplyr::mutate(cell = rownames(.)) %>%
  dplyr::arrange(factor(clust, levels = desired_order)) %>%
  pull(cell)


# Convert clust to a factor with specified order
seurat_obj$clust <- factor(seurat_obj$clust, levels = desired_order)


#############################################################################
#	t-SNE plot (cardiomyocytes)
#############################################################################

# Basic DimPlot using existing colors
p <- DimPlot(seurat_obj,
             reduction = "tsne",
             group.by = "clust",
             label = FALSE,                 # No labels on plot
             pt.size = 0.5) +               # Adjust point size
  theme_classic() +
  labs(x = "t-SNE 1", y = "t-SNE 2") +     # <-- Add axis titles
  theme(legend.position = "right",         # Legend on the side
        #axis.line = element_blank(),
        #axis.ticks = element_blank(),
        #axis.text = element_blank(),
        axis.title = element_blank())      # Clean up axis


# Preserve Current Colors
current_colors <- Seurat::DiscretePalette(length(unique(seurat_obj$clust)))
#names(current_colors) <- levels(seurat_obj$clust)

cell_type_colors <- c("midnightblue", "darkorange")
names(cell_type_colors) <- c("0", "1")


p <- p + scale_color_manual(values = cell_type_colors)

# Legend Improvements
p <- p + guides(color = guide_legend(
                  title = "Cluster",      # Legend title
                  override.aes = list(size = 3), # Legend dot size
                  ncol = 1))                # One column (vertical)

p <- p + coord_fixed()   # Ensures equal scaling of x and y axes

p <- p + ggtitle("t-SNE of cardiomyocytes")

p <- p + labs(x = "t-SNE 1", y = "t-SNE 2")


# Save the plot at 300–600 dpi (journal-ready):
ggsave("tSNE_cardiomyocytes.png", p, width = 7, height = 5, dpi = 600)  # PNG (raster)
ggsave("tSNE_cardiomyocytes.pdf", p, width = 7, height = 5)      


# Basic DimPlot using existing colors (split samples)
p <- DimPlot(seurat_obj,
             reduction = "tsne",
             group.by = "clust",
             split.by = "condition",   # Splits plot by sample,
             label = FALSE,                 # No labels on plot
             pt.size = 0.5) +               # Adjust point size
  theme_classic() +
  labs(x = "t-SNE 1", y = "t-SNE 2") +     # <-- Add axis titles
  theme(legend.position = "right",         # Legend on the side
        #axis.line = element_blank(),
        #axis.ticks = element_blank(),
        #axis.text = element_blank(),
        axis.title = element_blank())      # Clean up axis


# Preserve Current Colors
current_colors <- Seurat::DiscretePalette(length(unique(seurat_obj$clust)))
#names(current_colors) <- levels(seurat_obj$clust)


p <- p + scale_color_manual(values = cell_type_colors)

# Legend Improvements
p <- p + guides(color = guide_legend(
                  title = "Cluster",      # Legend title
                  override.aes = list(size = 3), # Legend dot size
                  ncol = 1))                # One column (vertical)

p <- p + coord_fixed()   # Ensures equal scaling of x and y axes

p <- p + ggtitle("t-SNE of cardiomyocytes")

p <- p + labs(x = "t-SNE 1", y = "t-SNE 2")


# Save the plot at 300–600 dpi (journal-ready):
ggsave("tSNE_cardiomyocytes_by_sample.png", p, width = 7, height = 5, dpi = 600)  # PNG (raster)
ggsave("tSNE_cardiomyocytes_by_sample.pdf", p, width = 7, height = 5)      




#############################################################################
#	Violin plots for selected genes (cardiomyocytes)
#############################################################################

genes <- c("MYL2", "HEY2", "FTH1", "FTL")


# Violin plot grouped by clusters
p <- VlnPlot(seurat_obj,
             features = genes,
             group.by = "clust",
             pt.size = 0,                    # no points
             cols = cell_type_colors,
             combine = TRUE,                 # one figure with facets
             ncol = 2)                       # force 2 columns (2x2 layout for 4 genes)


# Adjust facet strip title size
p <- p + theme(strip.text.x = element_text(size = 6))

ggsave("violin_cardiomyocytes_clusters.png", p, width = 4, height = 5, dpi = 600)   # PNG
ggsave("violin_cardiomyocytes_clusters.pdf", p, width = 4, height = 5)              # PDF (vector)



# Violin plot grouped by clusters (split by sample)

# Define custom colors for one annotation
ann_colors <- list(
  cluster = cell_type_colors,
  condition = c("VEGF" = "gold", "PDGFBB" = "firebrick")
)


p <- VlnPlot(seurat_obj,
             features = genes,
             group.by = "clust",
             split.by = "condition",   # Splits plot by sample,
             pt.size = 0,                    # no points
             cols = ann_colors$condition,
             combine = TRUE,                 # one figure with facets
             ncol = 2)                       # force 2 columns (2x2 layout for 4 genes)

# Adjust facet strip title size
p <- p + theme(strip.text.x = element_text(size = 6))

p <- p + 
  scale_fill_manual(
    values = ann_colors$condition,  # pick colors for samples
    name = "Condition"                                  # legend title
  ) +
  theme(legend.position = "right")   


ggsave("violin_cardiomyocytes_samples.png", p, width = 6, height = 5, dpi = 600)   # PNG
ggsave("violin_cardiomyocytes_samples.pdf", p, width = 6, height = 5)              # PDF (vector)





#############################################################################
#
#	Filtered data (S1 + 3, Epicardium/EPDCs)
#
#############################################################################


##### Set working dir to MS folder

setwd("~/data/heart_organoids/scRNA-seq/docs/MS/Figures")


##### Read in data

countMatrices <- readRDS("/mnt/scratch/home/jacek/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_epicardium_EPDCs/countMatrices.rds")

metadata <- readRDS("/mnt/scratch/home/jacek/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_epicardium_EPDCs/metadata.rds")

dimred <- readRDS("/mnt/scratch/home/jacek/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_epicardium_EPDCs/dimred.rds")

hvgs <- readRDS("/mnt/scratch/home/jacek/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_epicardium_EPDCs/hvgs.rds")


names(countMatrices)

cell_type <- metadata$cell_type
length(cell_type)

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

clust <- metadata$clust_res_0.4
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
  counts = counts,
  assay = "RNA",
  names.field = 1L,
  names.delim = "_",
  meta.data = NULL,
  project = "S1_3",
)


# Make sure dimnames match the Seurat object
all(rownames(counts) == rownames(seurat_obj))
all(colnames(counts) == colnames(seurat_obj))



# Set the "data" slot for the "RNA" assay
seurat_obj[["RNA"]] <- SetAssayData(seurat_obj[["RNA"]], slot = "data", new.data = data.normalized)



# Download barcode whitelist

whitelist_1 = read.table(gzfile("/mnt/scratch/home/jacek/data/heart_organoids/scRNA-seq/Sample1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"),sep="\t")
whitelist_2 = read.table(gzfile("/mnt/scratch/home/jacek/data/heart_organoids/scRNA-seq/Sample2/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"),sep="\t")
whitelist_3 = read.table(gzfile("/mnt/scratch/home/jacek/data/heart_organoids/scRNA-seq/Sample3/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"),sep="\t")

# Combine all barcodes lists
whitelist <- rbind(whitelist_1, whitelist_2, whitelist_3)

# Number of cells
n_cells <- ncol(seurat_obj)

# Assign to cells in your matrix
valid_barcodes <- unique(unlist(whitelist))[1:n_cells]
colnames(seurat_obj) <- valid_barcodes


# 2. Add variable features
VariableFeatures(seurat_obj) <- top_hvgs


# 3. Add dimensionality reductions

# Add PCA
# Make sure it is a numeric matrix with rownames = cell names

pca <- as.matrix(pca)
rownames(pca) <- colnames(seurat_obj)  # Ensure rows are cells

seurat_obj[["pca"]] <- Seurat::CreateDimReducObject(
  embeddings = pca,
  key = "PC_",
  assay = DefaultAssay(seurat_obj)
)


# Add UMAP

umap <- as.matrix(umap)

# Set rownames of umap to match Seurat cell names (i.e., colnames of the count matrix)
rownames(umap) <- colnames(seurat_obj)

seurat_obj[["umap"]] <- Seurat::CreateDimReducObject(
  embeddings = umap,
  key = "UMAP_",
  assay = DefaultAssay(seurat_obj)
)


# Add t-SNE

tsne <- as.matrix(tsne)

# Set rownames of tsne to match Seurat cell names (i.e., colnames of the count matrix)
rownames(tsne) <- colnames(seurat_obj)

# Now you can safely create and add the tsne slot
seurat_obj[["tsne"]] <- Seurat::CreateDimReducObject(
  embeddings = tsne,
  key = "tSNE_",
  assay = DefaultAssay(seurat_obj)
)


# Add PCA
#seurat_obj[["pca"]] <- CreateDimReducObject(embeddings = as.matrix(pca), key = "PC_", assay = DefaultAssay(seurat_obj))

# Add t-SNE
#seurat_obj[["tsne"]] <- CreateDimReducObject(embeddings = as.matrix(tsne), key = "tSNE_", assay = DefaultAssay(seurat_obj))

# Add UMAP
#seurat_obj[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(umap), key = "UMAP_", assay = DefaultAssay(seurat_obj))



# 4. Add Clusters, Dataset etc to Seurat Object

# Add the cluster information as metadata
seurat_obj$dataset <- dataset
seurat_obj$clust <- clust
seurat_obj$cell_type <- cell_type
seurat_obj$percentage_mt_genes <- percentage_mt_genes
seurat_obj$doublets_score <- doublets_score
seurat_obj$total_features <- total_features
seurat_obj$library_size <- library_size
seurat_obj$library_size_normalization <- library_size_normalization



head(seurat_obj)


# Merge clusters 0 and 1 into one cluster and re-name all clusters accordingly

# Extract the current cluster identities
# Subset object to include specific clusters
levels(x = seurat_obj)
Idents(object = seurat_obj) <- 'clust'
levels(x = seurat_obj)

current_clusters <- Idents(seurat_obj)

# Convert to a character vector so we can edit labels
new_clusters <- as.character(current_clusters)

# Merge clusters 0 ans 1
new_clusters[new_clusters %in% c("0", "1")] <- "0_1"

# Rename all clusters
# Create a mapping for the new names
name_mapping <- c(
  "0_1" = "0",
  "2"   = "1",
  "3"   = "2",
  "4"   = "3"
)

# Apply mapping
new_clusters <- name_mapping[new_clusters]
names(new_clusters) <- names(current_clusters) 

seurat_obj$clust <- new_clusters

# Set the updated identities back to Seurat
Idents(object = seurat_obj) <- 'clust'

# Check result
table(Idents(seurat_obj))


# Update cell type names

annotation <- seurat_obj@meta.data

annotation[ annotation$dataset == "S1", c("dataset") ] <- "VEGF"
annotation[ annotation$dataset == "S3", c("dataset") ] <- "PDGFBB"

seurat_obj <- AddMetaData(seurat_obj, annotation$dataset, col.name ="condition")
seurat_obj$condition <- factor(seurat_obj$condition, levels = c("VEGF", "PDGFBB"))


# Define desired order of cell types
desired_order <- c(0,1,2,3)


# Get cell (column) names ordered by cluster
cell_order <- annotation %>%
  dplyr::mutate(cell = rownames(.)) %>%
  dplyr::arrange(factor(clust, levels = desired_order)) %>%
  pull(cell)


# Convert clust to a factor with specified order
seurat_obj$clust <- factor(seurat_obj$clust, levels = desired_order)


#############################################################################
#	t-SNE plot (Epicardium/EPDCs)
#############################################################################

# Basic DimPlot using existing colors
p <- DimPlot(seurat_obj,
             reduction = "tsne",
             group.by = "clust",
             label = FALSE,                 # No labels on plot
             pt.size = 0.5) +               # Adjust point size
  theme_classic() +
  labs(x = "t-SNE 1", y = "t-SNE 2") +     # <-- Add axis titles
  theme(legend.position = "right",         # Legend on the side
        #axis.line = element_blank(),
        #axis.ticks = element_blank(),
        #axis.text = element_blank(),
        axis.title = element_blank())      # Clean up axis


# Preserve Current Colors
current_colors <- Seurat::DiscretePalette(length(unique(seurat_obj$clust)))
#names(current_colors) <- levels(seurat_obj$clust)

cell_type_colors <- c("darkolivegreen3", "darkgoldenrod2", "wheat4", "hotpink3")
names(cell_type_colors) <- c("0", "1", "2", "3")


p <- p + scale_color_manual(values = cell_type_colors)

# Legend Improvements
p <- p + guides(color = guide_legend(
                  title = "Cluster",      # Legend title
                  override.aes = list(size = 3), # Legend dot size
                  ncol = 1))                # One column (vertical)

p <- p + coord_fixed()   # Ensures equal scaling of x and y axes

p <- p + ggtitle("t-SNE of Epicardium/EPDCs")

p <- p + labs(x = "t-SNE 1", y = "t-SNE 2")


# Save the plot at 300–600 dpi (journal-ready):
ggsave("tSNE_epicardium_EPDCs.png", p, width = 7, height = 5, dpi = 600)  # PNG (raster)
ggsave("tSNE_epicardium_EPDCs.pdf", p, width = 7, height = 5)      


# Basic DimPlot using existing colors (split samples)
p <- DimPlot(seurat_obj,
             reduction = "tsne",
             group.by = "clust",
             split.by = "condition",   # Splits plot by sample,
             label = FALSE,                 # No labels on plot
             pt.size = 0.5) +               # Adjust point size
  theme_classic() +
  labs(x = "t-SNE 1", y = "t-SNE 2") +     # <-- Add axis titles
  theme(legend.position = "right",         # Legend on the side
        #axis.line = element_blank(),
        #axis.ticks = element_blank(),
        #axis.text = element_blank(),
        axis.title = element_blank())      # Clean up axis


# Preserve Current Colors
current_colors <- Seurat::DiscretePalette(length(unique(seurat_obj$clust)))
#names(current_colors) <- levels(seurat_obj$clust)


p <- p + scale_color_manual(values = cell_type_colors)

# Legend Improvements
p <- p + guides(color = guide_legend(
                  title = "Cluster",      # Legend title
                  override.aes = list(size = 3), # Legend dot size
                  ncol = 1))                # One column (vertical)

p <- p + coord_fixed()   # Ensures equal scaling of x and y axes

p <- p + ggtitle("t-SNE of Epicardium/EPDCs")

p <- p + labs(x = "t-SNE 1", y = "t-SNE 2")


# Save the plot at 300–600 dpi (journal-ready):
ggsave("tSNE_epicardium_EPDCs_by_sample.png", p, width = 7, height = 6, dpi = 600)  # PNG (raster)
ggsave("tSNE_epicardium_EPDCs_by_sample.pdf", p, width = 7, height = 6)      




#############################################################################
#	Violin plots for selected genes (Epicardium/EPDCs)
#############################################################################

#########################
genes <- c("PPARG", "GPAM", "CALB1", "FASN", "ADIPOR2")
#########################

# Violin plot grouped by clusters
p <- VlnPlot(seurat_obj,
             features = genes,
             group.by = "clust",
             pt.size = 0,                    # no points
             cols = cell_type_colors,
             combine = TRUE,                 # one figure with facets
             ncol = 3)                       # force 3 columns


# Adjust facet strip title size
p <- p + theme(strip.text.x = element_text(size = 6))

ggsave("violin_epicardium_EPDCs_clusters_geneset_1.png", p, width = 6, height = 5, dpi = 600)   # PNG
ggsave("violin_epicardium_EPDCs_clusters_geneset_1.pdf", p, width = 6, height = 5)              # PDF (vector)



# Violin plot grouped by clusters (split by sample)

# Define custom colors for one annotation
ann_colors <- list(
  cluster = cell_type_colors,
  condition = c("VEGF" = "gold", "PDGFBB" = "firebrick")
)


p <- VlnPlot(seurat_obj,
             features = genes,
             group.by = "clust",
             split.by = "condition",   # Splits plot by sample,
             pt.size = 0,                    # no points
             cols = ann_colors$condition,
             combine = TRUE,                 # one figure with facets
             ncol = 3)                       # force 3 columns

# Adjust facet strip title size
p <- p + theme(strip.text.x = element_text(size = 6))

p <- p + 
  scale_fill_manual(
    values = ann_colors$condition,  # pick colors for samples
    name = "Condition"                                  # legend title
  ) +
  theme(legend.position = "bottom")   


ggsave("violin_epicardium_EPDCs_samples_geneset_1.png", p, width = 6, height = 5, dpi = 600)   # PNG
ggsave("violin_epicardium_EPDCs_samples_geneset_1.pdf", p, width = 6, height = 5)              # PDF (vector)




#########################
genes <- c("TNNT2", "MYH7", "TTN")
#########################

# Violin plot grouped by clusters
p <- VlnPlot(seurat_obj,
             features = genes,
             group.by = "clust",
             pt.size = 0,                    # no points
             cols = cell_type_colors,
             combine = TRUE,                 # one figure with facets
             ncol = 3)                       # force 3 columns


# Adjust facet strip title size
p <- p + theme(strip.text.x = element_text(size = 6))

ggsave("violin_epicardium_EPDCs_clusters_geneset_2.png", p, width = 6, height = 3, dpi = 600)   # PNG
ggsave("violin_epicardium_EPDCs_clusters_geneset_2.pdf", p, width = 6, height = 3)              # PDF (vector)



# Violin plot grouped by clusters (split by sample)

# Define custom colors for one annotation
ann_colors <- list(
  cluster = cell_type_colors,
  condition = c("VEGF" = "gold", "PDGFBB" = "firebrick")
)


p <- VlnPlot(seurat_obj,
             features = genes,
             group.by = "clust",
             split.by = "condition",   # Splits plot by sample,
             pt.size = 0,                    # no points
             cols = ann_colors$condition,
             combine = TRUE,                 # one figure with facets
             ncol = 3)                       # force 3 columns

# Adjust facet strip title size
p <- p + theme(strip.text.x = element_text(size = 6))

p <- p + 
  scale_fill_manual(
    values = ann_colors$condition,  # pick colors for samples
    name = "Condition"                                  # legend title
  ) +
  theme(legend.position = "right")   


ggsave("violin_epicardium_EPDCs_samples_geneset_2.png", p, width = 8, height = 3, dpi = 600)   # PNG
ggsave("violin_epicardium_EPDCs_samples_geneset_2.pdf", p, width = 8, height = 3)              # PDF (vector)




#########################
genes <- c("ANXA4", "SLC4A4", "SEMA3C", "PTPRF", "PTK2B", "C3", "ALDH1A2", "CLDN1", "MST1", "CFB")
#########################

# Violin plot grouped by clusters
p <- VlnPlot(seurat_obj,
             features = genes,
             group.by = "clust",
             pt.size = 0,                    # no points
             cols = cell_type_colors,
             combine = TRUE,                 # one figure with facets
             ncol = 5)                       # force 5 columns


# Adjust facet strip title size
p <- p + theme(strip.text.x = element_text(size = 6))

ggsave("violin_epicardium_EPDCs_clusters_geneset_3.png", p, width = 10, height = 4, dpi = 600)   # PNG
ggsave("violin_epicardium_EPDCs_clusters_geneset_3.pdf", p, width = 10, height = 4)              # PDF (vector)



# Violin plot grouped by clusters (split by sample)

# Define custom colors for one annotation
ann_colors <- list(
  cluster = cell_type_colors,
  condition = c("VEGF" = "gold", "PDGFBB" = "firebrick")
)


p <- VlnPlot(seurat_obj,
             features = genes,
             group.by = "clust",
             split.by = "condition",   # Splits plot by sample,
             pt.size = 0,                    # no points
             cols = ann_colors$condition,
             combine = TRUE,                 # one figure with facets
             ncol = 5)                       # force 5 columns

# Adjust facet strip title size
p <- p + theme(strip.text.x = element_text(size = 6))

p <- p + 
  scale_fill_manual(
    values = ann_colors$condition,  # pick colors for samples
    name = "Condition"                                  # legend title
  ) +
  theme(legend.position = "bottom")   


ggsave("violin_epicardium_EPDCs_samples_geneset_3.png", p, width = 12, height = 4, dpi = 600)   # PNG
ggsave("violin_epicardium_EPDCs_samples_geneset_3.pdf", p, width = 12, height = 4)              # PDF (vector)
#############################################################################



#########################
genes <- c("CKB", "DSP", "ANXA4", "SLC4A4", "SEMA3C", "PTPRF", "PTK2B", "WWC1", "CDH1", "C3", "ALDH1A2", "KRT19", "CLDN1", "MST1", "CFB")
#########################

# Violin plot grouped by clusters
p <- VlnPlot(seurat_obj,
             features = genes,
             group.by = "clust",
             pt.size = 0,                    # no points
             cols = cell_type_colors,
             combine = TRUE,                 # one figure with facets
             ncol = 3)                       # force 3 columns


# Adjust facet strip title size
p <- p + theme(strip.text.x = element_text(size = 6))

ggsave("violin_epicardium_EPDCs_clusters_geneset_4.png", p, width = 5, height = 10, dpi = 600)   # PNG
ggsave("violin_epicardium_EPDCs_clusters_geneset_4.pdf", p, width = 5, height = 10)              # PDF (vector)



# Violin plot grouped by clusters (split by sample)

# Define custom colors for one annotation
ann_colors <- list(
  cluster = cell_type_colors,
  condition = c("VEGF" = "gold", "PDGFBB" = "firebrick")
)


p <- VlnPlot(seurat_obj,
             features = genes,
             group.by = "clust",
             split.by = "condition",   # Splits plot by sample,
             pt.size = 0,                    # no points
             cols = ann_colors$condition,
             combine = TRUE,                 # one figure with facets
             ncol = 3)                       # force 3 columns

# Adjust facet strip title size
p <- p + theme(strip.text.x = element_text(size = 6))

p <- p + 
  scale_fill_manual(
    values = ann_colors$condition,  # pick colors for samples
    name = "Condition"                                  # legend title
  ) +
  theme(legend.position = "bottom")   


ggsave("violin_epicardium_EPDCs_samples_geneset_4.png", p, width = 8, height = 10, dpi = 600)   # PNG
ggsave("violin_epicardium_EPDCs_samples_geneset_4.pdf", p, width = 8, height = 10)              # PDF (vector)





#############################################################################
#
#	Filtered data (S1 + 3, Epicardium/EPDCs + fibroblasts, trajectory inference)
#
#############################################################################


##### Set working dir to MS folder

setwd("~/data/heart_organoids/scRNA-seq/docs/MS/Figures")


##### Read in data

countMatrices <- readRDS("/mnt/scratch/home/jacek/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_epicardium_EPDCs_fibroblasts/countMatrices.rds")

metadata <- readRDS("/mnt/scratch/home/jacek/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_epicardium_EPDCs_fibroblasts/metadata.rds")

dimred <- readRDS("/mnt/scratch/home/jacek/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_epicardium_EPDCs_fibroblasts/dimred.rds")

hvgs <- readRDS("/mnt/scratch/home/jacek/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_epicardium_EPDCs_fibroblasts/hvgs.rds")


# Also load results from trajectory inference analysis
trajectory <- readRDS("/mnt/scratch/home/jacek/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_epicardium_EPDCs_fibroblasts/trajectory_inference/heart_organoid_S1_3_epicardium_EPDCs_fibroblasts_trajectory_inference.rds")

head(trajectory$Pseudotime)



names(countMatrices)

cell_type <- metadata$cell_type
length(cell_type)

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

clust <- metadata$clust_res_0.1
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
  counts = counts,
  assay = "RNA",
  names.field = 1L,
  names.delim = "_",
  meta.data = NULL,
  project = "S1_3",
)


# Make sure dimnames match the Seurat object
all(rownames(counts) == rownames(seurat_obj))
all(colnames(counts) == colnames(seurat_obj))



# Set the "data" slot for the "RNA" assay
seurat_obj[["RNA"]] <- SetAssayData(seurat_obj[["RNA"]], slot = "data", new.data = data.normalized)



# Download barcode whitelist

whitelist_1 = read.table(gzfile("/mnt/scratch/home/jacek/data/heart_organoids/scRNA-seq/Sample1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"),sep="\t")
whitelist_2 = read.table(gzfile("/mnt/scratch/home/jacek/data/heart_organoids/scRNA-seq/Sample2/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"),sep="\t")
whitelist_3 = read.table(gzfile("/mnt/scratch/home/jacek/data/heart_organoids/scRNA-seq/Sample3/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"),sep="\t")

# Combine all barcodes lists
whitelist <- rbind(whitelist_1, whitelist_2, whitelist_3)

# Number of cells
n_cells <- ncol(seurat_obj)

# Assign to cells in your matrix
valid_barcodes <- unique(unlist(whitelist))[1:n_cells]
colnames(seurat_obj) <- valid_barcodes


# 2. Add variable features
VariableFeatures(seurat_obj) <- top_hvgs


# 3. Add dimensionality reductions

# Add PCA
# Make sure it is a numeric matrix with rownames = cell names

pca <- as.matrix(pca)
rownames(pca) <- colnames(seurat_obj)  # Ensure rows are cells

seurat_obj[["pca"]] <- Seurat::CreateDimReducObject(
  embeddings = pca,
  key = "PC_",
  assay = DefaultAssay(seurat_obj)
)


# Add UMAP

umap <- as.matrix(umap)

# Set rownames of umap to match Seurat cell names (i.e., colnames of the count matrix)
rownames(umap) <- colnames(seurat_obj)

seurat_obj[["umap"]] <- Seurat::CreateDimReducObject(
  embeddings = umap,
  key = "UMAP_",
  assay = DefaultAssay(seurat_obj)
)


# Add t-SNE

tsne <- as.matrix(tsne)

# Set rownames of tsne to match Seurat cell names (i.e., colnames of the count matrix)
rownames(tsne) <- colnames(seurat_obj)

# Now you can safely create and add the tsne slot
seurat_obj[["tsne"]] <- Seurat::CreateDimReducObject(
  embeddings = tsne,
  key = "tSNE_",
  assay = DefaultAssay(seurat_obj)
)


# Add PCA
#seurat_obj[["pca"]] <- CreateDimReducObject(embeddings = as.matrix(pca), key = "PC_", assay = DefaultAssay(seurat_obj))

# Add t-SNE
#seurat_obj[["tsne"]] <- CreateDimReducObject(embeddings = as.matrix(tsne), key = "tSNE_", assay = DefaultAssay(seurat_obj))

# Add UMAP
#seurat_obj[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(umap), key = "UMAP_", assay = DefaultAssay(seurat_obj))



# 4. Add Clusters, Dataset etc to Seurat Object

# Add the cluster information as metadata
seurat_obj$dataset <- dataset
seurat_obj$clust <- clust
seurat_obj$cell_type <- cell_type
seurat_obj$percentage_mt_genes <- percentage_mt_genes
seurat_obj$doublets_score <- doublets_score
seurat_obj$total_features <- total_features
seurat_obj$library_size <- library_size
seurat_obj$library_size_normalization <- library_size_normalization



head(seurat_obj)


# Update cell type names

annotation <- seurat_obj@meta.data

annotation[ annotation$dataset == "S1", c("dataset") ] <- "VEGF"
annotation[ annotation$dataset == "S3", c("dataset") ] <- "PDGFBB"

seurat_obj <- AddMetaData(seurat_obj, annotation$dataset, col.name ="condition")
seurat_obj$condition <- factor(seurat_obj$condition, levels = c("VEGF", "PDGFBB"))


# Define desired order of cell types
desired_order <- c(0,1)

# Convert clust to a factor with specified order
seurat_obj$clust <- factor(seurat_obj$clust, levels = desired_order)


# Add pseudioime info
seurat_obj <- AddMetaData(seurat_obj, trajectory$Pseudotime, col.name ="pseudotime")

annotation <- seurat_obj@meta.data


# Get cell (column) names ordered by cluster
cell_order <- annotation %>%
  dplyr::mutate(cell = rownames(.)) %>%
  dplyr::arrange(pseudotime) %>%
  pull(cell)



#############################################################################
#	t-SNE plot (epicardium/EPDCs + fibroblasts)
#############################################################################

# Basic DimPlot using existing colors
p <- DimPlot(seurat_obj,
             reduction = "tsne",
             group.by = "clust",
             label = FALSE,                 # No labels on plot
             pt.size = 0.5) +               # Adjust point size
  theme_classic() +
  labs(x = "t-SNE 1", y = "t-SNE 2") +     # <-- Add axis titles
  theme(legend.position = "right",         # Legend on the side
        #axis.line = element_blank(),
        #axis.ticks = element_blank(),
        #axis.text = element_blank(),
        axis.title = element_blank())      # Clean up axis


# Preserve Current Colors
current_colors <- Seurat::DiscretePalette(length(unique(seurat_obj$clust)))
#names(current_colors) <- levels(seurat_obj$clust)

cell_type_colors <- c("mediumpurple2", "yellow3")
names(cell_type_colors) <- c("0", "1")


p <- p + scale_color_manual(values = cell_type_colors)

# Legend Improvements
p <- p + guides(color = guide_legend(
                  title = "Cluster",      # Legend title
                  override.aes = list(size = 3), # Legend dot size
                  ncol = 1))                # One column (vertical)

p <- p + coord_fixed()   # Ensures equal scaling of x and y axes

p <- p + ggtitle("t-SNE of epicardium/EPDCs + fibroblasts")

p <- p + labs(x = "t-SNE 1", y = "t-SNE 2")


# Save the plot at 300–600 dpi (journal-ready):
ggsave("tSNE_epicardium_EPDCs_fibroblasts.png", p, width = 7, height = 5, dpi = 600)  # PNG (raster)
ggsave("tSNE_epicardium_EPDCs_fibroblasts.pdf", p, width = 7, height = 5)      


# Basic DimPlot using existing colors (split samples)
p <- DimPlot(seurat_obj,
             reduction = "tsne",
             group.by = "clust",
             split.by = "condition",   # Splits plot by sample,
             label = FALSE,                 # No labels on plot
             pt.size = 0.5) +               # Adjust point size
  theme_classic() +
  labs(x = "t-SNE 1", y = "t-SNE 2") +     # <-- Add axis titles
  theme(legend.position = "right",         # Legend on the side
        #axis.line = element_blank(),
        #axis.ticks = element_blank(),
        #axis.text = element_blank(),
        axis.title = element_blank())      # Clean up axis


# Preserve Current Colors
current_colors <- Seurat::DiscretePalette(length(unique(seurat_obj$clust)))
#names(current_colors) <- levels(seurat_obj$clust)


p <- p + scale_color_manual(values = cell_type_colors)

# Legend Improvements
p <- p + guides(color = guide_legend(
                  title = "Cluster",      # Legend title
                  override.aes = list(size = 3), # Legend dot size
                  ncol = 1))                # One column (vertical)

p <- p + coord_fixed()   # Ensures equal scaling of x and y axes

p <- p + ggtitle("t-SNE of epicardium/EPDCs + fibroblasts")

p <- p + labs(x = "t-SNE 1", y = "t-SNE 2")


# Save the plot at 300–600 dpi (journal-ready):
ggsave("tSNE_epicardium_EPDCs_fibroblasts_by_sample.png", p, width = 7, height = 5, dpi = 600)  # PNG (raster)
ggsave("tSNE_epicardium_EPDCs_fibroblasts_by_sample.pdf", p, width = 7, height = 5)      



#############################################################################
#	Dotplot for selected genes (epicardium/EPDCs + fibroblasts)
############################################################################


####################
genes <- c("ANXA4", "SLC4A4", "SEMA3C", "PTPRF", "PTK2B", "WWC1", "C3", "ALDH1A2", "KRT19", "CLDN1", "MST1", "CFB", "GPC3", "BICC1", "FLRT2", "DAPK1", "SLIT3", "ARID5B", "TBX18", "WT1", "DCN", "PTN", "LRRC17", "BCHE", "IGFBP7", "TIMP1", "PRDX4", "KRT18", "KRT8", "GSTM3")
####################


# Reverse the y-axis order (top-to-bottom)
seurat_obj$clust <- factor(seurat_obj$clust, levels = rev(desired_order))

# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "lightgrey", high = "blue") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#p <- p + coord_flip()


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts.png", p, width = 10, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts.pdf", p, width = 10, height = 4)  



# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "lightgrey", high = "red") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_grey_red.png", p, width = 10, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_grey_red.pdf", p, width = 10, height = 4)  


# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "blue", high = "red") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_blue_red.png", p, width = 10, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_blue_red.pdf", p, width = 10, height = 4)  



# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "royalblue", high = "red3") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_blue_red_2.png", p, width = 10, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_blue_red_2.pdf", p, width = 10, height = 4)  



# Generate dot plot grouped by cluster
dotplot_colors <- viridis::viridis(2)
p <- DotPlot(seurat_obj,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = dotplot_colors[1], high = dotplot_colors[2]) +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_blue_yellow.png", p, width = 10, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_blue_yellow.pdf", p, width = 10, height = 4)  



# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "lightgrey", high = "darkgoldenrod2") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_darkgoldenrod.png", p, width = 10, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_darkgoldenrod.pdf", p, width = 10, height = 4)  



# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "khaki", high = "darkred") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_yellow_darkred.png", p, width = 10, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_yellow_darkred.pdf", p, width = 10, height = 4)  




####################
genes <- c("CKB", "DSP", "ANXA4", "SLC4A4", "SEMA3C", "PTPRF", "PTK2B", "WWC1", "CDH1", "C3", "ALDH1A2", "KRT19", "CLDN1", "MST1", "CFB")
####################


# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "lightgrey", high = "darkgoldenrod2") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  labs(x = "Gene", y = "Cluster") +  # Change axis titles
  coord_flip() +
  theme(
    axis.text.x = element_text(size = 8, angle = 0, hjust = 1),  # tick labels
    axis.text.y = element_text(size = 8),                       # tick labels
    axis.title.x = element_text(size = 10),                     # axis title
    axis.title.y = element_text(size = 10),                     # axis title
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.4, "cm")
  )


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_geneset_1_darkgoldenrod.png", p, width = 3, height = 3.5, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_geneset_1_darkgoldenrod.pdf", p, width = 3, height = 3.5)  


# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "khaki", high = "darkred") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  labs(x = "Gene", y = "Cluster") +  # Change axis titles
  coord_flip() +
  theme(
    axis.text.x = element_text(size = 8, angle = 0, hjust = 1),  # tick labels
    axis.text.y = element_text(size = 8),                       # tick labels
    axis.title.x = element_text(size = 10),                     # axis title
    axis.title.y = element_text(size = 10),                     # axis title
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.4, "cm")
  )


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_geneset_1_yellow_darkred.png", p, width = 3, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_geneset_1_yellow_darkred.pdf", p, width = 3, height = 4)  




####################
genes <- c("GPC3", "BICC1", "FLRT2", "DAPK1", "SLIT3", "ARID5B", "TBX18", "WT1")
####################


# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "lightgrey", high = "darkgoldenrod2") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_geneset_2_darkgoldenrod.png", p, width = 5, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_geneset_2_darkgoldenrod.pdf", p, width = 5, height = 4)  



# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "khaki", high = "darkred") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_geneset_2_yellow_darkred.png", p, width = 5, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_geneset_2_yellow_darkred.pdf", p, width = 5, height = 4)  




####################
genes <- c("DCN", "PTN", "LRRC17", "BCHE")
####################


# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "lightgrey", high = "darkgoldenrod2") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_geneset_3_darkgoldenrod.png", p, width = 4, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_geneset_3_darkgoldenrod.pdf", p, width = 4, height = 4)  



# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "khaki", high = "darkred") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_geneset_3_yellow_darkred.png", p, width = 4, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_geneset_3_yellow_darkred.pdf", p, width = 4, height = 4)  



####################
genes <- c("IGFBP7", "TIMP1", "PRDX4", "KRT18", "KRT8", "GSTM3")
####################


# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "lightgrey", high = "darkgoldenrod2") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_geneset_4_darkgoldenrod.png", p, width = 5, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_geneset_4_darkgoldenrod.pdf", p, width = 5, height = 4)  



# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "khaki", high = "darkred") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_geneset_4_yellow_darkred.png", p, width = 5, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_geneset_4_yellow_darkred.pdf", p, width = 5, height = 4)  
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_yellow_darkred.pdf", p, width = 10, height = 4)  




####################
genes <- c("GPC3", "BICC1", "FLRT2", "DAPK1", "SLIT3", "ARID5B", "TBX18", "WT1", "DCN", "PTN", "LRRC17", "BCHE", "IGFBP7", "TIMP1", "PRDX4", "KRT18", "KRT8", "GSTM3")
####################


# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "lightgrey", high = "darkgoldenrod2") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  labs(x = "Gene", y = "Cluster") +  # Change axis titles
  coord_flip() +
  theme(
    axis.text.x = element_text(size = 8, angle = 0, hjust = 1),  # tick labels
    axis.text.y = element_text(size = 8),                       # tick labels
    axis.title.x = element_text(size = 10),                     # axis title
    axis.title.y = element_text(size = 10),                     # axis title
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.4, "cm")
  )


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_geneset_5_darkgoldenrod.png", p, width = 3, height = 4.5, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_geneset_5_darkgoldenrod.pdf", p, width = 3, height = 4.5)  


# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "khaki", high = "darkred") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  labs(x = "Gene", y = "Cluster") +  # Change axis titles
  coord_flip() +
  theme(
    axis.text.x = element_text(size = 8, angle = 0, hjust = 1),  # tick labels
    axis.text.y = element_text(size = 8),                       # tick labels
    axis.title.x = element_text(size = 10),                     # axis title
    axis.title.y = element_text(size = 10),                     # axis title
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.4, "cm")
  )


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_geneset_5_yellow_darkred.png", p, width = 3, height = 4.5, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_geneset_5_yellow_darkred.pdf", p, width = 3, height = 4.5)  
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_yellow_darkred.pdf", p, width = 10, height = 4)  



# Reverse the y-axis order again
seurat_obj$clust <- factor(seurat_obj$clust, levels = desired_order)



#############################################################################
#	Ridgeline plots for selected genes (epicardium/EPDCs + fibroblasts)
############################################################################

# Ridgeline plots showing the distribution of pseudotime values for cells expressing each selected gene. For each gene (y-axis), only cells with non-zero expression are included, and their pseudotime densities (x-axis) are visualized. Ridges are color-coded by cluster identity, enabling comparison of the temporal distribution of gene expression across clusters.

# Ridgeline plots of pseudotime distributions for cells expressing selected genes (>0 expression), stratified by cluster.

####################
genes <- c("ANXA4", "SLC4A4", "SEMA3C", "PTPRF", "PTK2B", "WWC1", "C3", "ALDH1A2", "KRT19", "CLDN1", "MST1", "CFB", "GPC3", "BICC1", "FLRT2", "DAPK1", "SLIT3", "ARID5B", "TBX18", "WT1", "DCN", "PTN", "LRRC17", "BCHE", "IGFBP7", "TIMP1", "PRDX4", "KRT18", "KRT8", "GSTM3")
####################


# Extract pseudotime, cluster, and gene expression
expr_df <- FetchData(seurat_obj, vars = c("pseudotime", "clust", genes))

# Reshape into long format
expr_long <- expr_df %>%
  tidyr::pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expression") %>%
  mutate(gene = factor(gene, levels = genes))  # enforce gene order

expr_long$gene <- factor(expr_long$gene, levels = rev(genes))  # reverse order
expr_long$clust <- factor(expr_long$clust, levels = sort(names(cell_type_colors)))  # or desired order
expr_filtered <- expr_long %>% filter(expression > 0)  # keep only expressing cells


# Plot ridgeline plots grouped by cluster
p <- ggplot(expr_filtered, aes(x = pseudotime, y = gene, height = ..density.., fill = clust)) +
  geom_density_ridges(scale = 1.5, rel_min_height = 0.01, alpha = 0.8) +
  scale_fill_manual(values = cell_type_colors) +
  theme_classic() +
  labs(x = "Pseudotime", y = "Genes", fill = "Cluster")

# Save the plot at 300–600 dpi (journal-ready):
ggsave("ridgeline_epicardium_EPDCs_fibroblasts.png", p, width = 5, height = 12, dpi = 600)  # PNG (raster)
ggsave("ridgeline_epicardium_EPDCs_fibroblasts.pdf", p, width = 5, height = 12)  



####################
genes <- c("CKB", "DSP", "ANXA4", "SLC4A4", "SEMA3C", "PTPRF", "PTK2B", "WWC1", "CDH1", "C3", "ALDH1A2", "KRT19", "CLDN1", "MST1", "CFB")
####################


# Extract pseudotime, cluster, and gene expression
expr_df <- FetchData(seurat_obj, vars = c("pseudotime", "clust", genes))

# Reshape into long format
expr_long <- expr_df %>%
  tidyr::pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expression") %>%
  mutate(gene = factor(gene, levels = genes))  # enforce gene order

expr_long$gene <- factor(expr_long$gene, levels = rev(genes))  # reverse order
expr_long$clust <- factor(expr_long$clust, levels = sort(names(cell_type_colors)))  # or desired order
expr_filtered <- expr_long %>% filter(expression > 0)  # keep only expressing cells


# Plot ridgeline plots grouped by cluster
p <- ggplot(expr_filtered, aes(x = pseudotime, y = gene, height = ..density.., fill = clust)) +
  geom_density_ridges(scale = 1.5, rel_min_height = 0.01, alpha = 0.8) +
  scale_fill_manual(values = cell_type_colors) +
  theme_classic() +
  labs(x = "Pseudotime", y = "Genes", fill = "Cluster")

# Save the plot at 300–600 dpi (journal-ready):
ggsave("ridgeline_epicardium_EPDCs_fibroblasts_geneset_1.png", p, width = 5, height = 6, dpi = 600)  # PNG (raster)
ggsave("ridgeline_epicardium_EPDCs_fibroblasts_geneset_1.pdf", p, width = 5, height = 6)  



####################
genes <- c("GPC3", "BICC1", "FLRT2", "DAPK1", "SLIT3", "ARID5B", "TBX18", "WT1")
####################


# Extract pseudotime, cluster, and gene expression
expr_df <- FetchData(seurat_obj, vars = c("pseudotime", "clust", genes))

# Reshape into long format
expr_long <- expr_df %>%
  tidyr::pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expression") %>%
  mutate(gene = factor(gene, levels = genes))  # enforce gene order

expr_long$gene <- factor(expr_long$gene, levels = rev(genes))  # reverse order
expr_long$clust <- factor(expr_long$clust, levels = sort(names(cell_type_colors)))  # or desired order
expr_filtered <- expr_long %>% filter(expression > 0)  # keep only expressing cells


# Plot ridgeline plots grouped by cluster
p <- ggplot(expr_filtered, aes(x = pseudotime, y = gene, height = ..density.., fill = clust)) +
  geom_density_ridges(scale = 1.5, rel_min_height = 0.01, alpha = 0.8) +
  scale_fill_manual(values = cell_type_colors) +
  theme_classic() +
  labs(x = "Pseudotime", y = "Genes", fill = "Cluster")

# Save the plot at 300–600 dpi (journal-ready):
ggsave("ridgeline_epicardium_EPDCs_fibroblasts_geneset_2.png", p, width = 5, height = 4, dpi = 600)  # PNG (raster)
ggsave("ridgeline_epicardium_EPDCs_fibroblasts_geneset_2.pdf", p, width = 5, height = 4)  




####################
genes <- c("DCN", "PTN", "LRRC17", "BCHE")
####################


# Extract pseudotime, cluster, and gene expression
expr_df <- FetchData(seurat_obj, vars = c("pseudotime", "clust", genes))

# Reshape into long format
expr_long <- expr_df %>%
  tidyr::pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expression") %>%
  mutate(gene = factor(gene, levels = genes))  # enforce gene order

expr_long$gene <- factor(expr_long$gene, levels = rev(genes))  # reverse order
expr_long$clust <- factor(expr_long$clust, levels = sort(names(cell_type_colors)))  # or desired order
expr_filtered <- expr_long %>% filter(expression > 0)  # keep only expressing cells


# Plot ridgeline plots grouped by cluster
p <- ggplot(expr_filtered, aes(x = pseudotime, y = gene, height = ..density.., fill = clust)) +
  geom_density_ridges(scale = 1.5, rel_min_height = 0.01, alpha = 0.8) +
  scale_fill_manual(values = cell_type_colors) +
  theme_classic() +
  labs(x = "Pseudotime", y = "Genes", fill = "Cluster")

# Save the plot at 300–600 dpi (journal-ready):
ggsave("ridgeline_epicardium_EPDCs_fibroblasts_geneset_3.png", p, width = 5, height = 3, dpi = 600)  # PNG (raster)
ggsave("ridgeline_epicardium_EPDCs_fibroblasts_geneset_3.pdf", p, width = 5, height = 3)  



####################
genes <- c("IGFBP7", "TIMP1", "PRDX4", "KRT18", "KRT8", "GSTM3")
####################


# Extract pseudotime, cluster, and gene expression
expr_df <- FetchData(seurat_obj, vars = c("pseudotime", "clust", genes))

# Reshape into long format
expr_long <- expr_df %>%
  tidyr::pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expression") %>%
  mutate(gene = factor(gene, levels = genes))  # enforce gene order

expr_long$gene <- factor(expr_long$gene, levels = rev(genes))  # reverse order
expr_long$clust <- factor(expr_long$clust, levels = sort(names(cell_type_colors)))  # or desired order
expr_filtered <- expr_long %>% filter(expression > 0)  # keep only expressing cells


# Plot ridgeline plots grouped by cluster
p <- ggplot(expr_filtered, aes(x = pseudotime, y = gene, height = ..density.., fill = clust)) +
  geom_density_ridges(scale = 1.5, rel_min_height = 0.01, alpha = 0.8) +
  scale_fill_manual(values = cell_type_colors) +
  theme_classic() +
  labs(x = "Pseudotime", y = "Genes", fill = "Cluster")

# Save the plot at 300–600 dpi (journal-ready):
ggsave("ridgeline_epicardium_EPDCs_fibroblasts_geneset_4.png", p, width = 5, height = 3, dpi = 600)  # PNG (raster)
ggsave("ridgeline_epicardium_EPDCs_fibroblasts_geneset_4.pdf", p, width = 5, height = 3)  





#############################################################################
#	t-SNE plots for selected genes (epicardium/EPDCs + fibroblasts)
#############################################################################


####################
genes <- c("ANXA4", "SLC4A4", "SEMA3C", "PTPRF", "PTK2B", "WWC1", "C3", "ALDH1A2", "KRT19", "CLDN1", "MST1", "CFB", "GPC3", "BICC1", "FLRT2", "DAPK1", "SLIT3", "ARID5B", "TBX18", "WT1", "DCN", "PTN", "LRRC17", "BCHE", "IGFBP7", "TIMP1", "PRDX4", "KRT18", "KRT8", "GSTM3")
####################

# Generate one plot per gene (faceted automatically if combine=TRUE)
p <- FeaturePlot(seurat_obj,
            features = genes,
            min.cutoff = "q05",
            max.cutoff = "q90", # Uses max.cutoff to avoid oversaturation from extreme outliers
            reduction = "tsne",     # use t-SNE coordinates
            #cols = c("gray90", "lightgray", "darkred"),  # low → high expression
            cols = c("khaki", "khaki", "darkred"),  # low → high expression
            combine = TRUE,         # one panel per gene
            pt.size = 0.01,
            order = TRUE) &       # <-- ensures high-expression plotted last
  theme(
    plot.title = element_text(size=12),
    axis.title = element_text(size=10),
    axis.text  = element_text(size=8),
    legend.text = element_text(size=8),
    legend.title = element_text(size=9)
  )


# Save the plot at 300–600 dpi (journal-ready):
ggsave("t-SNE_epicardium_EPDCs_fibroblasts.png", p, width = 10, height = 14, dpi = 600)  # PNG (raster)
ggsave("t-SNE_epicardium_EPDCs_fibroblasts.pdf", p, width = 10, height = 14)  



####################
genes <- c("CKB", "DSP", "ANXA4", "SLC4A4", "SEMA3C", "PTPRF", "PTK2B", "WWC1", "CDH1", "C3", "ALDH1A2", "KRT19", "CLDN1", "MST1", "CFB")
####################

# Generate one plot per gene (faceted automatically if combine=TRUE)
p <- FeaturePlot(seurat_obj,
            features = genes,
            min.cutoff = "q05",
            max.cutoff = "q90", # Uses max.cutoff to avoid oversaturation from extreme outliers
            reduction = "tsne",     # use t-SNE coordinates
            #cols = c("gray90", "lightgray", "darkred"),  # low → high expression
            cols = c("khaki", "khaki", "darkred"),  # low → high expression
            combine = TRUE,         # one panel per gene
            pt.size = 0.01,
            order = TRUE) &       # <-- ensures high-expression plotted last
  theme(
    plot.title = element_text(size=12),
    axis.title = element_text(size=10),
    axis.text  = element_text(size=8),
    legend.text = element_text(size=8),
    legend.title = element_text(size=9)
  )


# Save the plot at 300–600 dpi (journal-ready):
ggsave("t-SNE_epicardium_EPDCs_fibroblasts_geneset_1.png", p, width = 12, height = 6, dpi = 600)  # PNG (raster)
ggsave("t-SNE_epicardium_EPDCs_fibroblasts_geneset_1.pdf", p, width = 12, height = 6)  


####################
genes <- c("GPC3", "BICC1", "FLRT2", "DAPK1", "SLIT3", "ARID5B", "TBX18", "WT1")
####################

# Generate one plot per gene (faceted automatically if combine=TRUE)
p <- FeaturePlot(seurat_obj,
            features = genes,
            min.cutoff = "q05",
            max.cutoff = "q90", # Uses max.cutoff to avoid oversaturation from extreme outliers
            reduction = "tsne",     # use t-SNE coordinates
            #cols = c("gray90", "lightgray", "darkred"),  # low → high expression
            cols = c("khaki", "khaki", "darkred"),  # low → high expression
            combine = TRUE,         # one panel per gene
            pt.size = 0.01,
            order = TRUE) &       # <-- ensures high-expression plotted last
  theme(
    plot.title = element_text(size=12),
    axis.title = element_text(size=10),
    axis.text  = element_text(size=8),
    legend.text = element_text(size=8),
    legend.title = element_text(size=9)
  )


# Save the plot at 300–600 dpi (journal-ready):
ggsave("t-SNE_epicardium_EPDCs_fibroblasts_geneset_2.png", p, width = 10, height = 6, dpi = 600)  # PNG (raster)
ggsave("t-SNE_epicardium_EPDCs_fibroblasts_geneset_2.pdf", p, width = 10, height = 6)  



####################
genes <- c("DCN", "PTN", "LRRC17", "BCHE")
####################

# Generate one plot per gene (faceted automatically if combine=TRUE)
p <- FeaturePlot(seurat_obj,
            features = genes,
            min.cutoff = "q05",
            max.cutoff = "q90", # Uses max.cutoff to avoid oversaturation from extreme outliers
            reduction = "tsne",     # use t-SNE coordinates
            #cols = c("gray90", "lightgray", "darkred"),  # low → high expression
            cols = c("khaki", "khaki", "darkred"),  # low → high expression
            combine = TRUE,         # one panel per gene
            pt.size = 0.01,
            order = TRUE) &       # <-- ensures high-expression plotted last
  theme(
    plot.title = element_text(size=12),
    axis.title = element_text(size=10),
    axis.text  = element_text(size=8),
    legend.text = element_text(size=8),
    legend.title = element_text(size=9)
  )


# Save the plot at 300–600 dpi (journal-ready):
ggsave("t-SNE_epicardium_EPDCs_fibroblasts_geneset_3.png", p, width = 10, height = 6, dpi = 600)  # PNG (raster)
ggsave("t-SNE_epicardium_EPDCs_fibroblasts_geneset_3.pdf", p, width = 10, height = 6)  


####################
genes <- c("IGFBP7", "TIMP1", "PRDX4", "KRT18", "KRT8", "GSTM3")
####################

# Generate one plot per gene (faceted automatically if combine=TRUE)
p <- FeaturePlot(seurat_obj,
            features = genes,
            min.cutoff = "q05",
            max.cutoff = "q90", # Uses max.cutoff to avoid oversaturation from extreme outliers
            reduction = "tsne",     # use t-SNE coordinates
            #cols = c("gray90", "lightgray", "darkred"),  # low → high expression
            cols = c("khaki", "khaki", "darkred"),  # low → high expression
            combine = TRUE,         # one panel per gene
            pt.size = 0.01,
            order = TRUE) &       # <-- ensures high-expression plotted last
  theme(
    plot.title = element_text(size=12),
    axis.title = element_text(size=10),
    axis.text  = element_text(size=8),
    legend.text = element_text(size=8),
    legend.title = element_text(size=9)
  )


# Save the plot at 300–600 dpi (journal-ready):
ggsave("t-SNE_epicardium_EPDCs_fibroblasts_geneset_4.png", p, width = 10, height = 8, dpi = 600)  # PNG (raster)
ggsave("t-SNE_epicardium_EPDCs_fibroblasts_geneset_4.pdf", p, width = 10, height = 8)  





#############################################################################
#	Heatmap for selected genes (epicardium/EPDCs + fibroblasts)
#############################################################################

####################
genes <- c("ANXA4", "SLC4A4", "SEMA3C", "PTPRF", "PTK2B", "WWC1", "C3", "ALDH1A2", "KRT19", "CLDN1", "MST1", "CFB", "GPC3", "BICC1", "FLRT2", "DAPK1", "SLIT3", "ARID5B", "TBX18", "WT1", "DCN", "PTN", "LRRC17", "BCHE", "IGFBP7", "TIMP1", "PRDX4", "KRT18", "KRT8", "GSTM3")
####################

# Extract scaled expression data
heatmap_data <- FetchData(seurat_obj, vars = genes)

# Reorder heatmap matrix columns and annotation rows
heatmap_data <- heatmap_data[ cell_order, ]

# Ensure scaling only happens on genes with non-zero variance
heatmap_data <- heatmap_data[apply(heatmap_data, 1, sd) > 0, ]

# Scale and remove rows with NA/Inf after scaling
heatmap_matrix <- t(scale(t(heatmap_data)))

# Remove rows with NA/Inf
heatmap_matrix <- t(heatmap_matrix[complete.cases(heatmap_matrix), ])



# Annotate cells
annotation <- seurat_obj@meta.data[, c("pseudotime", "clust", "condition", "total_features", "library_size", "percentage_mt_genes")]
colnames(annotation) <- c("pseudotime", "cluster", "condition", "total_features", "library_size", "percentage_mt_genes")

# Make sure that data and annotation are of the same length
heatmap_matrix <- heatmap_matrix[ , colnames(heatmap_matrix) %in% rownames(annotation) , drop = FALSE ]
annotation <- annotation[ rownames(annotation) %in% colnames(heatmap_matrix), , drop = FALSE ]

# Define color palette
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)


# Define custom colors for one annotation
ann_colors <- list(
  cluster = cell_type_colors,
  condition = c("VEGF" = "gold", "PDGFBB" = "firebrick"),
  pseudotime <- viridis::viridis(length((annotation$pseudotime)))
)


p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "none")

png("heatmap_epicardium_EPDCs_fibroblasts_annot.png", width = 3500, height = 2000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())


p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "column")

png("heatmap_epicardium_EPDCs_fibroblasts_scaled_annot.png", width = 3500, height = 2000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())



# Annotate cells
annotation <- seurat_obj@meta.data[, c("pseudotime", "clust")]
colnames(annotation) <- c("pseudotime", "cluster")

p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "none")

png("heatmap_epicardium_EPDCs_fibroblasts.png", width = 3500, height = 1500, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())


p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "column")

png("heatmap_epicardium_EPDCs_fibroblasts_scaled.png", width = 3500, height = 1500, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())



# Define color palette
heatmap_colors <- viridis::viridis(100)

p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "none")

png("heatmap_epicardium_EPDCs_fibroblasts_blue_yellow.png", width = 3500, height = 1500, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())


p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "column")

png("heatmap_epicardium_EPDCs_fibroblasts_blue_yellow_scaled.png", width = 3500, height = 1500, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())


# Define color palette
heatmap_colors <- colorRampPalette(c("royalblue4", "white", "red4"))(100)

p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "none")

png("heatmap_epicardium_EPDCs_fibroblasts_blue_red_2.png", width = 3500, height = 1500, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())


p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "column")

png("heatmap_epicardium_EPDCs_fibroblasts_blue_red_2_scaled.png", width = 3500, height = 1500, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())



p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = TRUE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "column")

png("heatmap_epicardium_EPDCs_fibroblasts_clust_blue_red_2_scaled.png", width = 3500, height = 1500, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())





####################
genes <- c("CKB", "DSP", "ANXA4", "SLC4A4", "SEMA3C", "PTPRF", "PTK2B", "WWC1", "CDH1", "C3", "ALDH1A2", "KRT19", "CLDN1", "MST1", "CFB")
####################

# Extract scaled expression data
heatmap_data <- FetchData(seurat_obj, vars = genes)

# Reorder heatmap matrix columns and annotation rows
heatmap_data <- heatmap_data[ cell_order, ]

# Ensure scaling only happens on genes with non-zero variance
heatmap_data <- heatmap_data[apply(heatmap_data, 1, sd) > 0, ]

# Scale and remove rows with NA/Inf after scaling
heatmap_matrix <- t(scale(t(heatmap_data)))

# Remove rows with NA/Inf
heatmap_matrix <- t(heatmap_matrix[complete.cases(heatmap_matrix), ])



# Annotate cells
annotation <- seurat_obj@meta.data[, c("pseudotime", "clust", "condition", "total_features", "library_size", "percentage_mt_genes")]
colnames(annotation) <- c("pseudotime", "cluster", "condition", "total_features", "library_size", "percentage_mt_genes")

# Make sure that data and annotation are of the same length
heatmap_matrix <- heatmap_matrix[ , colnames(heatmap_matrix) %in% rownames(annotation) , drop = FALSE ]
annotation <- annotation[ rownames(annotation) %in% colnames(heatmap_matrix), , drop = FALSE ]

# Define color palette
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)


# Define custom colors for one annotation
ann_colors <- list(
  cluster = cell_type_colors,
  condition = c("VEGF" = "gold", "PDGFBB" = "firebrick"),
  pseudotime <- viridis::viridis(length((annotation$pseudotime)))
)



p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "none")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_1_annot.png", width = 3500, height = 2000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())


p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "column")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_1_scaled_annot.png", width = 3500, height = 2000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())



# Annotate cells
annotation <- seurat_obj@meta.data[, c("pseudotime", "clust")]
colnames(annotation) <- c("pseudotime", "cluster")

p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "none")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_1.png", width = 3500, height = 1000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())


p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "column")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_1_scaled.png", width = 3500, height = 1000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())



# Define color palette
heatmap_colors <- viridis::viridis(100)

p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "none")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_1_blue_yellow.png", width = 3500, height = 1000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())


p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "column")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_1_blue_yellow_scaled.png", width = 3500, height = 1000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())




# Define color palette
heatmap_colors <- colorRampPalette(c("royalblue4", "white", "red4"))(100)

p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "none")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_1_blue_red_2.png", width = 3500, height = 1000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())


p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "column")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_1_blue_red_2_scaled.png", width = 3500, height = 1000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())


p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = TRUE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "column")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_1_clust_blue_red_2_scaled.png", width = 3500, height = 1000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())




####################
genes <- c("GPC3", "BICC1", "FLRT2", "DAPK1", "SLIT3", "ARID5B", "TBX18", "WT1")
####################

# Extract scaled expression data
heatmap_data <- FetchData(seurat_obj, vars = genes)

# Reorder heatmap matrix columns and annotation rows
heatmap_data <- heatmap_data[ cell_order, ]

# Ensure scaling only happens on genes with non-zero variance
heatmap_data <- heatmap_data[apply(heatmap_data, 1, sd) > 0, ]

# Scale and remove rows with NA/Inf after scaling
heatmap_matrix <- t(scale(t(heatmap_data)))

# Remove rows with NA/Inf
heatmap_matrix <- t(heatmap_matrix[complete.cases(heatmap_matrix), ])



# Annotate cells
annotation <- seurat_obj@meta.data[, c("pseudotime", "clust", "condition", "total_features", "library_size", "percentage_mt_genes")]
colnames(annotation) <- c("pseudotime", "cluster", "condition", "total_features", "library_size", "percentage_mt_genes")

# Make sure that data and annotation are of the same length
heatmap_matrix <- heatmap_matrix[ , colnames(heatmap_matrix) %in% rownames(annotation) , drop = FALSE ]
annotation <- annotation[ rownames(annotation) %in% colnames(heatmap_matrix), , drop = FALSE ]

# Define color palette
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)


# Define custom colors for one annotation
ann_colors <- list(
  cluster = cell_type_colors,
  condition = c("VEGF" = "gold", "PDGFBB" = "firebrick"),
  pseudotime <- viridis::viridis(length((annotation$pseudotime)))
)



p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "none")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_2_annot.png", width = 3500, height = 2000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())


p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "column")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_2_scaled_annot.png", width = 3500, height = 2000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())



# Annotate cells
annotation <- seurat_obj@meta.data[, c("pseudotime", "clust")]
colnames(annotation) <- c("pseudotime", "cluster")

p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "none")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_2.png", width = 3500, height = 1000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())


p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "column")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_2_scaled.png", width = 3500, height = 1000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())



# Define color palette
heatmap_colors <- viridis::viridis(100)

p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "none")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_2_blue_yellow.png", width = 3500, height = 1000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())


p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "column")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_2_blue_yellow_scaled.png", width = 3500, height = 1000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())




# Define color palette
heatmap_colors <- colorRampPalette(c("royalblue4", "white", "red4"))(100)

p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "none")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_2_blue_red_2.png", width = 3500, height = 1000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())


p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "column")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_2_blue_red_2_scaled.png", width = 3500, height = 1000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())


p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = TRUE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "column")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_2_clust_blue_red_2_scaled.png", width = 3500, height = 1000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())




####################
genes <- c("DCN", "PTN", "LRRC17", "BCHE")
####################

# Extract scaled expression data
heatmap_data <- FetchData(seurat_obj, vars = genes)

# Reorder heatmap matrix columns and annotation rows
heatmap_data <- heatmap_data[ cell_order, ]

# Ensure scaling only happens on genes with non-zero variance
heatmap_data <- heatmap_data[apply(heatmap_data, 1, sd) > 0, ]

# Scale and remove rows with NA/Inf after scaling
heatmap_matrix <- t(scale(t(heatmap_data)))

# Remove rows with NA/Inf
heatmap_matrix <- t(heatmap_matrix[complete.cases(heatmap_matrix), ])



# Annotate cells
annotation <- seurat_obj@meta.data[, c("pseudotime", "clust", "condition", "total_features", "library_size", "percentage_mt_genes")]
colnames(annotation) <- c("pseudotime", "cluster", "condition", "total_features", "library_size", "percentage_mt_genes")

# Make sure that data and annotation are of the same length
heatmap_matrix <- heatmap_matrix[ , colnames(heatmap_matrix) %in% rownames(annotation) , drop = FALSE ]
annotation <- annotation[ rownames(annotation) %in% colnames(heatmap_matrix), , drop = FALSE ]

# Define color palette
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)


# Define custom colors for one annotation
ann_colors <- list(
  cluster = cell_type_colors,
  condition = c("VEGF" = "gold", "PDGFBB" = "firebrick"),
  pseudotime <- viridis::viridis(length((annotation$pseudotime)))
)



p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "none")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_3_annot.png", width = 3500, height = 2000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())


p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "column")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_3_scaled_annot.png", width = 3500, height = 2000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())



# Annotate cells
annotation <- seurat_obj@meta.data[, c("pseudotime", "clust")]
colnames(annotation) <- c("pseudotime", "cluster")

p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "none")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_3.png", width = 3500, height = 1000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())


p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "column")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_3_scaled.png", width = 3500, height = 1000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())



# Define color palette
heatmap_colors <- viridis::viridis(100)

p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "none")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_3_blue_yellow.png", width = 3500, height = 1000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())


p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "column")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_3_blue_yellow_scaled.png", width = 3500, height = 1000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())




# Define color palette
heatmap_colors <- colorRampPalette(c("royalblue4", "white", "red4"))(100)

p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "none")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_3_blue_red_2.png", width = 3500, height = 1000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())


p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "column")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_3_blue_red_2_scaled.png", width = 3500, height = 1000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())



p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = TRUE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "column")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_3_clust_blue_red_2_scaled.png", width = 3500, height = 1000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())



####################
genes <- c("IGFBP7", "TIMP1", "PRDX4", "KRT18", "KRT8", "GSTM3")
####################

# Extract scaled expression data
heatmap_data <- FetchData(seurat_obj, vars = genes)

# Reorder heatmap matrix columns and annotation rows
heatmap_data <- heatmap_data[ cell_order, ]

# Ensure scaling only happens on genes with non-zero variance
heatmap_data <- heatmap_data[apply(heatmap_data, 1, sd) > 0, ]

# Scale and remove rows with NA/Inf after scaling
heatmap_matrix <- t(scale(t(heatmap_data)))

# Remove rows with NA/Inf
heatmap_matrix <- t(heatmap_matrix[complete.cases(heatmap_matrix), ])



# Annotate cells
annotation <- seurat_obj@meta.data[, c("pseudotime", "clust", "condition", "total_features", "library_size", "percentage_mt_genes")]
colnames(annotation) <- c("pseudotime", "cluster", "condition", "total_features", "library_size", "percentage_mt_genes")

# Make sure that data and annotation are of the same length
heatmap_matrix <- heatmap_matrix[ , colnames(heatmap_matrix) %in% rownames(annotation) , drop = FALSE ]
annotation <- annotation[ rownames(annotation) %in% colnames(heatmap_matrix), , drop = FALSE ]

# Define color palette
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)


# Define custom colors for one annotation
ann_colors <- list(
  cluster = cell_type_colors,
  condition = c("VEGF" = "gold", "PDGFBB" = "firebrick"),
  pseudotime <- viridis::viridis(length((annotation$pseudotime)))
)



p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "none")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_4_annot.png", width = 3500, height = 2000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())


p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "column")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_4_scaled_annot.png", width = 3500, height = 2000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())



# Annotate cells
annotation <- seurat_obj@meta.data[, c("pseudotime", "clust")]
colnames(annotation) <- c("pseudotime", "cluster")

p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "none")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_4.png", width = 3500, height = 1000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())


p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "column")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_4_scaled.png", width = 3500, height = 1000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())



# Define color palette
heatmap_colors <- viridis::viridis(100)

p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "none")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_4_blue_yellow.png", width = 3500, height = 1000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())


p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "column")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_4_blue_yellow_scaled.png", width = 3500, height = 1000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())




# Define color palette
heatmap_colors <- colorRampPalette(c("royalblue4", "white", "red4"))(100)

p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "none")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_4_blue_red_2.png", width = 3500, height = 1000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())


p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "column")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_4_blue_red_2_scaled.png", width = 3500, height = 1000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())


p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = TRUE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "column")

png("heatmap_epicardium_EPDCs_fibroblasts_geneset_4_clust_blue_red_2_scaled.png", width = 3500, height = 1000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())





#############################################################################
#	Expression pseudotime plots (epicardium/EPDCs + fibroblasts)
#############################################################################

genes1 <- c("CKB", "DSP", "ANXA4", "SLC4A4", "SEMA3C", "PTPRF", "PTK2B", "WWC1", "CDH1", "C3", "ALDH1A2", "KRT19", "CLDN1", "MST1", "CFB")


# Extract scaled expression for these genes
expr_mat <- FetchData(seurat_obj, vars = genes1)

# Combine with pseudotime
df <- cbind(pseudotime = annotation$pseudotime, expr_mat) %>%
  tidyr::pivot_longer(cols = all_of(genes1), 
                      names_to = "genes1", values_to = "expression")

p <- ggplot(df, aes(x = pseudotime, y = expression, color = genes1)) +
  geom_smooth(se = FALSE, method = "loess", span = 1) +
  theme_classic() +
  labs(x = "Pseudotime", y = "Expression", color = "Gene") +
  theme(legend.position = "right")


png("expression_vs_pseudotime_geneset_1_EPDCs_fibroblasts.png", width = 3500, height = 1500, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())



# Faceted Plots
# If overlapping curves are too busy, facet each gene separately:
p <- ggplot(df, aes(x = pseudotime, y = expression)) +
  geom_smooth(se = FALSE, color = "darkred", method = "loess") +
  facet_wrap(~ genes1, scales = "fixed", ncol = 3) +  # FIXED y-axis
  theme_bw() +
  labs(x = "Pseudotime", y = "Expression")


png("expression_vs_pseudotime_geneset_1_per_gene_EPDCs_fibroblasts.png", width = 3500, height = 2500, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())



genes2 <- c("GPC3", "BICC1", "FLRT2", "DAPK1", "SLIT3", "ARID5B", "TBX18", "WT1")


# Extract scaled expression for these genes
expr_mat <- FetchData(seurat_obj, vars = genes2)

# Combine with pseudotime
df <- cbind(pseudotime = annotation$pseudotime, expr_mat) %>%
  tidyr::pivot_longer(cols = all_of(genes2), 
                      names_to = "genes2", values_to = "expression")

p <- ggplot(df, aes(x = pseudotime, y = expression, color = genes2)) +
  geom_smooth(se = FALSE, method = "loess", span = 1) +
  theme_classic() +
  labs(x = "Pseudotime", y = "Expression", color = "Gene") +
  theme(legend.position = "right")


png("expression_vs_pseudotime_geneset_2_EPDCs_fibroblasts.png", width = 3500, height = 1500, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())



# Faceted Plots
# If overlapping curves are too busy, facet each gene separately:
p <- ggplot(df, aes(x = pseudotime, y = expression)) +
  geom_smooth(se = FALSE, color = "darkred", method = "loess") +
  facet_wrap(~ genes2, scales = "fixed", ncol = 3) +  # FIXED y-axis
  theme_bw() +
  labs(x = "Pseudotime", y = "Expression")


png("expression_vs_pseudotime_geneset_2_per_gene_EPDCs_fibroblasts.png", width = 3500, height = 2500, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())




genes3 <- c("MAL2", "AQP1", "SCG5", "SAT1", "FLRT3", "GJA1", "BCHE", "KRT8", "GSTM3", "KRT18", "TIMP1")


# Extract scaled expression for these genes
expr_mat <- FetchData(seurat_obj, vars = genes3)

# Combine with pseudotime
df <- cbind(pseudotime = annotation$pseudotime, expr_mat) %>%
  tidyr::pivot_longer(cols = all_of(genes3), 
                      names_to = "genes3", values_to = "expression")

p <- ggplot(df, aes(x = pseudotime, y = expression, color = genes3)) +
  geom_smooth(se = FALSE, method = "loess", span = 1) +
  theme_classic() +
  labs(x = "Pseudotime", y = "Expression", color = "Gene") +
  theme(legend.position = "right")


png("expression_vs_pseudotime_geneset_3_EPDCs_fibroblasts.png", width = 3500, height = 1500, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())



# Faceted Plots
# If overlapping curves are too busy, facet each gene separately:
p <- ggplot(df, aes(x = pseudotime, y = expression)) +
  geom_smooth(se = FALSE, color = "darkred", method = "loess") +
  facet_wrap(~ genes3, scales = "fixed", ncol = 3) +  # FIXED y-axis
  theme_bw() +
  labs(x = "Pseudotime", y = "Expression")


png("expression_vs_pseudotime_geneset_3_per_gene_EPDCs_fibroblasts.png", width = 3500, height = 2500, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())





genes4 <- c("LOX", "FBLN1", "THY1", "TCF21", "ASPN", "TGFBI", "MMP2", "LUM", "IGFBP7", "PRDX4", "TIMP1", "DCN", "PTN", "LRRC17")


# Extract scaled expression for these genes
expr_mat <- FetchData(seurat_obj, vars = genes4)

# Combine with pseudotime
df <- cbind(pseudotime = annotation$pseudotime, expr_mat) %>%
  tidyr::pivot_longer(cols = all_of(genes4), 
                      names_to = "genes4", values_to = "expression")

p <- ggplot(df, aes(x = pseudotime, y = expression, color = genes4)) +
  geom_smooth(se = FALSE, method = "loess", span = 1) +
  theme_classic() +
  labs(x = "Pseudotime", y = "Expression", color = "Gene") +
  theme(legend.position = "right")


png("expression_vs_pseudotime_geneset_4_EPDCs_fibroblasts.png", width = 3500, height = 1500, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())



# Faceted Plots
# If overlapping curves are too busy, facet each gene separately:
p <- ggplot(df, aes(x = pseudotime, y = expression)) +
  geom_smooth(se = FALSE, color = "darkred", method = "loess") +
  facet_wrap(~ genes4, scales = "fixed", ncol = 3) +  # FIXED y-axis
  theme_bw() +
  labs(x = "Pseudotime", y = "Expression")


png("expression_vs_pseudotime_geneset_4_per_gene_EPDCs_fibroblasts.png", width = 3500, height = 2500, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())





#############################################################################
#
#	Filtered data (S1 + 3, Epicardium/EPDCs + fibroblasts + proliferative + mural + endothelial cells, trajectory inference)
#
#############################################################################


##### Set working dir to MS folder

setwd("~/data/heart_organoids/scRNA-seq/docs/MS/Figures")


##### Read in data

countMatrices <- readRDS("/mnt/scratch/home/jacek/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells/countMatrices.rds")

metadata <- readRDS("/mnt/scratch/home/jacek/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells/metadata.rds")

dimred <- readRDS("/mnt/scratch/home/jacek/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells/dimred.rds")

hvgs <- readRDS("/mnt/scratch/home/jacek/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells/hvgs.rds")


# Also load results from trajectory inference analysis
trajectory <- readRDS("/mnt/scratch/home/jacek/applications/scstudio_Jacek/clustering/tokens/heart_organoid_S1_3_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells/trajectory_inference/heart_organoid_S1_3_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_trajectory_inference.rds")

head(trajectory$Pseudotime)



names(countMatrices)

cell_type <- metadata$cell_type
length(cell_type)

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

clust <- metadata$clust_res_0.1
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
  counts = counts,
  assay = "RNA",
  names.field = 1L,
  names.delim = "_",
  meta.data = NULL,
  project = "S1_3",
)


# Make sure dimnames match the Seurat object
all(rownames(counts) == rownames(seurat_obj))
all(colnames(counts) == colnames(seurat_obj))



# Set the "data" slot for the "RNA" assay
seurat_obj[["RNA"]] <- SetAssayData(seurat_obj[["RNA"]], slot = "data", new.data = data.normalized)



# Download barcode whitelist

whitelist_1 = read.table(gzfile("/mnt/scratch/home/jacek/data/heart_organoids/scRNA-seq/Sample1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"),sep="\t")
whitelist_2 = read.table(gzfile("/mnt/scratch/home/jacek/data/heart_organoids/scRNA-seq/Sample2/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"),sep="\t")
whitelist_3 = read.table(gzfile("/mnt/scratch/home/jacek/data/heart_organoids/scRNA-seq/Sample3/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"),sep="\t")

# Combine all barcodes lists
whitelist <- rbind(whitelist_1, whitelist_2, whitelist_3)

# Number of cells
n_cells <- ncol(seurat_obj)

# Assign to cells in your matrix
valid_barcodes <- unique(unlist(whitelist))[1:n_cells]
colnames(seurat_obj) <- valid_barcodes


# 2. Add variable features
VariableFeatures(seurat_obj) <- top_hvgs


# 3. Add dimensionality reductions

# Add PCA
# Make sure it is a numeric matrix with rownames = cell names

pca <- as.matrix(pca)
rownames(pca) <- colnames(seurat_obj)  # Ensure rows are cells

seurat_obj[["pca"]] <- Seurat::CreateDimReducObject(
  embeddings = pca,
  key = "PC_",
  assay = DefaultAssay(seurat_obj)
)


# Add UMAP

umap <- as.matrix(umap)

# Set rownames of umap to match Seurat cell names (i.e., colnames of the count matrix)
rownames(umap) <- colnames(seurat_obj)

seurat_obj[["umap"]] <- Seurat::CreateDimReducObject(
  embeddings = umap,
  key = "UMAP_",
  assay = DefaultAssay(seurat_obj)
)


# Add t-SNE

tsne <- as.matrix(tsne)

# Set rownames of tsne to match Seurat cell names (i.e., colnames of the count matrix)
rownames(tsne) <- colnames(seurat_obj)

# Now you can safely create and add the tsne slot
seurat_obj[["tsne"]] <- Seurat::CreateDimReducObject(
  embeddings = tsne,
  key = "tSNE_",
  assay = DefaultAssay(seurat_obj)
)


# Add PCA
#seurat_obj[["pca"]] <- CreateDimReducObject(embeddings = as.matrix(pca), key = "PC_", assay = DefaultAssay(seurat_obj))

# Add t-SNE
#seurat_obj[["tsne"]] <- CreateDimReducObject(embeddings = as.matrix(tsne), key = "tSNE_", assay = DefaultAssay(seurat_obj))

# Add UMAP
#seurat_obj[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(umap), key = "UMAP_", assay = DefaultAssay(seurat_obj))



# 4. Add Clusters, Dataset etc to Seurat Object

# Add the cluster information as metadata
seurat_obj$dataset <- dataset
seurat_obj$clust <- clust
seurat_obj$cell_type <- cell_type
seurat_obj$percentage_mt_genes <- percentage_mt_genes
seurat_obj$doublets_score <- doublets_score
seurat_obj$total_features <- total_features
seurat_obj$library_size <- library_size
seurat_obj$library_size_normalization <- library_size_normalization



head(seurat_obj)


# Update cell type names

annotation <- seurat_obj@meta.data

annotation[ annotation$dataset == "S1", c("dataset") ] <- "VEGF"
annotation[ annotation$dataset == "S3", c("dataset") ] <- "PDGFBB"

seurat_obj <- AddMetaData(seurat_obj, annotation$dataset, col.name ="condition")
seurat_obj$condition <- factor(seurat_obj$condition, levels = c("VEGF", "PDGFBB"))


# Define desired order of cell types
desired_order <- c(3,0,2,1,4)

# Convert clust to a factor with specified order
seurat_obj$clust <- factor(seurat_obj$clust, levels = desired_order)


# Add pseudioime info
seurat_obj <- AddMetaData(seurat_obj, trajectory$Pseudotime, col.name ="pseudotime")

annotation <- seurat_obj@meta.data


# Get cell (column) names ordered by cluster
cell_order <- annotation %>%
  dplyr::mutate(cell = rownames(.)) %>%
  dplyr::arrange(pseudotime) %>%
  pull(cell)



#############################################################################
#	t-SNE plot (Epicardium/EPDCs + fibroblasts + proliferative + mural + endothelial cells)
#############################################################################

# Basic DimPlot using existing colors
p <- DimPlot(seurat_obj,
             reduction = "tsne",
             group.by = "clust",
             label = FALSE,                 # No labels on plot
             pt.size = 0.5) +               # Adjust point size
  theme_classic() +
  labs(x = "t-SNE 1", y = "t-SNE 2") +     # <-- Add axis titles
  theme(legend.position = "right",         # Legend on the side
        #axis.line = element_blank(),
        #axis.ticks = element_blank(),
        #axis.text = element_blank(),
        axis.title = element_blank())      # Clean up axis


# Preserve Current Colors
current_colors <- Seurat::DiscretePalette(length(unique(seurat_obj$clust)))
#names(current_colors) <- levels(seurat_obj$clust)

cell_type_colors <- c("darkolivegreen3", "violetred3", "slategray2", "snow4", "darkgoldenrod3")
names(cell_type_colors) <- c("0", "1", "2", "3", "4")


#cell_type_colors <- c("yellow2", "mediumpurple2", "darkolivegreen3", "lightpink2", "coral4", "burlywood4", "midnightblue")
#names(cell_type_colors) <- c("Epicardium/EPDCs", "Fibroblasts", "Proliferative fibroblasts", "Mural cells", "Endothelial cells", "Cardiomyocytes", "Endoderm-derived epithelial cells", "Hepatic-derived cells")




p <- p + scale_color_manual(values = cell_type_colors)

# Legend Improvements
p <- p + guides(color = guide_legend(
                  title = "Cluster",      # Legend title
                  override.aes = list(size = 3), # Legend dot size
                  ncol = 1))                # One column (vertical)

p <- p + coord_fixed()   # Ensures equal scaling of x and y axes

p <- p + ggtitle("t-SNE of Epicardium/EPDCs + fibroblasts + proliferative + mural + endothelial cells")

p <- p + labs(x = "t-SNE 1", y = "t-SNE 2")


# Save the plot at 300–600 dpi (journal-ready):
ggsave("tSNE_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells.png", p, width = 7, height = 5, dpi = 600)  # PNG (raster)
ggsave("tSNE_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells.pdf", p, width = 7, height = 5)      


# Basic DimPlot using existing colors (split samples)
p <- DimPlot(seurat_obj,
             reduction = "tsne",
             group.by = "clust",
             split.by = "condition",   # Splits plot by sample,
             label = FALSE,                 # No labels on plot
             pt.size = 0.5) +               # Adjust point size
  theme_classic() +
  labs(x = "t-SNE 1", y = "t-SNE 2") +     # <-- Add axis titles
  theme(legend.position = "right",         # Legend on the side
        #axis.line = element_blank(),
        #axis.ticks = element_blank(),
        #axis.text = element_blank(),
        axis.title = element_blank())      # Clean up axis


# Preserve Current Colors
current_colors <- Seurat::DiscretePalette(length(unique(seurat_obj$clust)))
#names(current_colors) <- levels(seurat_obj$clust)


p <- p + scale_color_manual(values = cell_type_colors)

# Legend Improvements
p <- p + guides(color = guide_legend(
                  title = "Cluster",      # Legend title
                  override.aes = list(size = 3), # Legend dot size
                  ncol = 1))                # One column (vertical)

p <- p + coord_fixed()   # Ensures equal scaling of x and y axes

p <- p + ggtitle("t-SNE of Epicardium/EPDCs + fibroblasts + proliferative + mural + endothelial cells")

p <- p + labs(x = "t-SNE 1", y = "t-SNE 2")


# Save the plot at 300–600 dpi (journal-ready):
ggsave("tSNE_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_by_sample.png", p, width = 7, height = 5, dpi = 600)  # PNG (raster)
ggsave("tSNE_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_by_sample.pdf", p, width = 7, height = 5)      



#############################################################################
#	Heatmap for selected genes (Epicardium/EPDCs + fibroblasts + proliferative + mural + endothelial cells)
#############################################################################

genes <- c("TK1", "ZWINT", "GMNN", "CENPU", "GGH", "CCNB2", "UBE2T", "MAD2L1", "HMGB3", "LSM4", "CKS2", "KPNA2", "TMSB15A", "RANBP1", "UBE2C", "NASP", "CENPW", "TOP2A", "HMGN2", "STMN1", "TUBA1B", "SEMA3C", "HMGA2", "NEGR1", "COL15A1", "CRABP2", "RGS5", "BAMBI", "CDH11", "FXYD5", "COL12A1", "THBS2", "MGP", "LAMA4", "MATN2", "LAMA2", "LAMC1", "RCN3", "OGN", "GSN", "ELN", "CFH", "FBN1", "PLAC9", "POSTN", "LUM", "CST3", "SPARC", "FBLN1", "DCN", "NTRK2", "SFRP2", "OSR1", "SERPINF1", "RCN3", "DUSP1", "LOX", "F10", "CPE", "C1R", "CRABP2", "ELN", "CFH", "BGN", "SPARC", "FBLN1")


# Extract scaled expression data
heatmap_data <- FetchData(seurat_obj, vars = genes)

# Reorder heatmap matrix columns and annotation rows
heatmap_data <- heatmap_data[ cell_order, ]

# Ensure scaling only happens on genes with non-zero variance
heatmap_data <- heatmap_data[apply(heatmap_data, 1, sd) > 0, ]

# Scale and remove rows with NA/Inf after scaling
heatmap_matrix <- t(scale(t(heatmap_data)))

# Remove rows with NA/Inf
heatmap_matrix <- t(heatmap_matrix[complete.cases(heatmap_matrix), ])



# Annotate cells
annotation <- seurat_obj@meta.data[, c("pseudotime", "clust", "condition", "total_features", "library_size", "percentage_mt_genes")]
colnames(annotation) <- c("pseudotime", "cluster", "condition", "total_features", "library_size", "percentage_mt_genes")

# Make sure that data and annotation are of the same length
heatmap_matrix <- heatmap_matrix[ , colnames(heatmap_matrix) %in% rownames(annotation) , drop = FALSE ]
annotation <- annotation[ rownames(annotation) %in% colnames(heatmap_matrix), , drop = FALSE ]

# Define color palette
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)


# Define custom colors for one annotation
ann_colors <- list(
  cluster = cell_type_colors,
  condition = c("VEGF" = "gold", "PDGFBB" = "firebrick"),
  pseudotime <- viridis::viridis(length((annotation$pseudotime)))
)



p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "none")

png("heatmap_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_annot.png", width = 3500, height = 2000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())


# Annotate cells
annotation <- seurat_obj@meta.data[, c("pseudotime", "clust", "condition")]
colnames(annotation) <- c("pseudotime", "cluster", "condition")

p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "none")

png("heatmap_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells.png", width = 3500, height = 1500, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())



# Define color palette
heatmap_colors <- viridis::viridis(100)

p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "none")

png("heatmap_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_blue_yellow.png", width = 3500, height = 1500, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())




# Define color palette
heatmap_colors <- colorRampPalette(c("royalblue4", "white", "red4"))(100)

p <- pheatmap(heatmap_matrix,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         #clustering_method = "ward.D2",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 10,
         color = heatmap_colors,
         scale = "none")

png("heatmap_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_blue_red_2.png", width = 3500, height = 1500, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())



#############################################################################
#	Expression pseudotime plots (Epicardium/EPDCs + fibroblasts + proliferative + mural + endothelial cells)
#############################################################################

genes1 <- c("TK1", "ZWINT", "GMNN", "CENPU", "GGH", "CCNB2", "UBE2T", "MAD2L1", "HMGB3", "LSM4", "CKS2", "KPNA2", "TMSB15A", "RANBP1", "UBE2C", "NASP", "CENPW", "TOP2A", "HMGN2", "STMN1", "TUBA1B")


# Extract scaled expression for these genes
expr_mat <- FetchData(seurat_obj, vars = genes1)


# Make sure that data and annotation are of the same length
annotation_genes1 <- annotation[ rownames(annotation) %in% rownames(expr_mat), , drop = FALSE ]

# Combine with pseudotime
df <- cbind(pseudotime = annotation_genes1 $pseudotime, expr_mat) %>%
  tidyr::pivot_longer(cols = all_of(genes1), 
                      names_to = "genes1", values_to = "expression")

p <- ggplot(df, aes(x = pseudotime, y = expression, color = genes1)) +
  geom_smooth(se = FALSE, method = "loess", span = 1) +
  theme_classic() +
  labs(x = "Pseudotime", y = "Expression", color = "Gene") +
  theme(legend.position = "right")


png("expression_vs_pseudotime_geneset_1_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells.png", width = 3500, height = 1500, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())



# Faceted Plots
# If overlapping curves are too busy, facet each gene separately:
p <- ggplot(df, aes(x = pseudotime, y = expression)) +
  geom_smooth(se = FALSE, color = "darkred", method = "loess") +
  facet_wrap(~ genes1, scales = "fixed", ncol = 3) +  # FIXED y-axis
  theme_bw() +
  labs(x = "Pseudotime", y = "Expression")


png("expression_vs_pseudotime_geneset_1_per_gene_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells.png", width = 3500, height = 3000, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())



# Subset object to include specific clusters
levels(x = seurat_obj)
Idents(object = seurat_obj) <- 'clust'
levels(x = seurat_obj)

seurat_obj_clust_0_1_4 <- subset(seurat_obj, idents = c("0", "1", "4"))


genes2 <- c("SEMA3C", "HMGA2", "NEGR1", "COL15A1", "CRABP2", "RGS5", "BAMBI", "CDH11", "FXYD5", "COL12A1", "THBS2", "MGP", "LAMA4", "MATN2", "LAMA2", "LAMC1", "RCN3", "OGN", "GSN", "ELN", "CFH", "FBN1", "PLAC9", "POSTN", "LUM", "CST3", "SPARC", "FBLN1", "DCN", "NTRK2", "SFRP2", "OSR1", "SERPINF1", "RCN3", "DUSP1", "LOX", "F10", "CPE", "C1R", "CRABP2", "ELN", "CFH", "BGN", "SPARC", "FBLN1")


# Extract scaled expression for these genes
expr_mat <- FetchData(seurat_obj_clust_0_1_4, vars = genes2)

# Make sure that data and annotation are of the same length
annotation_genes2 <- annotation[ rownames(annotation) %in% rownames(expr_mat), , drop = FALSE ]


# Combine with pseudotime
df <- cbind(pseudotime = annotation_genes2 $pseudotime, expr_mat) %>%
  tidyr::pivot_longer(cols = all_of(genes2), 
                      names_to = "genes2", values_to = "expression")

p <- ggplot(df, aes(x = pseudotime, y = expression, color = genes2)) +
  geom_smooth(se = FALSE, method = "loess", span = 1) +
  theme_classic() +
  labs(x = "Pseudotime", y = "Expression", color = "Gene") +
  theme(legend.position = "right")


png("expression_vs_pseudotime_geneset_2_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells.png", width = 3500, height = 1500, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())


# Faceted Plots
# If overlapping curves are too busy, facet each gene separately:
p <- ggplot(df, aes(x = pseudotime, y = expression)) +
  geom_smooth(se = FALSE, color = "darkred", method = "loess") +
  facet_wrap(~ genes2, scales = "fixed", ncol = 3) +  # FIXED y-axis
  theme_bw() +
  labs(x = "Pseudotime", y = "Expression")


png("expression_vs_pseudotime_geneset_2_per_gene_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells.png", width = 3500, height = 4500, pointsize = 10, units = "px", res = 300)
p
invisible(dev.off())




#############################################################################
#
#	clust: 3 -> 0 -> 2 -> 1 -> 4 (all)
#
############################################################################


#############################################################################
#	Dotplot for selected genes (epicardium/EPDCs + fibroblasts + proliferative + mural + endothelial cells)
############################################################################


# Reverse the y-axis order (top-to-bottom)
seurat_obj$clust <- factor(seurat_obj$clust, levels = rev(desired_order))



####################
genes <- c("TK1", "ZWINT", "GMNN", "CENPU", "GGH", "CCNB2", "UBE2T", "MAD2L1", "HMGB3", "LSM4", "CKS2", "KPNA2", "TMSB15A", "RANBP1", "UBE2C", "NASP", "CENPW", "TOP2A", "HMGN2", "STMN1", "TUBA1B", "MGP", "MATN2", "RCN3", "FBN1", "LUM", "CST3", "OSR1", "CRABP2", "CFH", "BGN", "SPARC", "FBLN1", "RGS5", "BAMBI", "CDH11", "FXYD5", "LAMC1", "LAMA4", "GSN", "DUSP1")
####################


# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "lightgrey", high = "blue") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#p <- p + coord_flip()


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells.png", p, width = 12, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells.pdf", p, width = 12, height = 4)  



# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "lightgrey", high = "red") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_grey_red.png", p, width = 12, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_grey_red.pdf", p, width = 12, height = 4)  


# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "blue", high = "red") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_blue_red.png", p, width = 12, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_blue_red.pdf", p, width = 12, height = 4)  



# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "royalblue", high = "red3") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_blue_red_2.png", p, width = 12, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_blue_red_2.pdf", p, width = 12, height = 4)  



# Generate dot plot grouped by cluster
dotplot_colors <- viridis::viridis(2)
p <- DotPlot(seurat_obj,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = dotplot_colors[1], high = dotplot_colors[2]) +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_blue_yellow.png", p, width = 12, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_blue_yellow.pdf", p, width = 12, height = 4)  



# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "lightgrey", high = "darkgoldenrod2") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_darkgoldenrod.png", p, width = 12, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_darkgoldenrod.pdf", p, width = 12, height = 4)  



# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "khaki", high = "darkred") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_yellow_darkred.png", p, width = 12, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_yellow_darkred.pdf", p, width = 12, height = 4)  



#############################################################################
#	Ridgeline plots for selected genes (epicardium/EPDCs + fibroblasts + proliferative + mural + endothelial cells)
############################################################################

# Ridgeline plots showing the distribution of pseudotime values for cells expressing each selected gene. For each gene (y-axis), only cells with non-zero expression are included, and their pseudotime densities (x-axis) are visualized. Ridges are color-coded by cluster identity, enabling comparison of the temporal distribution of gene expression across clusters.

# Ridgeline plots of pseudotime distributions for cells expressing selected genes (>0 expression), stratified by cluster.


# Extract pseudotime, cluster, and gene expression
expr_df <- FetchData(seurat_obj, vars = c("pseudotime", "clust", genes))

# Reshape into long format
expr_long <- expr_df %>%
  tidyr::pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expression") %>%
  mutate(gene = factor(gene, levels = genes))  # enforce gene order

expr_long$gene <- factor(expr_long$gene, levels = rev(genes))  # reverse order
expr_long$clust <- factor(expr_long$clust, levels = sort(names(cell_type_colors)))  # or desired order
expr_filtered <- expr_long %>% filter(expression > 0)  # keep only expressing cells


# Plot ridgeline plots grouped by cluster
p <- ggplot(expr_filtered, aes(x = pseudotime, y = gene, height = ..density.., fill = clust)) +
  geom_density_ridges(scale = 1.5, rel_min_height = 0.01, alpha = 0.8) +
  scale_fill_manual(values = cell_type_colors) +
  theme_classic() +
  labs(x = "Pseudotime", y = "Genes", fill = "Cluster")

# Save the plot at 300–600 dpi (journal-ready):
ggsave("ridgeline_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells.png", p, width = 5, height = 10, dpi = 600)  # PNG (raster)
ggsave("ridgeline_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells.pdf", p, width = 5, height = 10)  




#############################################################################
#	t-SNE plots for selected genes (epicardium/EPDCs + fibroblasts + proliferative + mural + endothelial cells)
#############################################################################


# Generate one plot per gene (faceted automatically if combine=TRUE)
p <- FeaturePlot(seurat_obj,
            features = genes,
            min.cutoff = "q05",
            max.cutoff = "q95", # Uses max.cutoff to avoid oversaturation from extreme outliers
            reduction = "tsne",     # use t-SNE coordinates
            #cols = c("gray90", "lightgray", "darkred"),  # low → high expression
            cols = c("khaki", "khaki", "darkred"),  # low → high expression
            combine = TRUE,         # one panel per gene
            pt.size = 0.01,
            order = TRUE) &       # <-- ensures high-expression plotted last
  theme(
    plot.title = element_text(size=12),
    axis.title = element_text(size=10),
    axis.text  = element_text(size=8),
    legend.text = element_text(size=8),
    legend.title = element_text(size=9)
  )


# Save the plot at 300–600 dpi (journal-ready):
ggsave("t-SNE_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells.png", p, width = 10, height = 16, dpi = 600)  # PNG (raster)
ggsave("t-SNE_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells.pdf", p, width = 10, height = 16)  








#############################################################################
#
#	clust: 3 -> 0 -> 2
#
############################################################################


#############################################################################
#	Dotplot for selected genes (epicardium/EPDCs + fibroblasts + proliferative + mural + endothelial cells)
############################################################################


# Subset object to include specific clusters
levels(x = seurat_obj)
Idents(object = seurat_obj) <- 'clust'
levels(x = seurat_obj)


# Subset data
seurat_obj_clust_3_0_2 <- subset(seurat_obj, idents = c("0", "2", "3"))


# Set the y-axis order (top-to-bottom)
seurat_obj_clust_3_0_2$clust <- factor(seurat_obj_clust_3_0_2$clust, levels = desired_order)



####################
genes <- c("TK1", "ZWINT", "GMNN", "CENPU", "GGH", "CCNB2", "UBE2T", "MAD2L1", "HMGB3", "LSM4", "CKS2", "KPNA2", "TMSB15A", "RANBP1", "UBE2C", "NASP", "CENPW", "TOP2A", "HMGN2", "STMN1", "TUBA1B")
####################


# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj_clust_3_0_2,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "lightgrey", high = "blue") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#p <- p + coord_flip()


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_2.png", p, width = 8, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_2.pdf", p, width = 8, height = 4)  



# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj_clust_3_0_2,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "lightgrey", high = "red") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_2_grey_red.png", p, width = 8, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_2_grey_red.pdf", p, width = 8, height = 4)  


# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj_clust_3_0_2,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "blue", high = "red") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_2_blue_red.png", p, width = 8, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_2_blue_red.pdf", p, width = 8, height = 4)  



# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj_clust_3_0_2,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "royalblue", high = "red3") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_2_blue_red_2.png", p, width = 8, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_2_blue_red_2.pdf", p, width = 8, height = 4)  



# Generate dot plot grouped by cluster
dotplot_colors <- viridis::viridis(2)
p <- DotPlot(seurat_obj_clust_3_0_2,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = dotplot_colors[1], high = dotplot_colors[2]) +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_2_blue_yellow.png", p, width = 8, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_2_blue_yellow.pdf", p, width = 8, height = 4)  



# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj_clust_3_0_2,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "lightgrey", high = "darkgoldenrod2") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_2_darkgoldenrod.png", p, width = 8, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_2_darkgoldenrod.pdf", p, width = 8, height = 4)  



# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj_clust_3_0_2,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "khaki", high = "darkred") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  labs(x = "Gene", y = "Cluster") +  # Change axis titles
  coord_flip() +
  theme(
    axis.text.x = element_text(size = 8, angle = 0, hjust = 1),  # tick labels
    axis.text.y = element_text(size = 8),                       # tick labels
    axis.title.x = element_text(size = 10),                     # axis title
    axis.title.y = element_text(size = 10),                     # axis title
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.4, "cm")
  )


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_2_yellow_darkred.png", p, width = 3.8, height = 5, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_2_yellow_darkred.pdf", p, width = 3.8, height = 5)  



#############################################################################
#	Ridgeline plots for selected genes (epicardium/EPDCs + fibroblasts + proliferative + mural + endothelial cells)
############################################################################

# Ridgeline plots showing the distribution of pseudotime values for cells expressing each selected gene. For each gene (y-axis), only cells with non-zero expression are included, and their pseudotime densities (x-axis) are visualized. Ridges are color-coded by cluster identity, enabling comparison of the temporal distribution of gene expression across clusters.

# Ridgeline plots of pseudotime distributions for cells expressing selected genes (>0 expression), stratified by cluster.


# Extract pseudotime, cluster, and gene expression
expr_df <- FetchData(seurat_obj_clust_3_0_2, vars = c("pseudotime", "clust", genes))

# Reshape into long format
expr_long <- expr_df %>%
  tidyr::pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expression") %>%
  mutate(gene = factor(gene, levels = genes))  # enforce gene order

expr_long$gene <- factor(expr_long$gene, levels = rev(genes))  # reverse order
expr_long$clust <- factor(expr_long$clust, levels = sort(names(cell_type_colors)))  # or desired order
expr_filtered <- expr_long %>% filter(expression > 0)  # keep only expressing cells


# Plot ridgeline plots grouped by cluster
p <- ggplot(expr_filtered, aes(x = pseudotime, y = gene, height = ..density.., fill = clust)) +
  geom_density_ridges(scale = 1.5, rel_min_height = 0.01, alpha = 0.8) +
  scale_fill_manual(values = cell_type_colors) +
  theme_classic() +
  labs(x = "Pseudotime", y = "Genes", fill = "Cluster")

# Save the plot at 300–600 dpi (journal-ready):
ggsave("ridgeline_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_2.png", p, width = 5, height = 7, dpi = 600)  # PNG (raster)
ggsave("ridgeline_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_2.pdf", p, width = 5, height = 7)  




#############################################################################
#	t-SNE plots for selected genes (epicardium/EPDCs + fibroblasts + proliferative + mural + endothelial cells)
#############################################################################


# Generate one plot per gene (faceted automatically if combine=TRUE)
p <- FeaturePlot(seurat_obj_clust_3_0_2,
            features = genes,
            min.cutoff = "q05",
            max.cutoff = "q95", # Uses max.cutoff to avoid oversaturation from extreme outliers
            reduction = "tsne",     # use t-SNE coordinates
            #cols = c("gray90", "lightgray", "darkred"),  # low → high expression
            cols = c("khaki", "khaki", "darkred"),  # low → high expression
            combine = TRUE,         # one panel per gene
            pt.size = 0.01,
            order = TRUE) &       # <-- ensures high-expression plotted last
  theme(
    plot.title = element_text(size=12),
    axis.title = element_text(size=10),
    axis.text  = element_text(size=8),
    legend.text = element_text(size=8),
    legend.title = element_text(size=9)
  )


# Save the plot at 300–600 dpi (journal-ready):
ggsave("t-SNE_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_2.png", p, width = 10, height = 10, dpi = 600)  # PNG (raster)
ggsave("t-SNE_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_2.pdf", p, width = 10, height = 10)  






#############################################################################
#
#	clust: 3 -> 0 -> 1 -> 4
#
############################################################################


#############################################################################
#	Dotplot for selected genes (epicardium/EPDCs + fibroblasts + proliferative + mural + endothelial cells)
############################################################################


# Subset object to include specific clusters
levels(x = seurat_obj)
Idents(object = seurat_obj) <- 'clust'
levels(x = seurat_obj)


# Subset data
seurat_obj_clust_3_0_1_4 <- subset(seurat_obj, idents = c("0", "1", "3", "4"))


# Set the y-axis order (top-to-bottom)
seurat_obj_clust_3_0_1_4$clust <- factor(seurat_obj_clust_3_0_1_4$clust, levels = desired_order)



####################
genes <- c("RGS5" , "BAMBI" , "CDH11" , "FXYD5" , "LAMC1" , "LAMA4" , "GSN" , "FBN1")
####################


# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj_clust_3_0_1_4,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "lightgrey", high = "blue") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#p <- p + coord_flip()


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4.png", p, width = 5, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4.pdf", p, width = 5, height = 4)  



# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj_clust_3_0_1_4,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "lightgrey", high = "red") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4_grey_red.png", p, width = 5, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4_grey_red.pdf", p, width = 5, height = 4)  


# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj_clust_3_0_1_4,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "blue", high = "red") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4_blue_red.png", p, width = 5, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4_blue_red.pdf", p, width = 5, height = 4)  



# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj_clust_3_0_1_4,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "royalblue", high = "red3") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4_blue_red_2.png", p, width = 5, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4_blue_red_2.pdf", p, width = 5, height = 4)  



# Generate dot plot grouped by cluster
dotplot_colors <- viridis::viridis(2)
p <- DotPlot(seurat_obj_clust_3_0_1_4,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = dotplot_colors[1], high = dotplot_colors[2]) +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4_blue_yellow.png", p, width = 5, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4_blue_yellow.pdf", p, width = 5, height = 4)  



# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj_clust_3_0_1_4,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "lightgrey", high = "darkgoldenrod2") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4_darkgoldenrod.png", p, width = 5, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4_darkgoldenrod.pdf", p, width = 5, height = 4)  



# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj_clust_3_0_1_4,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "khaki", high = "darkred") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  labs(x = "Gene", y = "Cluster") +  # Change axis titles
  coord_flip() +
  theme(
    axis.text.x = element_text(size = 8, angle = 0, hjust = 1),  # tick labels
    axis.text.y = element_text(size = 8),                       # tick labels
    axis.title.x = element_text(size = 10),                     # axis title
    axis.title.y = element_text(size = 10),                     # axis title
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.4, "cm")
  )


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4_yellow_darkred.png", p, width = 4, height = 2.5, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4_yellow_darkred.pdf", p, width = 4, height = 2.5)  



#############################################################################
#	Ridgeline plots for selected genes (epicardium/EPDCs + fibroblasts + proliferative + mural + endothelial cells)
############################################################################

# Ridgeline plots showing the distribution of pseudotime values for cells expressing each selected gene. For each gene (y-axis), only cells with non-zero expression are included, and their pseudotime densities (x-axis) are visualized. Ridges are color-coded by cluster identity, enabling comparison of the temporal distribution of gene expression across clusters.

# Ridgeline plots of pseudotime distributions for cells expressing selected genes (>0 expression), stratified by cluster.


# Extract pseudotime, cluster, and gene expression
expr_df <- FetchData(seurat_obj_clust_3_0_1_4, vars = c("pseudotime", "clust", genes))

# Reshape into long format
expr_long <- expr_df %>%
  tidyr::pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expression") %>%
  mutate(gene = factor(gene, levels = genes))  # enforce gene order

expr_long$gene <- factor(expr_long$gene, levels = rev(genes))  # reverse order
expr_long$clust <- factor(expr_long$clust, levels = sort(names(cell_type_colors)))  # or desired order
expr_filtered <- expr_long %>% filter(expression > 0)  # keep only expressing cells


# Plot ridgeline plots grouped by cluster
p <- ggplot(expr_filtered, aes(x = pseudotime, y = gene, height = ..density.., fill = clust)) +
  geom_density_ridges(scale = 1.5, rel_min_height = 0.01, alpha = 0.8) +
  scale_fill_manual(values = cell_type_colors) +
  theme_classic() +
  labs(x = "Pseudotime", y = "Genes", fill = "Cluster")

# Save the plot at 300–600 dpi (journal-ready):
ggsave("ridgeline_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4.png", p, width = 5, height = 3, dpi = 600)  # PNG (raster)
ggsave("ridgeline_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4.pdf", p, width = 5, height = 3)  




#############################################################################
#	t-SNE plots for selected genes (epicardium/EPDCs + fibroblasts + proliferative + mural + endothelial cells)
#############################################################################


# Generate one plot per gene (faceted automatically if combine=TRUE)
p <- FeaturePlot(seurat_obj_clust_3_0_1_4,
            features = genes,
            min.cutoff = "q05",
            max.cutoff = "q95", # Uses max.cutoff to avoid oversaturation from extreme outliers
            reduction = "tsne",     # use t-SNE coordinates
            #cols = c("gray90", "lightgray", "darkred"),  # low → high expression
            cols = c("khaki", "khaki", "darkred"),  # low → high expression
            combine = TRUE,         # one panel per gene
            pt.size = 0.01,
            order = TRUE) &       # <-- ensures high-expression plotted last
  theme(
    plot.title = element_text(size=12),
    axis.title = element_text(size=10),
    axis.text  = element_text(size=8),
    legend.text = element_text(size=8),
    legend.title = element_text(size=9)
  )


# Save the plot at 300–600 dpi (journal-ready):
ggsave("t-SNE_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4.png", p, width = 10, height = 6, dpi = 600)  # PNG (raster)
ggsave("t-SNE_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4.pdf", p, width = 10, height = 6)  



#############################################################################
#	Dotplot for selected genes (epicardium/EPDCs + fibroblasts + proliferative + mural + endothelial cells)
############################################################################


# Subset object to include specific clusters
levels(x = seurat_obj)
Idents(object = seurat_obj) <- 'clust'
levels(x = seurat_obj)


# Subset data
seurat_obj_clust_3_0_1_4 <- subset(seurat_obj, idents = c("0", "1", "3", "4"))


# Set the y-axis order (top-to-bottom)
seurat_obj_clust_3_0_1_4$clust <- factor(seurat_obj_clust_3_0_1_4$clust, levels = desired_order)



####################
genes <- c("MGP", "MATN2", "RCN3", "FBN1", "LUM", "CST3", "OSR1", "CRABP2", "CFH", "BGN")
####################



# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj_clust_3_0_1_4,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "lightgrey", high = "blue") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#p <- p + coord_flip()


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4_geneset_2.png", p, width = 8, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4_geneset_2.pdf", p, width = 8, height = 4)  



# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj_clust_3_0_1_4,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "lightgrey", high = "red") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4_grey_red_geneset_2.png", p, width = 8, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4_grey_red_geneset_2.pdf", p, width = 8, height = 4)  


# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj_clust_3_0_1_4,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "blue", high = "red") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4_blue_red_geneset_2.png", p, width = 8, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4_blue_red_geneset_2.pdf", p, width = 8, height = 4)  



# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj_clust_3_0_1_4,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "royalblue", high = "red3") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4_blue_red_2_geneset_2.png", p, width = 8, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4_blue_red_2_geneset_2.pdf", p, width = 8, height = 4)  



# Generate dot plot grouped by cluster
dotplot_colors <- viridis::viridis(2)
p <- DotPlot(seurat_obj_clust_3_0_1_4,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = dotplot_colors[1], high = dotplot_colors[2]) +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4_blue_yellow_geneset_2.png", p, width = 8, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4_blue_yellow_geneset_2.pdf", p, width = 8, height = 4)  



# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj_clust_3_0_1_4,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "lightgrey", high = "darkgoldenrod2") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4_darkgoldenrod_geneset_2.png", p, width = 8, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4_darkgoldenrod_geneset_2.pdf", p, width = 8, height = 4)  



# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj_clust_3_0_1_4,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "khaki", high = "darkred") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  labs(x = "Gene", y = "Cluster") +  # Change axis titles
  coord_flip() +
  theme(
    axis.text.x = element_text(size = 8, angle = 0, hjust = 1),  # tick labels
    axis.text.y = element_text(size = 8),                       # tick labels
    axis.title.x = element_text(size = 10),                     # axis title
    axis.title.y = element_text(size = 10),                     # axis title
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.4, "cm")
  )


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4_yellow_darkred_geneset_2.png", p, width = 4, height = 3, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4_yellow_darkred_geneset_2.pdf", p, width = 4, height = 3)  



#############################################################################
#	Ridgeline plots for selected genes (epicardium/EPDCs + fibroblasts + proliferative + mural + endothelial cells)
############################################################################

# Ridgeline plots showing the distribution of pseudotime values for cells expressing each selected gene. For each gene (y-axis), only cells with non-zero expression are included, and their pseudotime densities (x-axis) are visualized. Ridges are color-coded by cluster identity, enabling comparison of the temporal distribution of gene expression across clusters.

# Ridgeline plots of pseudotime distributions for cells expressing selected genes (>0 expression), stratified by cluster.


# Extract pseudotime, cluster, and gene expression
expr_df <- FetchData(seurat_obj_clust_3_0_1_4, vars = c("pseudotime", "clust", genes))

# Reshape into long format
expr_long <- expr_df %>%
  tidyr::pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expression") %>%
  mutate(gene = factor(gene, levels = genes))  # enforce gene order

expr_long$gene <- factor(expr_long$gene, levels = rev(genes))  # reverse order
expr_long$clust <- factor(expr_long$clust, levels = sort(names(cell_type_colors)))  # or desired order
expr_filtered <- expr_long %>% filter(expression > 0)  # keep only expressing cells


# Plot ridgeline plots grouped by cluster
p <- ggplot(expr_filtered, aes(x = pseudotime, y = gene, height = ..density.., fill = clust)) +
  geom_density_ridges(scale = 1.5, rel_min_height = 0.01, alpha = 0.8) +
  scale_fill_manual(values = cell_type_colors) +
  theme_classic() +
  labs(x = "Pseudotime", y = "Genes", fill = "Cluster")

# Save the plot at 300–600 dpi (journal-ready):
ggsave("ridgeline_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4_geneset_2.png", p, width = 5, height = 6, dpi = 600)  # PNG (raster)
ggsave("ridgeline_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4_geneset_2.pdf", p, width = 5, height = 6)  




#############################################################################
#	t-SNE plots for selected genes (epicardium/EPDCs + fibroblasts + proliferative + mural + endothelial cells)
#############################################################################


# Generate one plot per gene (faceted automatically if combine=TRUE)
p <- FeaturePlot(seurat_obj_clust_3_0_1_4,
            features = genes,
            min.cutoff = "q05",
            max.cutoff = "q95", # Uses max.cutoff to avoid oversaturation from extreme outliers
            reduction = "tsne",     # use t-SNE coordinates
            #cols = c("gray90", "lightgray", "darkred"),  # low → high expression
            cols = c("khaki", "khaki", "darkred"),  # low → high expression
            combine = TRUE,         # one panel per gene
            pt.size = 0.01,
            order = TRUE) &       # <-- ensures high-expression plotted last
  theme(
    plot.title = element_text(size=12),
    axis.title = element_text(size=10),
    axis.text  = element_text(size=8),
    legend.text = element_text(size=8),
    legend.title = element_text(size=9)
  )


# Save the plot at 300–600 dpi (journal-ready):
ggsave("t-SNE_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4_geneset_2.png", p, width = 10, height = 10, dpi = 600)  # PNG (raster)
ggsave("t-SNE_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_4_geneset_2.pdf", p, width = 10, height = 10)  




#############################################################################
#
#	clust: 3 -> 0 -> 1
#
############################################################################


#############################################################################
#	Dotplot for selected genes (epicardium/EPDCs + fibroblasts + proliferative + mural + endothelial cells)
############################################################################


# Subset object to include specific clusters
levels(x = seurat_obj)
Idents(object = seurat_obj) <- 'clust'
levels(x = seurat_obj)


# Subset data
seurat_obj_clust_3_0_1 <- subset(seurat_obj, idents = c("0", "1", "3"))


# Reverse the y-axis order (top-to-bottom)
seurat_obj_clust_3_0_1$clust <- factor(seurat_obj_clust_3_0_1$clust, levels = rev(desired_order))



####################
genes <- c("MGP", "MATN2", "RCN3", "FBN1", "LUM", "CST3", "OSR1", "CRABP2", "CFH", "BGN", "SPARC", "FBLN1")
####################


# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj_clust_3_0_1,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "lightgrey", high = "blue") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#p <- p + coord_flip()


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1.png", p, width = 6, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1.pdf", p, width = 6, height = 4)  



# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj_clust_3_0_1,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "lightgrey", high = "red") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_grey_red.png", p, width = 6, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_grey_red.pdf", p, width = 6, height = 4)  


# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj_clust_3_0_1,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "blue", high = "red") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_blue_red.png", p, width = 6, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_blue_red.pdf", p, width = 6, height = 4)  



# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj_clust_3_0_1,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "royalblue", high = "red3") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_blue_red_2.png", p, width = 6, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_blue_red_2.pdf", p, width = 6, height = 4)  



# Generate dot plot grouped by cluster
dotplot_colors <- viridis::viridis(2)
p <- DotPlot(seurat_obj_clust_3_0_1,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = dotplot_colors[1], high = dotplot_colors[2]) +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_blue_yellow.png", p, width = 6, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_blue_yellow.pdf", p, width = 6, height = 4)  



# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj_clust_3_0_1,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "lightgrey", high = "darkgoldenrod2") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_darkgoldenrod.png", p, width = 6, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_darkgoldenrod.pdf", p, width = 6, height = 4)  



# Generate dot plot grouped by cluster
p <- DotPlot(seurat_obj_clust_3_0_1,
        features = genes,
        group.by = "clust",
        scale = FALSE) +
  scale_color_gradient(low = "khaki", high = "darkred") +  # color scale
  scale_size(range = c(1, 6)) +                            # dot size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Save the plot at 300–600 dpi (journal-ready):
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_yellow_darkred.png", p, width = 6, height = 4, dpi = 600)  # PNG (raster)
ggsave("dot_plot_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1_yellow_darkred.pdf", p, width = 6, height = 4)  



#############################################################################
#	Ridgeline plots for selected genes (epicardium/EPDCs + fibroblasts + proliferative + mural + endothelial cells)
############################################################################

# Ridgeline plots showing the distribution of pseudotime values for cells expressing each selected gene. For each gene (y-axis), only cells with non-zero expression are included, and their pseudotime densities (x-axis) are visualized. Ridges are color-coded by cluster identity, enabling comparison of the temporal distribution of gene expression across clusters.

# Ridgeline plots of pseudotime distributions for cells expressing selected genes (>0 expression), stratified by cluster.


# Extract pseudotime, cluster, and gene expression
expr_df <- FetchData(seurat_obj_clust_3_0_1, vars = c("pseudotime", "clust", genes))

# Reshape into long format
expr_long <- expr_df %>%
  tidyr::pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expression") %>%
  mutate(gene = factor(gene, levels = genes))  # enforce gene order

expr_long$gene <- factor(expr_long$gene, levels = rev(genes))  # reverse order
expr_long$clust <- factor(expr_long$clust, levels = sort(names(cell_type_colors)))  # or desired order
expr_filtered <- expr_long %>% filter(expression > 0)  # keep only expressing cells


# Plot ridgeline plots grouped by cluster
p <- ggplot(expr_filtered, aes(x = pseudotime, y = gene, height = ..density.., fill = clust)) +
  geom_density_ridges(scale = 1.5, rel_min_height = 0.01, alpha = 0.8) +
  scale_fill_manual(values = cell_type_colors) +
  theme_classic() +
  labs(x = "Pseudotime", y = "Genes", fill = "Cluster")

# Save the plot at 300–600 dpi (journal-ready):
ggsave("ridgeline_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1.png", p, width = 5, height = 4, dpi = 600)  # PNG (raster)
ggsave("ridgeline_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1.pdf", p, width = 5, height = 4)  




#############################################################################
#	t-SNE plots for selected genes (epicardium/EPDCs + fibroblasts + proliferative + mural + endothelial cells)
#############################################################################


# Generate one plot per gene (faceted automatically if combine=TRUE)
p <- FeaturePlot(seurat_obj_clust_3_0_1,
            features = genes,
            min.cutoff = "q05",
            max.cutoff = "q95", # Uses max.cutoff to avoid oversaturation from extreme outliers
            reduction = "tsne",     # use t-SNE coordinates
            #cols = c("gray90", "lightgray", "darkred"),  # low → high expression
            cols = c("khaki", "khaki", "darkred"),  # low → high expression
            combine = TRUE,         # one panel per gene
            pt.size = 0.01,
            order = TRUE) &       # <-- ensures high-expression plotted last
  theme(
    plot.title = element_text(size=12),
    axis.title = element_text(size=10),
    axis.text  = element_text(size=8),
    legend.text = element_text(size=8),
    legend.title = element_text(size=9)
  )


# Save the plot at 300–600 dpi (journal-ready):
ggsave("t-SNE_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1.png", p, width = 10, height = 6, dpi = 600)  # PNG (raster)
ggsave("t-SNE_epicardium_EPDCs_fibroblasts_proliferative_mural_endothelial_cells_clust_3_0_1.pdf", p, width = 10, height = 6)  




