############################################################
# Install Required Packages
############################################################

# R packages
install.packages("Seurat")
install.packages("tidyverse")
install.packages("patchwork")
install.packages("reticulate")
install.packages("tidyr")
install.packages('BiocManager')
BiocManager::install('EnhancedVolcano')

# For h5ad support
reticulate::install_miniconda()

############################################################
# Load Libraries
############################################################

library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(reticulate)
library(Matrix)
library(patchwork)
library(EnhancedVolcano)

# convert .h5ad to .h5seurat 
Convert( "dataset/human_immune_cells.h5ad", dest = "h5seurat", overwrite = TRUE) 
# I faced many issues while loading h5seurat file, so I used Python methods (Compatibility issues)


############################################################
# Install Python Packages - To Read .h5ad file
############################################################

py_install("anndata", method = "auto", pip = TRUE)
py_install("scanpy")
py_install("numpy")
py_install("scipy")

# To check for AnnData packages
py_module_available("anndata")

############################################################
# Read .h5ad File Using Python
############################################################

ad <- import("anndata")
adata <- ad$read_h5ad("dataset/human_immune_cells.h5ad")

# Inspect AnnData object
adata

############################################################
# Extract Data From AnnData
############################################################

# Extract matrix
counts <- adata$X

# Convert to R sparse matrix
counts <- as(as.matrix(counts), "dgCMatrix")

# Transpose matrix (AnnData: Cells x Genes → Seurat: Genes x Cells)
counts <- t(counts)

# Extract metadata
cell_meta <- py_to_r(adata$obs)
gene_meta <- py_to_r(adata$var)


# To Check for Matrix Orientation
dim(counts)

# Assign row and column names
rownames(counts) <- gene_meta$gene_symbols
colnames(counts) <- rownames(cell_meta)

############################################################
# Create Seurat Object
############################################################

seurat_obj <- CreateSeuratObject(
  counts = counts,
  meta.data = cell_meta
)


############################################################
# Calculate QC Metrics (Mitochondrial %)
############################################################

# Identify mitochondrial genes
mito_genes <- as.character(gene_meta$gene_symbols[gene_meta$mito == TRUE])
mito_genes <- trimws(mito_genes)

# Calculate percent.mt
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
  seurat_obj,
  features = mito_genes
)


# Save backup
saveRDS(seurat_obj, "script/human_immune_seurat_backup.rds")


############################################################
# QC Visualization (Before Filtering)
############################################################

# Violin plots
VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

# QC Scatter Plots
plot1 <- FeatureScatter(
  seurat_obj,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)

plot2 <- FeatureScatter(
  seurat_obj,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
)

# Combine plots
plot1 + plot2

# Save Figure
png("figures/Fig2_QC_scatter.png", width = 1000, height = 500)
plot1 + plot2
dev.off()

############################################################
# Filter Cells (QC Thresholds)
############################################################

seurat_obj <- subset(
  seurat_obj,
  subset =
    nFeature_RNA > 200 &
    nFeature_RNA < 4500 &
    percent.mt < 12
)

############################################################
# QC Visualization (After Filtering)
############################################################

VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

############################################################
# Normalization & Feature Selection
############################################################

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst")

############################################################
# Scaling & PCA
############################################################

seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, npcs = 30)

############################################################
# PCA Diagnostics
############################################################

# Elbow Plot
ElbowPlot(seurat_obj)

png("figures/Fig3_PCA_Elbow.png", width = 700, height = 500)
ElbowPlot(seurat_obj)
dev.off()

DimHeatmap(seurat_obj, dims = 1:15)

png("figures/Fig3_PCA_Heatmap.png", width = 1000, height = 800)
DimHeatmap(seurat_obj, dims = 1:15)
dev.off()

# PCA Loadings
VizDimLoadings(seurat_obj, dims = 1:15, reduction = "pca")

png("figures/Fig3_PCA_Loadings.png", width = 1000, height = 1700)
VizDimLoadings(seurat_obj, dims = 1:15, reduction = "pca")
dev.off()

# PCA Scatter
DimPlot(seurat_obj, reduction = "pca", label = TRUE)

png("figures/Fig3_PCA_Scatter.png", width = 800, height = 600)
DimPlot(seurat_obj, reduction = "pca")
dev.off()


############################################################
# Clustering and UMAP
############################################################

# Different PCs
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:15)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)

# Different Resolutions
seurat_obj <- FindClusters(seurat_obj, resolution = 0.3)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.8)


# To Check number of clusters
length(unique(Idents(seurat_obj)))

# To Run UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:15)

# To UMAP Colored by Clusters Visualization and Save
DimPlot(seurat_obj, reduction = "umap", label = TRUE)

png("figures/Fig4_UMAP_clusters_dims1_15_res0_5.png", width = 650, height = 400)
DimPlot(seurat_obj, reduction = "umap",  label = TRUE)
dev.off()

# Find Markers for 9 Clusters (Dims 1:15, res = 0.5)
markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# To View top markers per clusters
head(markers)

# To Save Tables of markers
write.csv(markers, "tables/cluster_markers.csv", row.names = FALSE)

# TO visualize and save markers dotplot
png("figures/Fig5_marker_dotplot.png", width=4000, height=750)
DotPlot(seurat_obj, features = unique(markers$gene))
dev.off()

# To Extract top 10 Markers per cluster
library(dplyr)

top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

write.csv(top10, "tables/top10_markers.csv")


# Scale data if not already done- after selecting variable features
seurat_obj <- ScaleData(
  seurat_obj,
  features = rownames(seurat_obj)
)

# Heatmap Visualization
heatmap_plot <- DoHeatmap(
  seurat_obj,
  features = unique(top10$gene),
  group.by = "seurat_clusters"
) + 
  NoLegend()

heatmap_plot

png("figures/Fig5_marker_heatmap.png", width=6000, height=2000)
heatmap_plot
dev.off()

# To visualize dotplot markers
dot_plot <- DotPlot(
  seurat_obj,
  features = unique(top10$gene)
) + RotatedAxis()

dot_plot

png("figures/Fig5_top10_marker_dotplot.png", width=2000, height=600)
dot_plot
dev.off()

# To list Top 10 genes
top10 %>% filter(cluster == 0)

############################################################
# To Visualize Feature Plots
############################################################
# Cluster 0
png("figures/Fig5_marker_featureplot_cluster_0.png", width=1800, height=800)
FeaturePlot(seurat_obj, features = c("FCER1G", "MYL9", "CD151", "PARVB", "CALD1", "C12orf75" , "TREML1", "CD9", "FLNA", "CXCL5"))
dev.off()

# Cluster 1
png("figures/Fig5_marker_featureplot_cluster_1.png", width=1800, height=800)
FeaturePlot(seurat_obj, features = c("INKA1", "LGALS12", "IFRD1", "HEXIM2", "FAM110A", "FRMD3" , "EGLN3", "PTGIR", "FRMD4B", "TBPL1"))
dev.off()

# Cluster 2
png("figures/Fig5_marker_featureplot_cluster_2.png", width=1800, height=800)
FeaturePlot(seurat_obj, features = c("S100A9", "S100A8", "LYZ", "S100A6", "S100A4", "TMSB10" , "MYL9", "CLU", "CST3", "F13A1"))
dev.off()

# Cluster 3
png("figures/Fig5_marker_featureplot_cluster_3.png", width=1800, height=800)
FeaturePlot(seurat_obj, features = c("IFIT1B", "GYPB", "AHSP", "SLC4A1", "HBM", "EPB42" , "SELENBP1", "CA1", "IFI27", "BPGM"))
dev.off()

# Cluster 4
png("figures/Fig5_marker_featureplot_cluster_4.png", width=1800, height=800)
FeaturePlot(seurat_obj, features = c("AVP", "TMEM246", "C1QTNF4", "AC026369.3", "HTR1F", "AC011139.1" , "ADGRG6", "CRHBP", "CRYGD", "BSPRY"))
dev.off()

# Cluster 5
png("figures/Fig5_marker_featureplot_cluster_5.png", width=1800, height=800)
FeaturePlot(seurat_obj, features = c("LINC00861", "RORA", "SYNE2", "AAK1", "ANKRD44", "SLFN5" , "EML4", "CLEC2D", "PCSK7", "CD2"))
dev.off()

# Cluster 6
png("figures/Fig5_marker_featureplot_cluster_6.png", width=1800, height=800)
FeaturePlot(seurat_obj, features = c("DNTT", "PRSS2", "ACY3", "CLNK", "KIAA0087", "SLC2A5" , "VPREB1", "NPTX2", "MME", "SCN3A"))
dev.off()

# Cluster 7
png("figures/Fig5_marker_featureplot_cluster_7.png", width=1800, height=800)
FeaturePlot(seurat_obj, features = c("TPSAB1", "TPSB2", "SLC27A2", "CNRIP1", "DEPTOR", "CPA3" , "RHEX", "CSF2RB", "HPGDS", "FSCN1"))
dev.off()

# Cluster 8
png("figures/Fig5_marker_featureplot_cluster_8.png", width=1800, height=800)
FeaturePlot(seurat_obj, features = c("EPX", "CDK15", "PRG2", "ENPP3", "CLC", "TPSD1" , "FAM83F", "FOXJ1", "MS4A2", "HDC"))
dev.off

# To UMAP Colored by Clusters Visualization and Cell Types
p_clusters <- DimPlot(seurat_obj, reduction = "umap", label = TRUE)
p_clusters

############################################################
# Assign names to clusters
############################################################

# Rename Clusters
new.ids <- c(
  "Platelets",
  "Classical Monocytes (CD14⁺ Monocytes)",
  "Activated / Inflammatory Classical Monocytes",
  "Erythrocytes (Red Blood Cells)",
  "Doublets / Contaminating non-hematopoietic cells",
  "T cells",
  "Pre-B cells / Immature B cells",
  "Mast Cells",
  "Eosinophils"
)

names(new.ids) <- levels(seurat_obj)

seurat_obj <- RenameIdents(seurat_obj, new.ids)

# Save into metadata
seurat_obj$celltype <- Idents(seurat_obj)

# To visualize UMAP Cell Types
p_celltypes <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "celltype",
  label = TRUE,
  label.size = 4,
  pt.size = 0.4
)
p_celltypes

# To combine 2 plots
library(patchwork)

p_clusters + p_celltypes

png("figures/Fig4_UMAP_Celltypes.png", width=1500, height=500)
p_clusters + p_celltypes
dev.off()

############################################################
# Differential Expression 
############################################################

# Classical Monocytes C1 vs Activated / Inflammatory Classical Monocytes C2

Idents(seurat_obj) <- "seurat_clusters"

mono_DE <- FindMarkers(
  seurat_obj,
  ident.1 = "1",
  ident.2 = "2",
  logfc.threshold = 0.25,
  min.pct = 0.25
)

head(mono_DE)

# To save DE Tables
write.csv(mono_DE, "tables/DE_Classical_Monocytes(C1)_vs_Activated Monocytes(C2).csv")

# To visualize EnhancedVolcano Plot
library(EnhancedVolcano)

volcanoplot_de <- EnhancedVolcano(
  mono_DE,
  lab = rownames(mono_DE),
  x = 'avg_log2FC',
  y = 'p_val_adj',
  pCutoff = 0.05,
  FCcutoff = 0.5,
  title = 'Classical Monocytes vs Activated Monocytes'
)
volcanoplot_de

png("figures/Fig6_DE_results.png", width=1700, height=1000)
volcanoplot_de
dev.off()


# To Extract Top 10 Differentially Expressed Genes
library(dplyr)

top_classical <- mono_DE %>%
  arrange(desc(avg_log2FC)) %>%
  head(5)

top_classical

top_activated <- mono_DE %>%
  arrange(avg_log2FC) %>%
  head(5)

top_activated

genes_to_plot <- c(
  rownames(top_classical)[1:3],
  rownames(top_activated)[1:3]
)

genes_to_plot

# To visualize feature plot for top DE genes
FeaturePlot(seurat_obj, features = genes_to_plot, order = TRUE, pt.size = 0.4, min.cutoff = "q10", max.cutoff = "q90")

# To Save Feature plot
png("figures/Fig6_DE_featureplot.png", width=800, height=600)
FeaturePlot(seurat_obj, features = genes_to_plot, order = TRUE, pt.size = 0.4, min.cutoff = "q10", max.cutoff = "q90")
dev.off()

# To visualize DE Volcano Plot
library(ggplot2)

mono_DE$gene <- rownames(mono_DE)

ggplot(mono_DE, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  ggtitle("Differential Expression: Classical vs Activated Monocytes")

png("figures/Fig6_DE_volcanoplot.png", width=800, height=600)
ggplot(mono_DE, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  ggtitle("Differential Expression: Classical vs Activated Monocytes")
dev.off()


# To visualized detailed Volcano plot 
library(dplyr)
library(tibble)
mono_DE <- mono_DE %>%
  rownames_to_column(var = "gene_name")

mono_DE <- mono_DE %>%
  mutate(
    significance = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0.25  ~ "Up in Classical",
      p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "Up in Activated",
      TRUE ~ "Not Significant"
    )
  )

top_genes <- mono_DE %>%
  filter(significance != "Not Significant") %>%
  arrange(desc(abs(avg_log2FC))) %>%
  slice_head(n = 10)

library(ggplot2)
library(ggrepel)

volcano_plot <- ggplot(mono_DE,
                       aes(x = avg_log2FC,
                           y = -log10(p_val_adj))) +
  
  # Points
  geom_point(aes(color = significance),
             alpha = 0.7,
             size = 1.8) +
  
  # Color scheme
  scale_color_manual(values = c(
    "Up in Classical" = "red",
    "Up in Activated" = "blue",
    "Not Significant" = "grey"
  )) +
  
  # Threshold lines
  geom_vline(xintercept = c(-0.25, 0.25),
             linetype = "dashed") +
  
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  
  # Labels for top genes
  geom_text_repel(
    data = top_genes,
    aes(label = gene),
    size = 3,
    max.overlaps = 20
  ) +
  
  theme_minimal() +
  labs(
    title = "Differential Expression: Classical vs Activated Monocytes",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Expression"
  )

volcano_plot

png("figures/Fig6_DE_Coloured_volcanoplot.png", width=800, height=600)
volcano_plot
dev.off()



# To Save Session Info
sink("session_info/sessionInfo.txt")
sessionInfo()
sink()
