#BINF 6110
#Assignment 4: Single cell Transcriptomics
#Submitted by Iroayo Toki
#April 13th, 2026
#---------------------
#Packages Used
library(dplyr)
library(Seurat)
library(celldex)
library(SingleR)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(enrichplot)
library(clusterProfiler)
library(org.Mm.eg.db)
library(CellChat)
#Load RDS

Seurat_obj <- readRDS("seurat_ass4.rds")

#Filtering by mitochondrial DNA percentage
Seurat_obj[["percent.mt"]] <- PercentageFeatureSet(Seurat_obj, pattern = "^mt-")
table(Seurat_obj[["percent.mt"]])

# Visualize QC metrics as a violin plot
VlnPlot(Seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

#Subsetting data with QC metrics
Seurat_obj<- subset(Seurat_obj, subset = nFeature_RNA > 200 & percent.mt < 20)
VlnPlot(Seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

#Normalization of the data and scaling
Seurat_obj <- NormalizeData(Seurat_obj, normalization.method = "LogNormalize")
Seurat_obj <- FindVariableFeatures(Seurat_obj, selection.method = "vst", nfeatures = 2000)
#Scale
Seurat_obj<- ScaleData(Seurat_obj)

#Running PCA to produce principal components that can be used to cluster our cells
Seurat_obj <- RunPCA(Seurat_obj, features = VariableFeatures(object = Seurat_obj))

# Visualize the dimensionality of the dataset
ElbowPlot(Seurat_obj)

# Clustering
Seurat_obj <- FindNeighbors(Seurat_obj, dims = 1:18)
Seurat_obj <- FindClusters(Seurat_obj, resolution = 0.5)

# UMAP with clusters displayed
Seurat_obj <- RunUMAP(Seurat_obj, dims = 1:18)
DimPlot(Seurat_obj, reduction = "umap", label = TRUE)


#Automatic annotation
#Calculate average expression
cluster_avg <- AverageExpression(
  Seurat_obj,
  group.by = "seurat_clusters"
)
#Load References
ref_imm <- ImmGenData()
ref_all <- MouseRNAseqData()

#Annotate
pred_imm <- SingleR(
  test = cluster_avg$RNA,
  ref = ref_imm,
  labels = ref_imm$label.main
)

pred_all <- SingleR(
  test = cluster_avg$RNA,
  ref = ref_all,
  labels = ref_all$label.main
)

#create Df with with labels
cluster_ids <- colnames(cluster_avg$RNA)

cluster_labels <- data.frame(
  cluster = cluster_ids,
  immune_label = pred_imm$labels,
  tissue_label = pred_all$labels
)

#Manual annotation to confirm
#Create list with upregulated markers present in at least 25% of cells and have at least 0.25 log2 fold change to keep relevant markers.
markers <- FindAllMarkers(
  Seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
saveRDS(markers, "SCseqmarkers.rds")

# Alternatively Load in markers from rds
markers <- readRDS("SCseqmarkers.rds")

#Top 5 markers per cluster based on LFC filtered P value
Top_5_markers <- markers %>%
  filter(p_val_adj < 0.05) %>%      
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5)      

grep("Agr", rownames(Seurat_obj), value = TRUE)
 
#Checking manually for each cluster
head(rownames(markers[markers$cluster == 22, ]), n = 5)
FeaturePlot(Seurat_obj, features = c("Mdga2"))


#After reconciling Automatic and manual annotation 
#Rename clusters
new.cluster.ids <- c("Neurons", "Olfactory sensory Neurons", "Macrophages", "Basal Epithelial cells(Respiratory)", "B cells", "Endothileal cells", "Developing Neurons", "Natural killer cells", "Fibroblasts", "Epithelial cells", "Neutrophils", "Monocytes", "Epithelial cells", "Respiratory epithelial cells", "Epithelial cells", "Myeloid Leukocytes(Monocytes Macrophages and DC)", "Epithelial duct cells", "Respiratory Epithelial cells", "Monocytes", "Vascular smooth muscle cells", "Neutrophils", "Proliferating cells", "Secretory epithelial cells", "Epithelial cells", "Immature neurons", "Tear Secreting epithelial cells", "Microvillar Epithelial cells", "Oligodendrocytes", "Osteoblasts", "B cells", "Ciliated Epithelial cells", "Fibroblast", "Proliferating cells", "Epithelial cells",  "Natural killer cells", "Epithelial cells(Hillock)",   "Progenitor cells")
names(new.cluster.ids) <- levels(Seurat_obj)
Seurat_obj <- RenameIdents(Seurat_obj, new.cluster.ids)
#Plot with new labels
DimPlot(Seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
Seurat_obj$Cell_types <- Idents(Seurat_obj)
head(Seurat_obj@meta.data)
table(Seurat_obj$organ_custom)


#Cell type of interest: Secretory epithelial cells
# Feature plots of identifier genes from the markers dataframe
head(rownames(markers[markers$cluster == 25, ]), n = 5)

FeaturePlot(Seurat_obj, features = c("Cyp2g1",  "Sec14l3"))
FeaturePlot(Seurat_obj, features = c("Cyp1a2", "Muc2"))


#Differential Expression
#Pseudo bulk by disease state, Sample id and cell type
pseudo_data<- AggregateExpression(Seurat_obj, assays = "RNA", return.seurat = T, group.by = c("disease__ontology_label", "mouse_id", "Cell_types"))

#Create variable that combines cell type and disease state
pseudo_data$celltype.disease <- paste(pseudo_data$Cell_types, pseudo_data$disease__ontology_label , sep = "_")

Idents(pseudo_data) <- "celltype.disease"
table(pseudo_data$celltype.disease)
#Comparing samples within The secretory epithelial cells clusters
bulk.secretion.de <- FindMarkers(object = pseudo_data, 
                            ident.1 = "Secretory epithelial cells_influenza", 
                            ident.2 = "Secretory epithelial cells_normal",
                            test.use = "DESeq2")
head(bulk.secretion.de, n = 30)


#Add combined variable to main data 
Seurat_obj$celltype.disease <- paste(Seurat_obj$Cell_types, Seurat_obj$disease__ontology_label , sep = "_")

#Visualization of genes highly expressed in infected mice
VlnPlot(
  subset(Seurat_obj, Cell_types == "Secretory epithelial cells"),
  features = c('mt-Atp6','mt-Atp8','Muc16'),
  group.by = "disease__ontology_label")

FeaturePlot(Seurat_obj, features = c('mt-Atp6','mt-Atp8','Muc16'))

#Volcano plot 
bulk.secretion.de$significant <- ifelse(bulk.secretion.de$p_val_adj < 0.05, ifelse(bulk.secretion.de$avg_log2FC > 0, "Up", "Down"), "Not Sig")
bulk.secretion.de <- na.omit(bulk.secretion.de)
Filtered_df <-  bulk.secretion.de %>% filter(p_val_adj< 0.05)
#Plot with significant genes labelled
ggplot(bulk.secretion.de, aes(x = avg_log2FC, y = -log10(p_val), color = significant)) +
  geom_point() +
  scale_color_manual(values = c("Down" = "blue", "Not Sig" = "gray", "Up" = "red")) +
  labs(x = "Log2 Fold Change", y = "-Log10 p-value", 
       title = "Volcano Plot: Influenza vs normal in Secretory epithelial cells") +
  theme(legend.position = "right")  +
  geom_text_repel(data = Filtered_df,
                  aes(label = rownames(Filtered_df)),
                  size = 3)
#Heatmap
#Lets use significant genes to create a heatmap across days of the infection to see how expression changes in our top expressed genes.
#Create subset object
Idents(Seurat_obj) <- "Cell_types"
cluster_cells <- WhichCells(Seurat_obj, idents = "Secretory epithelial cells")
subset_obj <- subset(Seurat_obj, cells = cluster_cells)
unique(subset_obj$time)

#Reorder day variable
subset_obj$time <- factor(subset_obj$time,
                         levels = c("Naive", "D02", "D05", "D08", "D14"))

# Average expression
avg_expr <- AverageExpression(
  subset_obj,
  features = Sig_genes,
  group.by = "time"
)$RNA

# Reorder columns to match factor levels
avg_expr <- avg_expr[, levels(subset_obj$time)]

# Build annotation
annotation_col <- data.frame(day = factor(colnames(avg_expr),
                                          levels = levels(subset_obj$time)))
rownames(annotation_col) <- colnames(avg_expr)

# Plot
pheatmap(avg_expr,
         annotation_col = annotation_col, cluster_cols = F)

#GO Overrepresentation analysis for biological process and Cell com
all_genes <- rownames(bulk.secretion.de)

ego_bp <- enrichGO(gene = Sig_genes,
                   universe = all_genes,
                   OrgDb = org.Mm.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,
                   readable = TRUE, keyType = "SYMBOL")
ego_cc <- enrichGO(gene = Sig_genes,
                   universe = all_genes,
                   OrgDb = org.Mm.eg.db,
                   ont = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,
                   readable = TRUE, keyType = "SYMBOL")
ego_bp <- simplify(ego_bp)
ego_cc <- simplify(ego_cc)

head(as.data.frame(ego_bp))
head(as.data.frame(ego_cc))
# Dot plots
dotplot(ego_bp, showCategory = 15, title = "GO Biological Process")

dotplot(ego_cc, showCategory = 15, title = "GO Cellular component")


# Bar plot
barplot(ego_bp, showCategory = 15, title = "GO Biological Process")

barplot(ego_cc, showCategory = 15, title = "GO Cellular component")

#Cell-cell communication between the secretory epithelial cells and other cell types
#Downsample for cellchat (500 per cluster)
Seurat_obj_cellchat <- subset(
  Seurat_obj,
  cells = unlist(
    tapply(Cells(Seurat_obj), Seurat_obj$Cell_types,
           function(x) sample(x, min(500, length(x))))
  )
)

# Prepare input
data.input <- GetAssayData(Seurat_obj_cellchat, layer = "data")
meta <- Seurat_obj_cellchat@meta.data

# Create CellChat object
cellchat <- createCellChat(object = data.input,
                           meta = meta,
                           group.by = "Cell_types")


# Load ligand-receptor database
CellChatDB <- CellChatDB.mouse  
cellchat@DB <- CellChatDB


# Subset relevant data
cellchat <- setIdent(cellchat, ident.use = "Cell_types")
cellchat <- subsetData(cellchat)

# Identify overexpressed genes + interactions
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
Backup <- cellchat
class(cellchat)
# Compute communication probabilities
cellchat <- computeCommunProb(cellchat, nboot = 50)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)


#Aggregate network and visualize

cellchat <- aggregateNet(cellchat)

# global interaction strength plot
netVisual_circle(cellchat@net$weight)

# bubble plot: interactions targeting secretory epithelial cells
#Too noisy
netVisual_bubble(cellchat,
                 targets.use = "Secretory epithelial cells",
                 remove.isolate = TRUE)
#Filter out top intersctions
celltype <- "Secretory epithelial cells"
df.net <- subsetCommunication(cellchat)

# filter cell type as sender OR receiver
df.sub <- df.net[
  df.net$source == celltype | df.net$target == celltype,
]
df.top <- df.sub[order(df.sub$prob, decreasing = TRUE), ][1:5, ]

#Plot
netVisual_bubble(
  cellchat,
  sources.use = unique(df.top$source),
  targets.use = celltype,
  signaling = unique(df.top$pathway_name)
)
