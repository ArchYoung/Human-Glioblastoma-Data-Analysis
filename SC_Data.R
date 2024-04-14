# custom plot size function
fig <- function(width, heigth){
  options(repr.plot.width = width, repr.plot.height = heigth)
}

# install seurat
# remotes::install_github("satijalab/seurat", ref = "release/4.0.0")

# load required libararies
library(Seurat) # scRNAseq analysis toolset
suppressPackageStartupMessages(library(dplyr)) # data manipulation
library(patchwork) # data visualization in grids
library(ggplot2)
library(scales)
# Load the dataset
# feature / cell matrix (filtered)
sc_data <- Read10X(data.dir = "C:/Users/YangHanming/Desktop/课project/archive")
cat("Original dimension of the count matrix:", dim(sc_data))

# Initialize the Seurat object with the raw (non-normalized data).
sc <- CreateSeuratObject(counts = sc_data, project = "GBM", min.cells = 3, min.features = 200)
sc

# Look at first three genes in the first thirty cells
sc_data[1:3,1:30]

# calculate mitochondrial percentage and asign this to the "meta.data" dataframe
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")

# check the dataframe with the new column, cells with the most mitochondrial contamination first
# bc this is a filtered dataset, it should be free of contamination
sc@meta.data %>%
  arrange(desc(percent.mt)) %>%
  head(n = 10)

# if you do need to filter, use the base subset() function
sc <- subset(sc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 0.05)

# custom plot size
fig(16,9)

# visualize QC metrics
VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.5)

cat("The dataset contains cells with", min(sc@meta.data$nFeature_RNA), "up to", max(sc@meta.data$nFeature_RNA), 
    "features (genes) and", min(sc@meta.data$nCount_RNA), "up to", max(sc@meta.data$nCount_RNA), "neuro-specific target gene copies.")

# this method divides the feature expression measures for each cell by the total expression of the cell
# then then result is multiplied by 10'000 (scale.factor) and log transformed (normalization.method)
sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 10000)

# return the 500 most variable features
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 500)

# Identify the 10 most variable genes
top10_genes <- head(VariableFeatures(sc), 10)

cat("The top 10 most variable genes among glioblastoma cells are:\n", top10_genes)


# plot variable features with and without labels of the top 10 genes
plot1 <- VariableFeaturePlot(sc)
plot2 <- LabelPoints(plot = plot1, points = top10_genes, repel = TRUE, xnudge = -0.8, ynudge = 0)
plot1 + plot2


all.genes <- rownames(sc)

sc <- ScaleData(sc, features = all.genes, verbose = FALSE)
sc <- ScaleData(sc, vars.to.regress = "percent.mt")

sc = RunPCA(sc,features = VariableFeatures(object = sc),npcs = 50, verbose = FALSE)
print(sc[["pca"]], dims = 1:5, nfeatures = 5)


VizDimLoadings(sc, dims = 1:6,
               ncol = 3, 
               reduction = "pca")


# get the percent variance explained by the individual principles components
mat <- GetAssayData(sc, assay = "RNA", slot = "scale.data")
pca <- sc[["pca"]]

# Get the total variance:
total_variance <- sum(matrixStats::rowVars(mat))

eigValues <- (pca@stdev)^2  ## EigenValues
varExplained <- eigValues / total_variance

# calculate the % variance explained by the first two PCs
cat("The first two principle components explain", 
    scales::label_percent()(round(sum(varExplained[c(1,2)]), 4)), "of the total variance in the dataset.")

# custom plot size
fig(14, 10)

# store varEplained in a dataframe alongside the PC number
var_explained_df <- tibble(pc_number = 1:length(varExplained), var_explained = cumsum(varExplained))

# load ggplot2 package
library(ggplot2)

# plot this dataframe
ggplot(var_explained_df, aes(x = pc_number, y = var_explained)) +
  geom_point(size = 2) +
  geom_line(size = 0.8) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Cumulative Variance Explained by the first 50 PCs",
       x = "PC",
       y = "Cumulative Variance explained") +
  theme_bw() +
  theme(plot.title = element_text(size = 22),
        plot.subtitle = element_text(size = 16),
        axis.text.x= element_text(size = 15),
        axis.text.y= element_text(size = 15), 
        axis.title=element_text(size = 18))



DimPlot(sc, reduction = "pca",dims = c(1,2))

DimHeatmap(sc, 
           dims = 1:15,  #几个PCs
           cells = 500, #多少cell
           ncol = 3, #图形展示几列
           balanced = TRUE)

sc <- JackStraw(sc, num.replicate = 300, prop.freq = 0.05)
sc <- ScoreJackStraw(sc, dims = 1:20)

JackStrawPlot(sc, dims = 1:20)


ElbowPlot(sc, ndims = 30)

sc <- FindNeighbors(sc, dims = 1:18)
sc <- FindClusters(sc, algorithm = 1, resolution = 0.2)

# check cluster identity of the first 10 cells
head(Idents(sc), 10)


# custom fig dim
fig(14, 10)

# PCA plot again with clusters encoded as colors
DimPlot(sc, reduction = "pca", dims = c(1,2))




sc <- RunUMAP(sc, dims = 1:18)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(sc, reduction = "umap")


sc <- RunTSNE(sc, dims = 1:18)
head(sc@reductions$tsne@cell.embeddings)
DimPlot(sc, reduction = "tsne")

plot5<-DimPlot(sc, reduction = "umap",label = TRUE)
plot6<-DimPlot(sc, reduction = "tsne",label = TRUE)
plot5 + plot6



cluster.markers <- FindMarkers(sc,ident.1 = 2,ident.2 = c(0,4), min.pct = 0.3)
head(cluster.markers, n = 5)

sc.markers <- FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
sc.markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)

FeaturePlot(sc, features = c("MAG","COL4A1","PDGFRB","FLT4"))

top10 <- sc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
DoHeatmap(sc, features = top10$gene)


library(SingleR)
load("C:/Users/YangHanming/Desktop/课project/HumanPrimaryCellAtlasData.Rdata")
load("C:/Users/YangHanming/Desktop/课project/BlueprintEncodeData.Rdata")
refdata = hpca.se
testdata = GetAssayData(sc,slot = "data")
clusters = sc@meta.data$seurat_clusters

cellpred = SingleR(test = testdata, ref = list(HPCA = hpca.se), labels = list(hpca.se$label.main),
                   method = "cluster",clusters = clusters,
                   assay.type.test = "logcounts", assay.type.ref = "logcounts")


celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = FALSE)    ###制作细胞类型的注释文件
sc@meta.data$singleR = celltype[match(clusters,celltype$ClusterID),'celltype']

p1 = DimPlot(sc, group.by="singleR", label=T, label.size=5, reduction='tsne')

p2 = DimPlot(sc, group.by="singleR", label=T, label.size=5, reduction='umap')

p3 <- p1+p2+ plot_layout(guides = 'collect')


print(plotScoreHeatmap(cellpred))


library(Seurat)
library(SeuratObject)
library(tidyseurat)


seurat_object %>%
  add_count(sc) %>%
  filter(PC_1 > 0 & n > 40)



pbmc_small_UMAP <-
  sc %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:15, n.components = 3L)

sc <- RunUMAP(sc,reduction = "pca", dims = 1:15, n.components = 3L)

p = sc %>%
plot_ly(
        x = ~`umap_1`,
        y = ~`umap_2`,
        z = ~`umap_3`,
        color = sc$singleR,
        size = 1,
        symbol = sc$singleR
)

factor(sc@meta.data$singleR)

sc_1 = sc_1$seurat_clusters
sc_1$RNA_snn_res.0.2 = factor(sc_1@meta.data$singleR)
sc_1@misc

sc_1@meta.data$seurat_clusters = factor(sc_1@meta.data$singleR)
sc_1@active.ident = factor(sc_1@meta.data$singleR)
names(sc_1@active.ident) = names(sc@active.ident)



sc.markers <- FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top10 <- sc_1.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
DoHeatmap(sc_1, features = top10$gene)


