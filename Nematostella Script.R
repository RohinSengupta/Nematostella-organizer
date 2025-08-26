#Loading required packages
library(dplyr)
library(Seurat)
library(devtools)
library(readxl)
library(pals)

#Reading filtered feature matrix
matrix.mtx <- Read10X(data.dir = "//primus.img.local/data/29_lab/Rohin/Nematostella single cell/Cole2024 Nematostella single cell/SC Libraries/24hr gastrula")

#For gene conversion
genes = read_excel(path = "//primus.img.local/data/29_lab/Rohin/Nematostella single cell/Cole2024 Nematostella single cell/Cole 2024 Annotation.xlsx",
                   sheet = 'genes.LUT') #from supplementary materials
genes <- as.data.frame(genes)

#generate some gene lists for filtering:
mito.genes <- grep(pattern = "mitochondrial", genes$annotation_notes)
mitochondria = genes$gene_short_name[mito.genes]

rownames(matrix.mtx) <- genes$gene.name
#Start of actual code
FMatrix <- CreateSeuratObject(counts = matrix.mtx, project = "24hr Gastrula", min.cells = 1, min.features = 1)
FMatrix

# Visualize QC metrics as a violin plot
VlnPlot(FMatrix, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(FMatrix, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
FMatrix <- NormalizeData(FMatrix, normalization.method = "LogNormalize", scale.factor = 5000)
FMatrix <- FindVariableFeatures(FMatrix, selection.method = "vst", nfeatures = 1000)

#Scale the data
all.genes <- rownames(FMatrix)
FMatrix <- ScaleData(FMatrix, features = all.genes)

#Linear Dimensional Reduction
FMatrix <- RunPCA(FMatrix, features = all.genes, ndims = 30)

#Clustering of Dataset
FMatrix <- FindNeighbors(FMatrix, dims = 1:30)
FMatrix <- FindClusters(FMatrix, resolution = 3.4)

#Non-Linear Dimensional Reduction (UMAP)
FMatrix <- RunUMAP(FMatrix, dims = 1:30, min.dist = 0.6, n.neighbors = 18)
DimPlot(FMatrix, reduction = "umap", label = TRUE, pt.size = 2, order = TRUE)

FeaturePlot(FMatrix, features = c("CHRD-like-1", "FoxA", "DPP2-like-1", "WntA", "Wnt2", "Wnt3", "Wnt4", "FoxQ2a", "OtxC", "Brachyury", "Six3-6", "SnailA", "SnailB", "Fgf8-17-like", "FGF8A", "FoxB", "DMBX1A"), pt.size = 3, order = TRUE)

#Rename Clusters
new.cluster.ids <- c("Aboral Ectoderm", "Aboral Ectoderm", "Oral Ectoderm", "Unidentified 1", "Aboral Ectoderm", "Aboral Ectoderm", "Aboral Ectoderm", "Oral Ectoderm", "Aboral Ectoderm", "Oral Ectoderm", "Aboral Ectoderm", "Oral Ectoderm", "Aboral Ectoderm", "Aboral Ectoderm", "Endomesoderm 1", "Oral Ectoderm", "Blastopore Lip", "Oral Ectoderm", "Aboral Ectoderm", "Endomesoderm 2", "Unidentified 4", "Aboral Ectoderm", "Blastopore Lip", "Unidentified 3", "Unidentified 2", "Unidentified 2", "Elav+ Aboral Ectoderm", "Oral Ectoderm", "Aboral Ectoderm", "Unidentified 2", "Unidentified 1", "Unidentified 2", "Endomesoderm 2", "Unidentified 2", "Aboral Ectoderm", "Aboral Ectoderm", "Unidentified 4")
names(new.cluster.ids) <- levels(FMatrix)
FMatrix <- RenameIdents(FMatrix, new.cluster.ids)

saveRDS(FMatrix,'nv_late_gastrula_v4.rds')
