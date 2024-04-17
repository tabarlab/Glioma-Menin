### Required Libraries
library(dplyr)
library(Seurat)
library(ggplot2)
library(harmony)
library(miloDE)
library(miloR)
library(scWGCNA)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(uwot)
library(scran)
library(future)
plan(multicore, workers = 10)
options(future.globals.maxSize = 50 * 1024 ^ 3) # for 50 Gb RAM

######## Read Seurat Object
RNA_all <- readRDS('./scRNA_YY_combined.rds')

######## Integrate Seurat Object Across Cell Lines
RNA_all <- NormalizeData(RNA_all)
RNA_all <- FindVariableFeatures(RNA_all, selection.method = "vst", nfeatures = 2000)
RNA_all <- ScaleData(RNA_all, verbose = FALSE)
RNA_all <- RunPCA(RNA_all, npcs = 30, verbose = FALSE)
RNA_all <- RunHarmony(RNA_all, c('cell.line'))
RNA_all <- RunUMAP(RNA_all, reduction = "harmony", dims = 1:30)
RNA_all <- FindNeighbors(RNA_all, reduction = "harmony", dims = 1:30)
RNA_all <- FindClusters(RNA_all, resolution = 0.3)

######## Convert Seurat Object to a SingleCellExperiment Object
sce <- as.SingleCellExperiment(RNA_all, assay = 'RNA')

sce <- logNormCounts(sce)

# Get Harmonized UMAP Coordinates
set.seed(32)
umaps <- as.data.frame(uwot::umap(reducedDim(sce , "HARMONY")))
# Store UMAP Coordinates in SCE Object for visualization
reducedDim(sce , "UMAP") = umaps

######## Extended Figure 6D
cell_line_colormap <- c('#ffc0cb', '#fdb509', '#a0522d', '#da70d6', '#40e0d0', '#193c66', '#008b00', '#a0522d', '#193c66')

umaps = cbind(as.data.frame(colData(sce)) , reducedDim(sce , "UMAP"))

d <- ggplot(umaps , aes(x = V1 , y = V2 , col = new.cell.line)) +
            geom_point(size = .1) + 
            scale_color_manual(values = cell_line_colormap) +
            theme_classic() + 
            labs(x = "UMAP-1", y = "UMAP-2") + theme(aspect.ratio=1)


RNA_2075_NSG <- CreateSeuratObject(counts = IN, assay = "RNA", min.cells = 3, min.features = 200)
RNA_2075_NSG[["percent.mt"]] <- PercentageFeatureSet(RNA_2075_NSG, pattern = "^MT-")
RNA_2075_NSG[["percent.rp"]] <- PercentageFeatureSet(RNA_2075_NSG, pattern="^(RPS|RPL)")

######## Adding HTO data as a new assay
RNA_2075_SCG_HTO_data <- LoadH5Seurat("./scRNA_YY/2075_SGC/YY-2075_SGC.h5seurat")
RNA_2075_SCG[["HTO"]] <- RNA_2075_SCG_HTO_data[['RNA']]

RNA_2075_NSG_HTO_data <- LoadH5Seurat("./scRNA_YY/2075_NSG/YY-2075_NSG.h5seurat")
RNA_2075_NSG[["HTO"]] <- RNA_2075_NSG_HTO_data[['RNA']]

######## Filter Low Quality Cells
VlnPlot(RNA_2075_SCG, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 4)
FeatureScatter(RNA_2075_SCG, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
RNA_2075_SCG <- subset(RNA_2075_SCG, subset = nFeature_RNA > 200 & nFeature_RNA < 7500  & percent.mt < 25)

VlnPlot(RNA_2075_NSG, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 4)
FeatureScatter(RNA_2075_NSG, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
RNA_2075_NSG <- subset(RNA_2075_NSG, subset = nFeature_RNA > 200 & nFeature_RNA < 7500  & percent.mt < 25)

######## Demultiplex cells based on HTO enrichment
RNA_2075_SCG <- HTODemux(RNA_2075_SCG, assay = "HTO", positive.quantile = 0.99)
RNA_2075_NSG <- HTODemux(RNA_2075_NSG, assay = "HTO", positive.quantile = 0.99)

######## Remove doublets and negative cells
RNA_2075_NSG <- subset(RNA_2075_NSG, idents = c("Negative", "Doublet"), invert = TRUE)
RNA_2075_SCG <- subset(RNA_2075_SCG, idents = c("Negative", "Doublet"), invert = TRUE)

######## Combine Seurat Objects and Subset Cell lines for Downstream Analysis 
SCG_names <- c('MSK-20', 'MSK-23', 'MSK-09', 'MSK-19', 'MSK-28')
NSG_names <- c('MSK-28T1-d4', 'MSK-19KO', 'MSK-18', 'MSK-04', 'MSK-28T2-d5')
names(NSG_names) <- levels(RNA_2075_NSG)
names(SCG_names) <- levels(RNA_2075_SCG)

RNA_2075_SCG <- RenameIdents(RNA_2075_SCG, SCG_names)
RNA_2075_NSG <- RenameIdents(RNA_2075_NSG, NSG_names)

RNA_2075_SCG$cell.line <- Idents(RNA_2075_SCG)
RNA_2075_NSG$cell.line <- Idents(RNA_2075_NSG)

RNA_2075_SCG$batch <- '1'
RNA_2075_NSG$batch <- '2'

RNA_all <- merge(x = RNA_2075_SCG,y=(list(RNA_2075_NSG)),
                 add.cell.ids = c("SGC", "NSG"))

RNA_all <- subset(RNA_all, idents = c("MSK-20", "MSK-23", "MSK-18", "MSK-19", "MSK-28", "MSK-19KO", "MSK-04"))

######## Add Menin Sensitivity Status to Metadata
levels(RNA_all)
# 'MSK-20' 'MSK-23' 'MSK-19' 'MSK-28' 'MSK-19KO' 'MSK-18' 'MSK-04'
# Sensitive: MSK-19, MSK-28, MSK-23, MSK-20. 
# Non-sensitive: MSK-18, MSK-04
MEN1.sensitivity <- c('Sensitive', 'Sensitive', 'Sensitive', 'Sensitive', 'KO', 'Non-Sensitive', 'Non-Sensitive')
names(MEN1.sensitivity) <- levels(RNA_all)
RNA_all <- RenameIdents(RNA_all, MEN1.sensitivity)
RNA_all$MEN1.sensitivity <- Idents(RNA_all)

######## Normalization, Scaling, PCA, UMAP Generation and Clustering ########

RNA_all <- NormalizeData(RNA_all)
RNA_all <- FindVariableFeatures(RNA_all, selection.method = "vst", nfeatures = 2000)
RNA_all <- ScaleData(RNA_all, verbose = FALSE)
RNA_all <- RunPCA(RNA_all, npcs = 30, verbose = FALSE)
RNA_all <- RunUMAP(RNA_all, reduction = "pca", dims = 1:30)
RNA_all <- FindNeighbors(RNA_all, reduction = "pca", dims = 1:30)
RNA_all <- FindClusters(RNA_all, resolution = 0.5)

saveRDS(RNA_all, file="scRNA_YY_combined.rds")

######## UMAP for Extended Figure 6A

p1 <- DimPlot(RNA_all, reduction = "umap", group.by="cell.line")

ggsave(file="./scRNA_merged_UMAP_cell_lines.pdf", plot=p1, width=10, height=10)

