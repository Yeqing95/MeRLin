setwd("~/Documents/Wistar/Haiyin/MeRLin/Data/Original data/")
library(Seurat)
library(data.table)
library(tibble)
library(DoubletFinder)
library(ggplot2)
library(dplyr)
library(vegan)

`%ni%` <- Negate(`%in%`)

#### CTRL ####
CTRL <- Read10X_h5(filename = "In-vitro/CTRL/filtered_feature_bc_matrix.h5")
CTRL <- CreateSeuratObject(counts = CTRL)
CTRL <- PercentageFeatureSet(CTRL, pattern = "^MT-", col.name = "percent_mito")

CTRL$Experiment <- "In-vitro"
CTRL$Condition <- "CTRL"

CTRL_bc <- fread(input = "In-vitro/CTRL/ctrl_bc.txt", header = T) %>% column_to_rownames(var = "CB")
CTRL <- AddMetaData(object = CTRL, metadata = CTRL_bc, col.name = "Barcode")

VlnPlot(object = CTRL, layer = "counts", features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3, pt.size = 0.1)

nCount_cutoff <- 50000
nFeature_cutoff <- 2000
pMito_cutoff <- 20

CTRL$QC <- "PASS"
CTRL$QC[which(CTRL$nCount_RNA > nCount_cutoff)] <- "High_nCount"
CTRL$QC[which(CTRL$nFeature_RNA < nFeature_cutoff)] <- "Low_nFeature"
CTRL$QC[which(CTRL$percent_mito > pMito_cutoff)] <- "High_pMito"
CTRL$QC[which(CTRL$nCount_RNA > nCount_cutoff & CTRL$nFeature_RNA < nFeature_cutoff)] <- "High_nCount; Low_nFeature"
CTRL$QC[which(CTRL$nCount_RNA > nCount_cutoff & CTRL$percent_mito > pMito_cutoff)] <- "High_nCount; High_pMito"
CTRL$QC[which(CTRL$nFeature_RNA < nFeature_cutoff & CTRL$percent_mito > pMito_cutoff)] <- "Low_nFeature; High_pMito"
CTRL$QC[which(CTRL$nCount_RNA > nCount_cutoff & CTRL$nFeature_RNA < nFeature_cutoff & CTRL$percent_mito > pMito_cutoff)] <- "High_nCount; Low_nFeature; High_pMito"

temp <- CTRL
temp <- NormalizeData(temp, verbose = FALSE)
temp <- FindVariableFeatures(temp, verbose = FALSE)
temp <- ScaleData(temp, verbose = FALSE)
temp <- RunPCA(temp, npcs = 20, verbose = FALSE)
temp <- RunUMAP(temp, dims = 1:10, verbose = FALSE)
DimPlot(object = temp, group.by = "QC")

temp <- FindNeighbors(temp)
temp <- FindClusters(temp, resolution = 0.3)
DimPlot(object = temp, group.by = "seurat_clusters", label = T, label.size = 5) + NoLegend()

sweep.res <- paramSweep(temp, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
best_pK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
nExp_poi <- round(0.075 * nrow(temp@meta.data))
temp <- doubletFinder(temp, PCs = 1:10, pN = 0.25, pK = best_pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(object = temp, group.by = "DF.classifications_0.25_0.19_1115")

keep_cell_id <- WhichCells(object = temp, expression = seurat_clusters %in% c(0,1,2,3,6) & DF.classifications_0.25_0.19_1115 == "Singlet" & QC == "PASS")
p <- DimPlot(temp, reduction = "umap", cells.highlight = keep_cell_id) + 
  scale_color_manual(values = c("grey", "#F8766D")) + NoLegend()
select_cells <- CellSelector(p)
keep_cell_id <- setdiff(keep_cell_id, select_cells)

CTRL <- subset(CTRL, cells = keep_cell_id)
VlnPlot(object = CTRL, layer = "counts", features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3, pt.size = 0.1)

saveRDS(object = CTRL, file = "In-vitro/CTRL.rds")

#### DT ####
DT <- Read10X_h5(filename = "In-vitro/DT/filtered_feature_bc_matrix.h5")
DT <- CreateSeuratObject(counts = DT)
DT <- PercentageFeatureSet(DT, pattern = "^MT-", col.name = "percent_mito")

DT$Experiment <- "In-vitro"
DT$Condition <- "DT"

DT_bc <- fread(input = "In-vitro/DT/DT_bc.txt", header = T) %>% column_to_rownames(var = "CB")
DT <- AddMetaData(object = DT, metadata = DT_bc, col.name = "Barcode")

VlnPlot(object = DT, layer = "counts", features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3, pt.size = 0.1)

nCount_cutoff <- 50000
nFeature_cutoff <- 2000
pMito_cutoff <- 20

DT$QC <- "PASS"
DT$QC[which(DT$nCount_RNA > nCount_cutoff)] <- "High_nCount"
DT$QC[which(DT$nFeature_RNA < nFeature_cutoff)] <- "Low_nFeature"
DT$QC[which(DT$percent_mito > pMito_cutoff)] <- "High_pMito"
DT$QC[which(DT$nCount_RNA > nCount_cutoff & DT$nFeature_RNA < nFeature_cutoff)] <- "High_nCount; Low_nFeature"
DT$QC[which(DT$nCount_RNA > nCount_cutoff & DT$percent_mito > pMito_cutoff)] <- "High_nCount; High_pMito"
DT$QC[which(DT$nFeature_RNA < nFeature_cutoff & DT$percent_mito > pMito_cutoff)] <- "Low_nFeature; High_pMito"
DT$QC[which(DT$nCount_RNA > nCount_cutoff & DT$nFeature_RNA < nFeature_cutoff & DT$percent_mito > pMito_cutoff)] <- "High_nCount; Low_nFeature; High_pMito"

temp <- DT
temp <- NormalizeData(temp, verbose = FALSE)
temp <- FindVariableFeatures(temp, verbose = FALSE)
temp <- ScaleData(temp, verbose = FALSE)
temp <- RunPCA(temp, npcs = 20, verbose = FALSE)
temp <- RunUMAP(temp, dims = 1:10, verbose = FALSE)
DimPlot(object = temp, group.by = "QC")

temp <- FindNeighbors(temp)
temp <- FindClusters(temp, resolution = 0.3)
DimPlot(object = temp, group.by = "seurat_clusters", label = T, label.size = 5) + NoLegend()

sweep.res <- paramSweep(temp, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
best_pK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
nExp_poi <- round(0.075 * nrow(temp@meta.data))
temp <- doubletFinder(temp, PCs = 1:10, pN = 0.25, pK = best_pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(object = temp, group.by = "DF.classifications_0.25_0.16_1431")

keep_cell_id <- WhichCells(object = temp, expression = seurat_clusters %in% c(0,1,2,3,5) & DF.classifications_0.25_0.16_1431 == "Singlet" & QC == "PASS")
p <- DimPlot(temp, reduction = "umap", cells.highlight = keep_cell_id) + 
  scale_color_manual(values = c("grey", "#F8766D")) + NoLegend()
select_cells <- CellSelector(p)
keep_cell_id <- setdiff(keep_cell_id, select_cells)

DT <- subset(DT, cells = keep_cell_id)
VlnPlot(object = DT, layer = "counts", features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3, pt.size = 0.1)

saveRDS(object = DT, file = "In-vitro/DT.rds")

#### D0 ####
D0 <- Read10X_h5(filename = "In-vivo/D0/filtered_feature_bc_matrix.h5")
D0_mice <- fread(input = "In-vivo/D0/D0_mouse_cells.txt", header = F)
D0 <- D0[,colnames(D0) %ni% D0_mice$V1]

D0 <- CreateSeuratObject(counts = D0)
D0 <- PercentageFeatureSet(D0, pattern = "^MT-", col.name = "percent_mito")

D0$Experiment <- "In-vivo"
D0$Condition <- "D0"

D0_bc <- fread(input = "In-vivo/D0/D0_bc.txt", header = T) %>% column_to_rownames(var = "CB")
D0 <- AddMetaData(object = D0, metadata = D0_bc, col.name = "Barcode")

VlnPlot(object = D0, layer = "counts", features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3, pt.size = 0.1)

nCount_cutoff <- 75000
nFeature_cutoff <- 1000
pMito_cutoff <- 20

temp <- D0
temp$QC <- "PASS"
temp$QC[which(temp$nCount_RNA > nCount_cutoff)] <- "High_nCount"
temp$QC[which(temp$nFeature_RNA < nFeature_cutoff)] <- "Low_nFeature"
temp$QC[which(temp$percent_mito > pMito_cutoff)] <- "High_pMito"
temp$QC[which(temp$nCount_RNA > nCount_cutoff & temp$nFeature_RNA < nFeature_cutoff)] <- "High_nCount; Low_nFeature"
temp$QC[which(temp$nCount_RNA > nCount_cutoff & temp$percent_mito > pMito_cutoff)] <- "High_nCount; High_pMito"
temp$QC[which(temp$nFeature_RNA < nFeature_cutoff & temp$percent_mito > pMito_cutoff)] <- "Low_nFeature; High_pMito"
temp$QC[which(temp$nCount_RNA > nCount_cutoff & temp$nFeature_RNA < nFeature_cutoff & temp$percent_mito > pMito_cutoff)] <- "High_nCount; Low_nFeature; High_pMito"

temp <- NormalizeData(temp, verbose = FALSE)
temp <- FindVariableFeatures(temp, verbose = FALSE)
temp <- ScaleData(temp, verbose = FALSE)
temp <- RunPCA(temp, npcs = 10, verbose = FALSE)
temp <- RunUMAP(temp, dims = 1:10, verbose = FALSE)
DimPlot(object = temp, group.by = "QC")

sweep.res <- paramSweep(temp, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
best_pK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric[c(-1:-3,-30:-32)])]))
nExp_poi <- round(0.075 * nrow(temp@meta.data))
temp <- doubletFinder(temp, PCs = 1:10, pN = 0.25, pK = best_pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
temp$"DF_classifications" <- temp[[paste("DF.classifications_0.25", best_pK, nExp_poi, sep = "_")]]
DimPlot(object = temp, group.by = "DF_classifications")

temp <- FindNeighbors(temp, reduction = "pca", dims = 1:10)
temp <- FindClusters(temp, resolution = 0.4)
DimPlot(object = temp, group.by = "seurat_clusters", label = T, label.size = 7)

keep_cell_id <- WhichCells(object = temp, expression = seurat_clusters %in% c(0,1,3,4,5,6) & DF_classifications == "Singlet")
p <- DimPlot(temp, reduction = "umap", cells.highlight = keep_cell_id) + 
  scale_color_manual(values = c("grey", "#F8766D")) + NoLegend()
select_cells <- CellSelector(p)
keep_cell_id <- setdiff(keep_cell_id, select_cells)

D0 <- subset(D0, cells = keep_cell_id)
VlnPlot(object = D0, layer = "counts", features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3, pt.size = 0.1)
D0 <- subset(D0, percent_mito <= 20)

saveRDS(object = D0, file = "In-vivo/D0_filtered.rds")

#summary
ncol(D0)
mean(D0$nFeature_RNA)
mean(D0$nCount_RNA)
mean(D0$percent_mito)
length(unique(D0$Barcode[which(!is.na(D0$Barcode))]))

D0_div <- table(as.vector(D0$Barcode[which(!is.na(D0$Barcode))]))
diversity(D0_div, index = "shannon")
diversity(D0_div, index = "shannon")/log(length(D0_div))
diversity(D0_div, index = "invsimpson")
(length(D0_div)-1)/log(sum(D0_div))

#### D21 ####
D21 <- Read10X_h5(filename = "In-vivo/D21/filtered_feature_bc_matrix.h5")
D21_mice <- fread(input = "In-vivo/D21/D21_mouse_cells.txt", header = F)
D21 <- D21[,colnames(D21) %ni% D21_mice$V1]

D21 <- CreateSeuratObject(counts = D21)
D21 <- PercentageFeatureSet(D21, pattern = "^MT-", col.name = "percent_mito")

D21$Experiment <- "In-vivo"
D21$Condition <- "D21"

D21_bc <- fread(input = "In-vivo/D21/D21_bc.txt", header = T) %>% column_to_rownames(var = "CB")
D21 <- AddMetaData(object = D21, metadata = D21_bc, col.name = "Barcode")

VlnPlot(object = D21, layer = "counts", features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3, pt.size = 0.1)

nCount_cutoff <- 75000
nFeature_cutoff <- 1000
pMito_cutoff <- 20

temp <- D21
temp$QC <- "PASS"
temp$QC[which(temp$nCount_RNA > nCount_cutoff)] <- "High_nCount"
temp$QC[which(temp$nFeature_RNA < nFeature_cutoff)] <- "Low_nFeature"
temp$QC[which(temp$percent_mito > pMito_cutoff)] <- "High_pMito"
temp$QC[which(temp$nCount_RNA > nCount_cutoff & temp$nFeature_RNA < nFeature_cutoff)] <- "High_nCount; Low_nFeature"
temp$QC[which(temp$nCount_RNA > nCount_cutoff & temp$percent_mito > pMito_cutoff)] <- "High_nCount; High_pMito"
temp$QC[which(temp$nFeature_RNA < nFeature_cutoff & temp$percent_mito > pMito_cutoff)] <- "Low_nFeature; High_pMito"
temp$QC[which(temp$nCount_RNA > nCount_cutoff & temp$nFeature_RNA < nFeature_cutoff & temp$percent_mito > pMito_cutoff)] <- "High_nCount; Low_nFeature; High_pMito"

temp <- NormalizeData(temp, verbose = FALSE)
temp <- FindVariableFeatures(temp, verbose = FALSE)
temp <- ScaleData(temp, verbose = FALSE)
temp <- RunPCA(temp, npcs = 10, verbose = FALSE)
temp <- RunUMAP(temp, dims = 1:10, verbose = FALSE)
DimPlot(object = temp, group.by = "QC")

sweep.res <- paramSweep(temp, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
best_pK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric[c(-1:-3,-30:-32)])]))
nExp_poi <- round(0.075 * nrow(temp@meta.data))
temp <- doubletFinder(temp, PCs = 1:10, pN = 0.25, pK = best_pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
temp$"DF_classifications" <- temp[[paste("DF.classifications_0.25", best_pK, nExp_poi, sep = "_")]]
DimPlot(object = temp, group.by = "DF_classifications")

temp <- FindNeighbors(temp, reduction = "pca", dims = 1:10)
temp <- FindClusters(temp, resolution = 0.4)
DimPlot(object = temp, group.by = "seurat_clusters", label = T, label.size = 7)

keep_cell_id <- WhichCells(object = temp, expression = seurat_clusters %in% c(0,1,2,5,6,7) & DF_classifications == "Singlet")
p <- DimPlot(temp, reduction = "umap", cells.highlight = keep_cell_id) + 
  scale_color_manual(values = c("grey", "#F8766D")) + NoLegend()
select_cells <- CellSelector(p)
keep_cell_id <- setdiff(keep_cell_id, select_cells)

D21 <- subset(D21, cells = keep_cell_id)
VlnPlot(object = D21, layer = "counts", features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3, pt.size = 0.1)
D21 <- subset(D21, percent_mito <= 20)

saveRDS(object = D21, file = "In-vivo/D21_filtered.rds")

#summary
ncol(D21)
mean(D21$nFeature_RNA)
mean(D21$nCount_RNA)
mean(D21$percent_mito)
length(unique(D21$Barcode[which(!is.na(D21$Barcode))]))

D21_div <- table(as.vector(D21$Barcode[which(!is.na(D21$Barcode))]))
diversity(D21_div, index = "shannon")
diversity(D21_div, index = "shannon")/log(length(D21_div))
diversity(D21_div, index = "invsimpson")
(length(D21_div)-1)/log(sum(D21_div))

#### EP ####
EP <- Read10X_h5(filename = "In-vivo/EP/filtered_feature_bc_matrix.h5")
EP_mice <- fread(input = "In-vivo/EP/EP_mouse_cells.txt", header = F)
EP <- EP[,colnames(EP) %ni% EP_mice$V1]

EP <- CreateSeuratObject(counts = EP)
EP <- PercentageFeatureSet(EP, pattern = "^MT-", col.name = "percent_mito")

EP$Experiment <- "In-vivo"
EP$Condition <- "EP"

EP_bc <- fread(input = "In-vivo/EP/EP_bc.txt", header = T) %>% column_to_rownames(var = "CB")
EP <- AddMetaData(object = EP, metadata = EP_bc, col.name = "Barcode")

VlnPlot(object = EP, layer = "counts", features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3, pt.size = 0.1)

nCount_cutoff <- 75000
nFeature_cutoff <- 1000
pMito_cutoff <- 20

temp <- EP
temp$QC <- "PASS"
temp$QC[which(temp$nCount_RNA > nCount_cutoff)] <- "High_nCount"
temp$QC[which(temp$nFeature_RNA < nFeature_cutoff)] <- "Low_nFeature"
temp$QC[which(temp$percent_mito > pMito_cutoff)] <- "High_pMito"
temp$QC[which(temp$nCount_RNA > nCount_cutoff & temp$nFeature_RNA < nFeature_cutoff)] <- "High_nCount; Low_nFeature"
temp$QC[which(temp$nCount_RNA > nCount_cutoff & temp$percent_mito > pMito_cutoff)] <- "High_nCount; High_pMito"
temp$QC[which(temp$nFeature_RNA < nFeature_cutoff & temp$percent_mito > pMito_cutoff)] <- "Low_nFeature; High_pMito"
temp$QC[which(temp$nCount_RNA > nCount_cutoff & temp$nFeature_RNA < nFeature_cutoff & temp$percent_mito > pMito_cutoff)] <- "High_nCount; Low_nFeature; High_pMito"

temp <- NormalizeData(temp, verbose = FALSE)
temp <- FindVariableFeatures(temp, verbose = FALSE)
temp <- ScaleData(temp, verbose = FALSE)
temp <- RunPCA(temp, npcs = 10, verbose = FALSE)
temp <- RunUMAP(temp, dims = 1:10, verbose = FALSE)
DimPlot(object = temp, group.by = "QC")

sweep.res <- paramSweep(temp, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
best_pK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric[c(-1:-3,-30:-32)])]))
nExp_poi <- round(0.075 * nrow(temp@meta.data))
temp <- doubletFinder(temp, PCs = 1:10, pN = 0.25, pK = best_pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
temp$"DF_classifications" <- temp[[paste("DF.classifications_0.25", best_pK, nExp_poi, sep = "_")]]
DimPlot(object = temp, group.by = "DF_classifications")

temp <- FindNeighbors(temp, reduction = "pca", dims = 1:10)
temp <- FindClusters(temp, resolution = 0.4)
DimPlot(object = temp, group.by = "seurat_clusters", label = T, label.size = 7)

keep_cell_id <- WhichCells(object = temp, expression = seurat_clusters %in% c(1,2,3,4,5,6,7,9) & DF_classifications == "Singlet")
p <- DimPlot(temp, reduction = "umap", cells.highlight = keep_cell_id) + 
  scale_color_manual(values = c("grey", "#F8766D")) + NoLegend()
select_cells <- CellSelector(p)
keep_cell_id <- setdiff(keep_cell_id, select_cells)

EP <- subset(EP, cells = keep_cell_id)
VlnPlot(object = EP, layer = "counts", features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3, pt.size = 0.1)
EP <- subset(EP, percent_mito <= 20)

saveRDS(object = EP, file = "In-vivo/EP_filtered.rds")

#summary
ncol(EP)
mean(EP$nFeature_RNA)
mean(EP$nCount_RNA)
mean(EP$percent_mito)
length(unique(EP$Barcode[which(!is.na(EP$Barcode))]))

EP_div <- table(as.vector(EP$Barcode[which(!is.na(EP$Barcode))]))
diversity(EP_div, index = "shannon")
diversity(EP_div, index = "shannon")/log(length(EP_div))
diversity(EP_div, index = "invsimpson")
(length(EP_div)-1)/log(sum(EP_div))
