setwd("~/Documents/Wistar/Haiyin/Barcode_v3/scRNAseq/")
library(Seurat)
library(stringr)
library(ggplot2)

Data <- Read10X_h5(filename = "D0/D0_filtered_feature_bc_matrix.h5")
Data <- CreateSeuratObject(counts = Data)

Data[["percent.mito"]] <- PercentageFeatureSet(Data, pattern = "^GRCh38-MT-")
Data[["percent.mice"]] <- PercentageFeatureSet(Data, pattern = "^GRCm39-")

write.table(x = Cells(subset(Data, percent.mice > 0.2)), file = "D0_mouse_cells.txt", quote = F, row.names = F, col.names = F)

ggplot(Data@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mito)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method = "lm") +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_bw()

VlnPlot(Data, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.mice"), ncol = 2, pt.size = 0.1)

Data <- subset(Data, nFeature_RNA > 1000 & percent.mito < 20)
VlnPlot(Data, layer = "counts", features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.mice"), ncol = 2, pt.size = 0.1)
