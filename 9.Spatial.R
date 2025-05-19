library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(readxl)
library(AUCell)
library(cowplot)

#load data
setwd("~/Documents/Wistar/Haiyin/MeRLin/Data/Original data/In-vivo/Spatial/T4_S2/")

data <- Read10X(data.dir = "matrix")
data <- CreateSeuratObject(data, min.cells = 1, min.features = 1)
img <- Read10X_Image(image.dir = "spatial")
img <- img[Cells(data)]
DefaultAssay(img) <- DefaultAssay(data)
data[["slice"]] <- img

data <- subset(data, cells = rownames(GetTissueCoordinates(data)))

mouse_genes <- data[1:15410,]
human_genes <- data[15411:36952,]
mouse_counts <- Matrix::colSums(mouse_genes@assays$RNA$counts)
human_counts <- Matrix::colSums(human_genes@assays$RNA$counts)
mouse_ratio <- mouse_counts / (mouse_counts + human_counts)
data[["percent_mouse"]] <- mouse_ratio

data$species <- "Human"
data$species[which(data$percent_mouse > 0.2)] <- "Mouse"

SpatialFeaturePlot(data, features = "percent_mouse")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent_mouse"))

SpatialDimPlot(data, group.by = "species")

data <- subset(data, percent_mouse <= 0.2)
data <- data[15411:36952,]

data <- NormalizeData(data)
data <- FindVariableFeatures(data)
data <- ScaleData(data)
data <- RunPCA(data, assay = "RNA")
ElbowPlot(data)

data <- RunUMAP(data, reduction = "pca", dims = 1:10)
#data <- FindNeighbors(data, dims = 1:10)
#data <- FindClusters(data, resolution = 0.2)

FeaturePlot(data, features = "MITF")
SpatialFeaturePlot(data, features = "MITF")

counts <- GetAssayData(object = data, assay = "RNA", layer = "data")

Signature <- read_excel(path = "~/Documents/Wistar/Haiyin/MeRLin/Data/Source data/MeRLin_signatures.xlsx", sheet = "Identified signatures")
Stress <- as.character(na.omit(Signature$Stress_like_signature))
NC <- as.character(na.omit(Signature$NC_like_signature))
Lipid <- as.character(na.omit(Signature$Lipid_metabolism))
PI3K <- as.character(na.omit(Signature$PI3K_signaling))
ECM <- as.character(na.omit(Signature$ECM_remodeling))
ME <- as.character(na.omit(Signature$Melanocytic_markers))

data$`Stress-like signature` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = Stress)))
data$`Neural crest-like signature` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = NC)))
data$`Lipid metabolism` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = Lipid)))
data$`PI3K signaling` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = PI3K)))
data$`ECM remodeling` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = ECM)))
data$`Melanocytic markers` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = ME)))

SpatialFeaturePlot(data, features = "Stress-like signature", image.alpha = 0.5)
SpatialFeaturePlot(data, features = "Neural crest-like signature", image.alpha = 0.5)
SpatialFeaturePlot(data, features = "Lipid metabolism", image.alpha = 0.5)
SpatialFeaturePlot(data, features = "PI3K signaling", image.alpha = 0.5)
SpatialFeaturePlot(data, features = "ECM remodeling", image.alpha = 0.5)
SpatialFeaturePlot(data, features = "Melanocytic markers", image.alpha = 0.5)

SpatialFeaturePlot(data, features = c("ALDOA", "BNIP3", "FAM162A", "P4HA1", "PDK1", "PGAM1", "PGK1", "SLC2A1"), image.alpha = 0.5)

df <- data.frame(row.names = colnames(data),
                 `Stress-like signature` = scale(data$`Stress-like signature`),
                 `Lipid metabolism` = scale(data$`Lipid metabolism`),
                 `PI3K signaling` = scale(data$`PI3K signaling`),
                 `ECM remodeling` = scale(data$`ECM remodeling`),
                 `Melanocytic markers` = scale(data$`Melanocytic markers`),
                 check.rows = F, check.names = F)

table(colnames(df)[apply(df, 1, which.max)])
data$group <- colnames(df)[apply(df, 1, which.max)]

library(cowplot)

p1 <- SpatialDimPlot(data, group.by = "group", image.alpha = 0.5, cols = c("Stress-like signature" = "#E64B35FF",
                                                                           "Lipid metabolism" = "#3C5488FF", 
                                                                           "PI3K signaling" = "#00A087FF", 
                                                                           "ECM remodeling" = "#7E6148FF", 
                                                                           "Melanocytic markers" = "grey"))

p2 <- SpatialDimPlot(data, image.alpha = 0.1, cells.highlight = WhichCells(data, expression = group == "Stress-like signature"), cols.highlight = c("#E64B35FF", "transparent")) + NoLegend() + ggtitle(label = "Stress-like cells")
p3 <- SpatialDimPlot(data, image.alpha = 0.1, cells.highlight = WhichCells(data, expression = group == "Lipid metabolism"), cols.highlight = c("#3C5488FF","transparent")) + NoLegend() + ggtitle(label = "Lipid metabolism cells")
p4 <- SpatialDimPlot(data, image.alpha = 0.1, cells.highlight = WhichCells(data, expression = group == "PI3K signaling"), cols.highlight = c("#00A087FF","transparent")) + NoLegend() + ggtitle(label = "PI3K signaling cells")
p5 <- SpatialDimPlot(data, image.alpha = 0.1, cells.highlight = WhichCells(data, expression = group == "ECM remodeling"), cols.highlight = c("#7E6148FF","transparent")) + NoLegend() + ggtitle(label = "ECM remodeling cells")
p6 <- SpatialDimPlot(data, image.alpha = 0.1, cells.highlight = WhichCells(data, expression = group == "Melanocytic markers"), cols.highlight = c("grey","transparent")) + NoLegend() + ggtitle(label = "Melanocytic cells")

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/scRNAseq/Spatial/T4_S1.tiff", width = 1300, height = 1570, res = 200)
plot_grid(p1, p2, p3, p4, p5, p6, ncol = 2)
dev.off()

library(SingleR)
library(SummarizedExperiment)

ref <- readRDS("~/Documents/Wistar/Haiyin/MeRLin/Data/Processed data/in-vivo/EP_Clonocluster.rds")
ref <- subset(ref, bc_group %in% c("Barcode group 1", "Barcode group 2", "Barcode group 3", "Barcode group 4", "Barcode group 5"))
ref <- NormalizeData(ref)

#common_genes <- intersect(rownames(ref), rownames(data))
#ref <- ref[common_genes,]
#data <- data[common_genes,]

ref_expr <- GetAssayData(ref, layer = "data")
ref_labels <- ref$bc_group
ref_se <- SummarizedExperiment(assays = list(logcounts = ref_expr))

query_expr <- GetAssayData(data, layer = "data")
query_se <- SummarizedExperiment(assays = list(logcounts = query_expr))

pred <- SingleR(test = query_se, ref = ref_se, labels = ref_labels)
data$transfer_labels <- pred$labels

SpatialDimPlot(data, group.by = "transfer_labels")

# dist to boundary
library(spatstat.geom)

coords <- as.data.frame(GetTissueCoordinates(data))
coords$label <- data$group

range_x <- range(coords$x)
range_y <- range(coords$y)
win <- owin(xrange = range_x, yrange = range_y)

ppp_obj <- ppp(x = coords$x, y = coords$y, window = win)

boundary_dists <- bdist.points(ppp_obj)
coords$boundary_dist <- boundary_dists
coords$label <- factor(x = coords$label, levels = c("Stress-like signature",
                                                    "Lipid metabolism",
                                                    "PI3K signaling",
                                                    "ECM remodeling",
                                                    "Melanocytic markers"))

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/scRNAseq/Spatial/WM4237 T4_S1/Dist_to_boundary.tiff", width = 1380, height = 1200, res = 300)
ggplot(data = coords, mapping = aes(x = label, y = boundary_dist, fill = label)) +
  geom_boxplot(width = 0.6, outliers = F) +
  geom_point(size = 0.3, position = position_jitter(width = 0.1)) +
  scale_fill_manual(values = c("Stress-like signature" = "#E64B35FF",
                               "Lipid metabolism" = "#3C5488FF", 
                               "PI3K signaling" = "#00A087FF", 
                               "ECM remodeling" = "#7E6148FF",
                               "Melanocytic markers" = "gray")) + 
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(x = NULL, y = "Mean dist to boundary")
dev.off()

mean_dist_per_label <- aggregate(boundary_dist ~ label, data = coords, FUN = mean)
mean_dist_per_label$label <- factor(x = mean_dist_per_label$label, levels = c("Stress-like signature",
                                                                              "Lipid metabolism",
                                                                              "PI3K signaling",
                                                                              "ECM remodeling",
                                                                              "Melanocytic markers"))

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/scRNAseq/Spatial/WM4237 T4_S1/Mean_dist_to_boundary.tiff", width = 1380, height = 1200, res = 300)
ggplot(data = mean_dist_per_label, mapping = aes(x = label, y = boundary_dist, fill = label)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = c("Stress-like signature" = "#E64B35FF",
                               "Lipid metabolism" = "#3C5488FF", 
                               "PI3K signaling" = "#00A087FF", 
                               "ECM remodeling" = "#7E6148FF",
                               "Melanocytic markers" = "gray")) + 
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(x = NULL, y = "Mean dist to boundary")
dev.off()

# spatial dependency score
library(spdep)

group <- "ECM remodeling"
label_A_coords <- coords[coords$label == group,]
coords_mat <- as.matrix(label_A_coords[, c("x", "y")])
score_vec <- data$`ECM remodeling`[which(data$group == group)]
nb <- knn2nb(knearneigh(coords_mat, k=5))
lw <- nb2listw(nb, style="W")
morans_I <- moran.test(score_vec, lw)

df <- data.frame(Group = factor(x = c("Stress-like signature", "Lipid metabolism", "PI3K signaling", "ECM remodeling", "Melanocytic markers"),
                                levels = c("Stress-like signature", "Lipid metabolism", "PI3K signaling", "ECM remodeling", "Melanocytic markers")),
                 morans_I = c(0.197056590, 0.441485139, 0.378686592, 0.078266108, 0.2846246367))

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/scRNAseq/Spatial/WM4237 T4_S1/Morans_Index.tiff", width = 1380, height = 1200, res = 300)
ggplot(data = df, mapping = aes(x = Group, y = morans_I, fill = Group)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = c("Stress-like signature" = "#E64B35FF",
                               "Lipid metabolism" = "#3C5488FF", 
                               "PI3K signaling" = "#00A087FF", 
                               "ECM remodeling" = "#7E6148FF",
                               "Melanocytic markers" = "gray")) + 
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(x = NULL, y = "Morans Index")
dev.off()

# cellchat
library(CellChat)

spatial.locs = GetTissueCoordinates(data, scale = NULL, cols = c("imagerow", "imagecol"))[,c(1,2)]

spot.size = 55
scalefactors = jsonlite::fromJSON(txt = file.path("spatial", 'scalefactors_json.json'))
conversion.factor = spot.size/scalefactors$spot_diameter_fullres
spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)

d.spatial <- computeCellDistance(coordinates = spatial.locs, ratio = spatial.factors$ratio, tol = spatial.factors$tol)
min(d.spatial[d.spatial!=0])

cellchat <- createCellChat(object = data, 
                           meta = data@meta.data, 
                           group.by = "group",
                           datatype = "spatial", 
                           coordinates = spatial.locs, 
                           spatial.factors = spatial.factors)

cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat, variable.both = F)

cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
                              distance.use = TRUE, interaction.range = 250, scale.distance = 0.1,
                              contact.dependent = TRUE, contact.range = 100)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")
netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Blues")


for (i in seq_along(cellchat@netP$pathways)){
  pathways.show <- cellchat@netP$pathways[i]
  par(mfrow=c(1,1), xpd = TRUE)
  
  p <- netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle", title.space = 10, vertex.label.cex = 1.2)
  
  tiff(filename = paste0("~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/scRNAseq/CellChat/",pathways.show,".tiff"), width = 1900, height = 1230, res = 200)
  print(p)
  dev.off()
}

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/scRNAseq/CellChat/WM4237 T4_S2/TGFb1 L-R pairs (Lipid).tiff", width = 1280, height = 900, res = 200)
netVisual_bubble(cellchat, signaling = "TGFb", dot.size.min = 6, dot.size.max = 12, sources.use = "Lipid metabolism", remove.isolate = FALSE, font.size = 14, angle.x = 45) + theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12))
dev.off()

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/scRNAseq/CellChat/WM4237 T4_S2/TGFb1 L-R pairs (Stress-like).tiff", width = 1350, height = 1000, res = 200)
netVisual_bubble(cellchat, signaling = "TGFb", dot.size.min = 6, dot.size.max = 12, sources.use = "Stress-like signature", remove.isolate = FALSE, font.size = 14, angle.x = 45) + theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12))
dev.off()


library(NMF)
library(ggalluvial)

selectK(cellchat, pattern = "outgoing")

nPatterns = 3

par(mfrow=c(1,1), xpd = TRUE)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/scRNAseq/CellChat/outgoing_k3.tiff", width = 1165, height = 820, res = 100)
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()

cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/scRNAseq/CellChat/incoming_k3.tiff", width = 1165, height = 820, res = 100)
netAnalysis_river(cellchat, pattern = "incoming")
dev.off()

