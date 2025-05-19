################################################ WM4007 ###################################################
library(Seurat)
library(readxl)
library(AUCell)

data <- readRDS("~/Documents/Wistar/Haiyin/WM4007/Data/WM4007/EP_DT_filtered.rds")
#data <- readRDS("~/Documents/Wistar/Haiyin/WM4007/Data/WM4007/D20_DT_filtered.rds")
#data <- readRDS("~/Documents/Wistar/Haiyin/WM4007/Data/WM4007/D13_DT_filtered.rds")

data <- NormalizeData(data, verbose = FALSE)
data <- FindVariableFeatures(data, verbose = FALSE)
data <- ScaleData(data, verbose = FALSE)
data <- RunPCA(data, npcs = 20, verbose = FALSE)
ElbowPlot(data)

data <- RunUMAP(data, dims = 1:15, verbose = FALSE)
counts <- GetAssayData(object = data, assay = "RNA", layer = "data")

Signature <- read_excel(path = "../../Source data/MeRLin_signatures.xlsx", sheet = "Sheet2")
Stress <- as.character(na.omit(Signature$Stress_like_signature))
NC <- as.character(na.omit(Signature$NC_like_signature))
Lipid <- as.character(na.omit(Signature$Lipid_metabolism))
PI3K <- as.character(na.omit(Signature$PI3K_signaling))
ECM <- as.character(na.omit(Signature$ECM_remodeling))
ME <- as.character(na.omit(Signature$Melanocytic_signature))

data$`Stress-like signature` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = Stress)))
data$`Neural crest-like signature` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = NC)))
data$`Lipid metabolism` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = Lipid)))
data$`PI3K signaling` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = PI3K)))
data$`ECM remodeling` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = ECM)))
data$`Melanocytic signature` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = ME)))

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Other models/WM4007 EP/Stress-like signature in WM4007 EP.tiff", width = 700, height = 500, res = 150)
FeaturePlot(data, reduction = "umap", features = "Stress-like signature", pt.size = 0.4, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Other models/WM4007 EP/Neural Crest-like signature in WM4007 EP.tiff", width = 700, height = 500, res = 150)
FeaturePlot(data, reduction = "umap", features = "Neural crest-like signature", pt.size = 0.4, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Other models/WM4007 EP/Lipid metabolism in WM4007 EP.tiff", width = 700, height = 500, res = 150)
FeaturePlot(data, reduction = "umap", features = "Lipid metabolism", pt.size = 0.4, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Other models/WM4007 EP/PI3K signaling in WM4007 EP.tiff", width = 700, height = 500, res = 150)
FeaturePlot(data, reduction = "umap", features = "PI3K signaling", pt.size = 0.4, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Other models/WM4007 EP/ECM remodeling in WM4007 EP.tiff", width = 700, height = 500, res = 150)
FeaturePlot(data, reduction = "umap", features = "ECM remodeling", pt.size = 0.4, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Other models/WM4007 EP/Melanocytic signature in WM4007 EP.tiff", width = 700, height = 500, res = 150)
FeaturePlot(data, reduction = "umap", features = "Melanocytic signature", pt.size = 0.4, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Other models/WM4007 EP/JUN in WM4007 EP.tiff", width = 700, height = 500, res = 150)
FeaturePlot(data, reduction = "umap", features = "JUN", pt.size = 0.4, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Other models/WM4007 EP/EDNRB in WM4007 EP.tiff", width = 700, height = 500, res = 150)
FeaturePlot(data, reduction = "umap", features = "EDNRB", pt.size = 0.4, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Other models/WM4007 EP/DCT in WM4007 EP.tiff", width = 700, height = 500, res = 150)
FeaturePlot(data, reduction = "umap", features = "DCT", pt.size = 0.4, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Other models/WM4007 EP/TYRP1 in WM4007 EP.tiff", width = 700, height = 500, res = 150)
FeaturePlot(data, reduction = "umap", features = "TYRP1", pt.size = 0.4, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Other models/WM4007 EP/FZD3 in WM4007 EP.tiff", width = 700, height = 500, res = 150)
FeaturePlot(data, reduction = "umap", features = "FZD3", pt.size = 0.4, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()


Signature <- read_excel(path = "~/Documents/Wistar/Haiyin/MeRLin/Data/Source data/MeRLin_signatures.xlsx", sheet = "Sheet1")
Group_1_signature <- as.character(na.omit(Signature$Group_1_signature))
Group_2_signature <- as.character(na.omit(Signature$Group_2_signature))
Group_3_signature <- as.character(na.omit(Signature$Group_3_signature))
Group_4_signature <- as.character(na.omit(Signature$Group_4_signature))
Group_5_signature <- as.character(na.omit(Signature$Group_5_signature))

data$Group_1_signature <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = Group_1_signature)))
data$Group_2_signature <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = Group_2_signature)))
data$Group_3_signature <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = Group_3_signature)))
data$Group_4_signature <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = Group_4_signature)))
data$Group_5_signature <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = Group_5_signature)))

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Other models/WM4007 EP/Group 1 signature in WM4007 EP.tiff", width = 700, height = 500, res = 150)
FeaturePlot(data, reduction = "umap", features = "Group_1_signature", cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Other models/WM4007 EP/Group 2 signature in WM4007 EP.tiff", width = 700, height = 500, res = 150)
FeaturePlot(data, reduction = "umap", features = "Group_2_signature", cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Other models/WM4007 EP/Group 3 signature in WM4007 EP.tiff", width = 700, height = 500, res = 150)
FeaturePlot(data, reduction = "umap", features = "Group_3_signature", cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Other models/WM4007 EP/Group 4 signature in WM4007 EP.tiff", width = 700, height = 500, res = 150)
FeaturePlot(data, reduction = "umap", features = "Group_4_signature", cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Other models/WM4007 EP/Group 5 signature in WM4007 EP.tiff", width = 700, height = 500, res = 150)
FeaturePlot(data, reduction = "umap", features = "Group_5_signature", cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

################################################ WM4380 ###################################################
library(Seurat)
library(readxl)
library(AUCell)

raw <- Read10X("/Volumes/herlynm/linux/ychen/scRNASeq/3rd/outs/filtered_feature_bc_matrix/")
data <- CreateSeuratObject(counts = raw$`Gene Expression`)
rownames(raw$`Antibody Capture`) <- c("Control","DT","Das","Das+DT")
data[['Hashtag']] <- CreateAssayObject(counts = raw$`Antibody Capture`)

data <- NormalizeData(data, assay = "Hashtag", normalization.method = "CLR")
data <- HTODemux(data, assay = "Hashtag", positive.quantile = 0.9)
table(data$Hashtag_classification.global)

Idents(data) <- "Hashtag_maxID"
RidgePlot(data, assay = "Hashtag", features = rownames(data[["Hashtag"]])[1:4], ncol = 2)

data <- subset(data, cells = WhichCells(data, expression = Hashtag_classification == "DT"))

data <- NormalizeData(data)
data <- FindVariableFeatures(data)
data <- ScaleData(data)
data <- RunPCA(data, assay = "RNA")
ElbowPlot(data)

data <- RunUMAP(data, reduction = "pca", dims = 1:15)

counts <- GetAssayData(object = data, assay = "RNA", layer = "data")

Signature <- read_excel(path = "../../Source data/MeRLin_signatures.xlsx", sheet = "Sheet2")
Stress <- as.character(na.omit(Signature$Stress_like_signature))
NC <- as.character(na.omit(Signature$NC_like_signature))
Lipid <- as.character(na.omit(Signature$Lipid_metabolism))
PI3K <- as.character(na.omit(Signature$PI3K_signaling))
ECM <- as.character(na.omit(Signature$ECM_remodeling))
ME <- as.character(na.omit(Signature$Melanocytic_signature))

data$`Stress-like signature` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = Stress)))
data$`Neural crest-like signature` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = NC)))
data$`Lipid metabolism` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = Lipid)))
data$`PI3K signaling` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = PI3K)))
data$`ECM remodeling` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = ECM)))
data$`Melanocytic signature` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = ME)))

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Other models/WM4380/Stress-like signature in WM4380 DT.tiff", width = 700, height = 500, res = 150)
FeaturePlot(data, reduction = "umap", features = "Stress-like signature", pt.size = 0.4, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Other models/WM4380/Neural Crest-like signature in WM4380 DT.tiff", width = 700, height = 500, res = 150)
FeaturePlot(data, reduction = "umap", features = "Neural crest-like signature", pt.size = 0.4, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Other models/WM4380/Lipid metabolism in WM4380 DT.tiff", width = 700, height = 500, res = 150)
FeaturePlot(data, reduction = "umap", features = "Lipid metabolism", pt.size = 0.4, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Other models/WM4380/PI3K signaling in WM4380 DT.tiff", width = 700, height = 500, res = 150)
FeaturePlot(data, reduction = "umap", features = "PI3K signaling", pt.size = 0.4, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Other models/WM4380/ECM remodeling in WM4380 DT.tiff", width = 700, height = 500, res = 150)
FeaturePlot(data, reduction = "umap", features = "ECM remodeling", pt.size = 0.4, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Other models/WM4380/Melanocytic signature in WM4380 DT.tiff", width = 700, height = 500, res = 150)
FeaturePlot(data, reduction = "umap", features = "Melanocytic signature", pt.size = 0.4, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Other models/WM4380/JUN in WM4380 DT.tiff", width = 700, height = 500, res = 150)
FeaturePlot(data, reduction = "umap", features = "JUN", pt.size = 0.4, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Other models/WM4380/EDNRB in WM4380 DT.tiff", width = 700, height = 500, res = 150)
FeaturePlot(data, reduction = "umap", features = "EDNRB", pt.size = 0.4, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Other models/WM4380/DCT in WM4380 DT.tiff", width = 700, height = 500, res = 150)
FeaturePlot(data, reduction = "umap", features = "DCT", pt.size = 0.4, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Other models/WM4380/TYRP1 in WM4380 DT.tiff", width = 700, height = 500, res = 150)
FeaturePlot(data, reduction = "umap", features = "TYRP1", pt.size = 0.4, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Other models/WM4380/FZD3 in WM4380 DT.tiff", width = 700, height = 500, res = 150)
FeaturePlot(data, reduction = "umap", features = "FZD3", pt.size = 0.4, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()


Signature <- read_excel(path = "~/Documents/Wistar/Haiyin/MeRLin/Data/Source data/MeRLin_signatures.xlsx", sheet = "Sheet1")
Group_1_signature <- as.character(na.omit(Signature$Group_1_signature))
Group_2_signature <- as.character(na.omit(Signature$Group_2_signature))
Group_3_signature <- as.character(na.omit(Signature$Group_3_signature))
Group_4_signature <- as.character(na.omit(Signature$Group_4_signature))
Group_5_signature <- as.character(na.omit(Signature$Group_5_signature))

data$Group_1_signature <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = Group_1_signature)))
data$Group_2_signature <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = Group_2_signature)))
data$Group_3_signature <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = Group_3_signature)))
data$Group_4_signature <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = Group_4_signature)))
data$Group_5_signature <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = Group_5_signature)))

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Other models/WM4380/Group 1 signature in WM4380 DT.tiff", width = 700, height = 500, res = 150)
FeaturePlot(data, reduction = "umap", features = "Group_1_signature", cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Other models/WM4380/Group 2 signature in WM4380 DT.tiff", width = 700, height = 500, res = 150)
FeaturePlot(data, reduction = "umap", features = "Group_2_signature", cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Other models/WM4380/Group 3 signature in WM4380 DT.tiff", width = 700, height = 500, res = 150)
FeaturePlot(data, reduction = "umap", features = "Group_3_signature", cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Other models/WM4380/Group 4 signature in WM4380 DT.tiff", width = 700, height = 500, res = 150)
FeaturePlot(data, reduction = "umap", features = "Group_4_signature", cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Other models/WM4380/Group 5 signature in WM4380 DT.tiff", width = 700, height = 500, res = 150)
FeaturePlot(data, reduction = "umap", features = "Group_5_signature", cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()


################################################ Grouping ##################################################
library(ggplot2)
library(ggpattern)
library(tidyr)
library(dplyr)
#library(pheatmap)

df <- data.frame(row.names = colnames(data),
                 `Stress-like signature` = scale(data$`Stress-like signature`),
                 `Lipid metabolism` = scale(data$`Lipid metabolism`),
                 `PI3K signaling` = scale(data$`PI3K signaling`),
                 `ECM remodeling` = scale(data$`ECM remodeling`),
                 `JUN positive` = data@assays$RNA$data["JUN",],
                 check.rows = F, check.names = F)

positive_values <- df$`JUN positive`[df$`JUN positive` > 0]
df$`JUN positive`[df$`JUN positive` == 0] <- -1
df$`JUN positive`[df$`JUN positive` > 0] <- scale(positive_values)

#softmax <- function(x) {
#  exp_x <- exp(x)
#  return(exp_x / sum(exp_x))
#}

#df2 <- t(apply(df, 1, softmax))
#colnames(df2) <- colnames(df)
#rownames(df2) <- rownames(df)

#pheatmap(t(df2), cluster_rows = F, cluster_cols = T, show_colnames = F)

table(colnames(df)[apply(df, 1, which.max)])

df <- data.frame(Model = c("WM4007 D13", "WM4007 D20", "WM4007 EP", "WM4237 EP", "WM4380 EP"),
                 `Stress-like signature` = c(1533,1762,2387,690,699),
                 `Lipid metabolism` = c(1918,2135,2584,399,687),
                 `PI3K signaling` = c(1963,2289,2641,430,685),
                 `ECM remodeling` = c(1722,2057,2333,485,765),
                 check.rows = F, check.names = F)

df_long <- pivot_longer(df, 
                        cols = -Model, 
                        names_to = "Signature", 
                        values_to = "Count")

df_long$Signature <- factor(df_long$Signature, levels = unique(names(df)[-1]))
df_long <- subset(df_long, Model != "WM4007 D20")

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/scRNAseq/Other models/Signature group proportion (color).tiff", width = 2000, height = 500, res = 200)
ggplot(df_long, aes(x = Model, y = Count, fill = Signature)) +
  geom_bar(position="dodge", stat="identity", width = 0.7) +
  scale_fill_manual(values = c("Stress-like signature" = "#E64B35FF",
                               "Lipid metabolism" = "#3C5488FF", 
                               "PI3K signaling" = "#00A087FF", 
                               "ECM remodeling" = "#7E6148FF")) + 
  facet_wrap(~Model, scales = "free", nrow = 1) +
  theme_classic(base_size = 12) +
  labs(x = NULL, y = "Cell count")
dev.off()

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/scRNAseq/Other models/Signature group proportion (texture).tiff", width = 2000, height = 500, res = 200)
ggplot(df_long, aes(x = Model, y = Count)) +
  geom_col_pattern(position = "dodge",
    aes(pattern = Signature, pattern_angle = Signature, pattern_spacing = Signature), 
    fill = 'white', colour = 'black', pattern_density = 0.5, pattern_fill = 'black', pattern_colour = 'darkgrey') +
  facet_wrap(~Model, scales = "free", nrow = 1) +
  theme_classic(base_size = 12) +
  labs(x = NULL, y = "Cell count")
dev.off()


