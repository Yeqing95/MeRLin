#### Marine signatures ####

setwd("~/Documents/Wistar/Haiyin/MeRLin/Data/Processed data/in-vivo/")
library(Seurat)
library(AUCell)
library(readxl)

data <- readRDS(file = "EP_Clonocluster.rds")
data <- subset(data, bc_group %in% c("Barcode group 1", "Barcode group 2", "Barcode group 3", "Barcode group 4", "Barcode group 5"))
counts <- GetAssayData(object = data, assay = "RNA", layer = "data")

Pathway <- read_excel(path = "../../Source data/Marine_signatures.xlsx", sheet = "Sheet2")
Stem <- toupper(as.character(na.omit(Pathway$`Stem-like`)))
NCSC <- toupper(as.character(na.omit(Pathway$`Neural Crest-like`)))
Mesenchymal <- toupper(as.character(na.omit(Pathway$`Mesenchymal-like`)))
RNA <- toupper(as.character(na.omit(Pathway$`RNA processing`)))
AP <- toupper(as.character(na.omit(Pathway$`Antigen Presentation`)))
Stress <- toupper(as.character(na.omit(Pathway$`Stress-like`)))
Melanocytic <- toupper(as.character(na.omit(Pathway$Melanocytic)))

data$`Stem-like` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = Stem)))
data$`Neural Crest-like` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = NCSC)))
data$`Mesenchymal-like` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = Mesenchymal)))
data$`RNA processing` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = RNA)))
data$`Antigen Presentation` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = AP)))
data$`Stress-like` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = Stress)))
data$Melanocytic <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = Melanocytic)))

FeaturePlot(data, reduction = "umap_wf6", features = "Stem-like", cols = c("lightgrey", "#EF3B2C"))
FeaturePlot(data, reduction = "umap_wf6", features = "Neural Crest-like", cols = c("lightgrey", "#EF3B2C"))
FeaturePlot(data, reduction = "umap_wf6", features = "Mesenchymal-like", cols = c("lightgrey", "#EF3B2C"))
FeaturePlot(data, reduction = "umap_wf6", features = "RNA processing", cols = c("lightgrey", "#EF3B2C"))
FeaturePlot(data, reduction = "umap_wf6", features = "Antigen Presentation", cols = c("lightgrey", "#EF3B2C"))
FeaturePlot(data, reduction = "umap_wf6", features = "Stress-like", cols = c("lightgrey", "#EF3B2C"))
FeaturePlot(data, reduction = "umap_wf6", features = "Melanocytic", cols = c("lightgrey", "#EF3B2C"))

# DEGs
markers <- FindAllMarkers(object = data, assay = "RNA", group.by = "bc_group", logfc.threshold = 0.8, min.pct = 0.3, only.pos = T)

group1_signature <- markers$gene[which(markers$cluster == "Barcode group 1" & markers$p_val_adj < 0.1)]
group2_signature <- markers$gene[which(markers$cluster == "Barcode group 2" & markers$p_val_adj < 0.1)]
group3_signature <- markers$gene[which(markers$cluster == "Barcode group 3" & markers$p_val_adj < 0.1)]
group4_signature <- markers$gene[which(markers$cluster == "Barcode group 4" & markers$p_val_adj < 0.1)]
group5_signature <- markers$gene[which(markers$cluster == "Barcode group 5" & markers$p_val_adj < 0.1)]

data$group1_signature <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = group1_signature)))
data$group2_signature <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = group2_signature)))
data$group3_signature <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = group3_signature)))
data$group4_signature <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = group4_signature)))
data$group5_signature <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = group5_signature)))

p2 <- FeaturePlot(data, reduction = "umap_wf6", features = "group1_signature", cols = c("lightgrey", "#EF3B2C"))
p3 <- FeaturePlot(data, reduction = "umap_wf6", features = "group2_signature", cols = c("lightgrey", "#EF3B2C"))
p4 <- FeaturePlot(data, reduction = "umap_wf6", features = "group3_signature", cols = c("lightgrey", "#EF3B2C"))
p5 <- FeaturePlot(data, reduction = "umap_wf6", features = "group4_signature", cols = c("lightgrey", "#EF3B2C"))
p6 <- FeaturePlot(data, reduction = "umap_wf6", features = "group5_signature", cols = c("lightgrey", "#EF3B2C"))

plot_grid(p1, p2, p3, p4, p5, p6, ncol = 2)

WriteXLS::WriteXLS(x = markers, ExcelFileName = "../../Source data/EP_DEGs_by_barcode_groups.xlsx", SheetNames = "res", row.names = F, col.names = T)

DE1 <- FindMarkers(object = data, assay = "RNA", group.by = "bc_group", ident.1 = c("Barcode group 1","Barcode group 2", "Barcode group 3", "Barcode group 4"), ident.2 = "Barcode group 5", logfc.threshold = 1, min.pct = 0.5, only.pos = T)
WriteXLS::WriteXLS(x = DE1, ExcelFileName = "../../Source data/EP_DEGs_persisters_vs_group5.xlsx", SheetNames = "res", row.names = T, col.names = T)

DE2 <- FindMarkers(object = data, assay = "RNA", group.by = "bc_group", ident.1 = "Barcode group 1", ident.2 = "Barcode group 5", logfc.threshold = 0.8, min.pct = 0.3, only.pos = T)
DE3 <- FindMarkers(object = data, assay = "RNA", group.by = "bc_group", ident.1 = "Barcode group 2", ident.2 = "Barcode group 5", logfc.threshold = 0.8, min.pct = 0.3, only.pos = T)
DE4 <- FindMarkers(object = data, assay = "RNA", group.by = "bc_group", ident.1 = "Barcode group 3", ident.2 = "Barcode group 5", logfc.threshold = 0.8, min.pct = 0.3, only.pos = T)
DE5 <- FindMarkers(object = data, assay = "RNA", group.by = "bc_group", ident.1 = "Barcode group 4", ident.2 = "Barcode group 5", logfc.threshold = 0.8, min.pct = 0.3, only.pos = T)

all <- unique(Reduce(f = union, x = list(rownames(DE2), rownames(DE3), rownames(DE4), rownames(DE5))))

WriteXLS::WriteXLS(x = DE2, ExcelFileName = "../../Source data/EP_DEGs_group1_vs_group5.xlsx", SheetNames = "res", row.names = T, col.names = T)
WriteXLS::WriteXLS(x = DE3, ExcelFileName = "../../Source data/EP_DEGs_group2_vs_group5.xlsx", SheetNames = "res", row.names = T, col.names = T)
WriteXLS::WriteXLS(x = DE4, ExcelFileName = "../../Source data/EP_DEGs_group3_vs_group5.xlsx", SheetNames = "res", row.names = T, col.names = T)
WriteXLS::WriteXLS(x = DE5, ExcelFileName = "../../Source data/EP_DEGs_group4_vs_group5.xlsx", SheetNames = "res", row.names = T, col.names = T)

DE6 <- FindMarkers(object = data, assay = "RNA", group.by = "bc_group", ident.1 = "Barcode group 3", ident.2 = "Barcode group 2", logfc.threshold = 0.5, min.pct = 0.3, only.pos = T)
WriteXLS::WriteXLS(x = DE6, ExcelFileName = "../../Source data/EP_DEGs_group3_vs_group2.xlsx", SheetNames = "res", row.names = T, col.names = T)

#### MeRLin signature ####
setwd("~/Documents/Wistar/Haiyin/MeRLin/Data/Processed data/in-vivo/")
library(Seurat)
library(AUCell)
library(readxl)

data <- readRDS(file = "EP_Clonocluster.rds")
data <- subset(data, bc_group %in% c("Barcode group 1", "Barcode group 2", "Barcode group 3", "Barcode group 4", "Barcode group 5"))
counts <- GetAssayData(object = data, assay = "RNA", layer = "data")

Signature <- read_excel(path = "../../Source data/MeRLin_signatures.xlsx", sheet = "Identified signatures")
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

tiff(filename = "../../../Figures/In-vivo/scRNAseq/ClonoCluster/EP/Stress-like signature in EP.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(data, reduction = "umap_wf6", features = "Stress-like signature", pt.size = 0.4, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/ClonoCluster/EP/Neural Crest-like signature in EP.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(data, reduction = "umap_wf6", features = "Neural crest-like signature", pt.size = 0.4, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/ClonoCluster/EP/Lipid metabolism in EP.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(data, reduction = "umap_wf6", features = "Lipid metabolism", pt.size = 0.4, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/ClonoCluster/EP/PI3K signaling in EP.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(data, reduction = "umap_wf6", features = "PI3K signaling", pt.size = 0.4, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/ClonoCluster/EP/ECM remodeling in EP.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(data, reduction = "umap_wf6", features = "ECM remodeling", pt.size = 0.4, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/ClonoCluster/EP/Melanocytic signature in EP.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(data, reduction = "umap_wf6", features = "Melanocytic markers", pt.size = 0.4, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/ClonoCluster/EP/JUN in EP.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(data, reduction = "umap_wf6", features = "JUN", pt.size = 0.4, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/ClonoCluster/EP/EDNRB in EP.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(data, reduction = "umap_wf6", features = "EDNRB", pt.size = 0.4, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/ClonoCluster/EP/FZD3 in EP.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(data, reduction = "umap_wf6", features = "FZD3", pt.size = 0.4, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/ClonoCluster/EP/TYRP1 in EP.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(data, reduction = "umap_wf6", features = "TYRP1", pt.size = 0.4, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/ClonoCluster/EP/DCT in EP.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(data, reduction = "umap_wf6", features = "DCT", pt.size = 0.4, cols = c("lightgrey", "#FFE680", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/ClonoCluster/EP/SERPINE2 in EP.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(data, reduction = "umap_wf6", features = "SERPINE2", pt.size = 0.4, min.cutoff = "q1", cols = c("lightgrey", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/ClonoCluster/EP/DUSP4 in EP.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(data, reduction = "umap_wf6", features = "DUSP4", pt.size = 0.4, min.cutoff = "q1", cols = c("lightgrey", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/ClonoCluster/EP/DUSP6 in EP.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(data, reduction = "umap_wf6", features = "DUSP6", pt.size = 0.4, min.cutoff = "q1", cols = c("lightgrey", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/ClonoCluster/EP/CSPG4 in EP.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(data, reduction = "umap_wf6", features = "CSPG4", pt.size = 0.4, min.cutoff = "q1", cols = c("lightgrey", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/ClonoCluster/EP/ITGA6 in EP.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(data, reduction = "umap_wf6", features = "ITGA6", pt.size = 0.4, min.cutoff = "q1", cols = c("lightgrey", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/ClonoCluster/EP/ITGA7 in EP.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(data, reduction = "umap_wf6", features = "ITGA7", pt.size = 0.4, min.cutoff = "q1", cols = c("lightgrey", "#FF6B6B"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/ClonoCluster/EP/FXYD3 in EP.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(data, reduction = "umap_wf6", features = "FXYD3", pt.size = 0.4, min.cutoff = "q1", cols = c("lightgrey", "#FF6B6B"))
dev.off()

library(msigdbr)
genesets <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "KEGG_LEGACY")
gene_sets <- split(x = genesets$gene_symbol, f = genesets$gs_name)
pathway <- gene_sets$KEGG
data$pathway <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = intersect(pathway, Group_3_signature))))
FeaturePlot(data, reduction = "umap_wf6", features = "pathway", cols = c("lightgrey", "#EF3B2C"))

################################################# D21 ####################################################
setwd("~/Documents/Wistar/Haiyin/MeRLin/Data/Processed data/in-vivo/")
library(Seurat)
library(readxl)
library(AUCell)

data <- readRDS(file = "D21_Clonocluster.rds")
DimPlot(object = data, group.by = "bc_group")
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

FeaturePlot(data, reduction = "umap_wf6", features = "Stress-like signature", pt.size = 0.4, cols = c("lightgrey", "#EF3B2C"))
FeaturePlot(data, reduction = "umap_wf6", features = "Neural crest-like signature", pt.size = 0.4, cols = c("lightgrey", "#EF3B2C"))
FeaturePlot(data, reduction = "umap_wf6", features = "Lipid metabolism", pt.size = 0.4, cols = c("lightgrey", "#EF3B2C"))
FeaturePlot(data, reduction = "umap_wf6", features = "PI3K signaling", pt.size = 0.4, cols = c("lightgrey", "#EF3B2C"))
FeaturePlot(data, reduction = "umap_wf6", features = "ECM remodeling", pt.size = 0.4, cols = c("lightgrey", "#EF3B2C"))
FeaturePlot(data, reduction = "umap_wf6", features = "Melanocytic signature", pt.size = 0.4, cols = c("lightgrey", "#EF3B2C"))

# DEGs
markers <- FindMarkers(object = data, assay = "RNA", group.by = "bc_group", ident.1 = c("Barcode group 1", "Barcode group 2"), logfc.threshold = 1, min.pct = 0.5, only.pos = T)

group1_signature <- markers$gene[which(markers$cluster == "Barcode group 1" & markers$p_val_adj < 0.05)]
data$group1_signature <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = group1_signature)))

FeaturePlot(data, reduction = "umap_wf6", features = "TYRP1", cols = c("lightgrey", "#EF3B2C"))
FeaturePlot(data, reduction = "umap_wf6", features = "group2_signature", cols = c("lightgrey", "#EF3B2C"))

################################################ EP Genes #################################################

setwd("~/Documents/Wistar/Haiyin/MeRLin/Data/Processed data/in-vivo/")
library(Seurat)
library(readxl)

data <- readRDS(file = "EP_Clonocluster.rds")
data <- subset(data, bc_group %in% c("Barcode group 1", "Barcode group 2", "Barcode group 3", "Barcode group 4", "Barcode group 5"))

Signature1 <- read_excel(path = "../../Source data/MeRLin_signatures.xlsx", sheet = "Identified signatures")
Stress <- as.character(na.omit(Signature1$Stress_like_signature))
NC <- as.character(na.omit(Signature1$NC_like_signature))
Lipid <- as.character(na.omit(Signature1$Lipid_metabolism))
PI3K <- as.character(na.omit(Signature1$PI3K_signaling))
ECM <- as.character(na.omit(Signature1$ECM_remodeling))

Signature2 <- read_excel(path = "../../Source data/MeRLin_signatures.xlsx", sheet = "Epi screen")
FDR <- as.character(na.omit(Signature2$`FDR<0.15`))
P <- as.character(na.omit(Signature2$`P<0.05`))

plot_every_gene <- function(object, gene_list){
  for(i in seq_along(gene_list)){
    gene <- gene_list[i]
    p <- FeaturePlot(object, reduction = "umap_wf6", pt.size = 0.4, min.cutoff = "q1", features = gene, cols = c("lightgrey", "#EF3B2C"))
    
    figure_name <- paste0("~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/scRNAseq/ClonoCluster/EP/Genes/", gene, "_in_EP.tiff")
    tiff(filename = figure_name, width = 700, height = 500, res = 150)
    print(p)
    dev.off()
  }
}

plot_every_gene(object = data, gene_list = Stress)
plot_every_gene(object = data, gene_list = NC)
plot_every_gene(object = data, gene_list = Lipid)
plot_every_gene(object = data, gene_list = PI3K)
plot_every_gene(object = data, gene_list = ECM)
plot_every_gene(object = data, gene_list = FDR)
plot_every_gene(object = data, gene_list = P)

Signature3 <- read_excel(path = "../../Source data/MeRLin_signatures.xlsx", sheet = "MR")
MR <- as.character(na.omit(Signature3$MR))
plot_every_gene(object = data, gene_list = MR)


View(table(data$Barcode, data$bc_group))
data$plot <- "Other"
data$plot[which(data$Barcode == "GTTGAACGACCACAA")] <- "GTTGAACGACCACAA" # Group1 top1
data$plot[which(data$Barcode == "GATGATCAACAAGGA")] <- "GATGATCAACAAGGA" # Group2 top1
data$plot[which(data$Barcode == "CTTGGAGGTCCACAA")] <- "CTTGGAGGTCCACAA" # Group3 top1
data$plot[which(data$Barcode == "GTACATCATGTAGCT")] <- "GTACATCATGTAGCT" # Group4 top1

data$plot <- factor(data$plot, levels = c("GATCTACCTCCAGAT", "GATGATCAACAAGGA", "CTTGGAGGTCCACAA", "GTACATCATGTAGCT", "Other"))

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/scRNAseq/ClonoCluster/EP/Top1 barcode in each group.tiff", width = 1800, height = 1000, res = 300)
DimPlot(data, reduction = "umap_wf6", group.by = "plot", pt.size = 0.4) + 
  scale_color_manual(values = c("GATCTACCTCCAGAT" = "#E64B35FF",
                                "GATGATCAACAAGGA" = "#3C5488FF", 
                                "CTTGGAGGTCCACAA" = "#00A087FF", 
                                "GTACATCATGTAGCT" = "#7E6148FF", 
                                "Other" = "grey"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/ClonoCluster/EP/Genes/MITF in EP.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(data, reduction = "umap_wf6", features = "MITF", pt.size = 0.4, min.cutoff = "q1", cols = c("lightgrey", "#EF3B2C"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/ClonoCluster/EP/Genes/MKI67 in EP.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(data, reduction = "umap_wf6", features = "MKI67", pt.size = 0.4, min.cutoff = "q1", cols = c("lightgrey", "#EF3B2C"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/ClonoCluster/EP/Genes/MCAM in EP.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(data, reduction = "umap_wf6", features = "MCAM", pt.size = 0.4, min.cutoff = "q1", cols = c("lightgrey", "#EF3B2C"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/ClonoCluster/EP/Genes/HEY1 in EP.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(data, reduction = "umap_wf6", features = "HEY1", pt.size = 0.4, min.cutoff = "q1", cols = c("lightgrey", "#EF3B2C"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/ClonoCluster/EP/Genes/SPP1 in EP.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(data, reduction = "umap_wf6", features = "SPP1", pt.size = 0.4, min.cutoff = "q1", cols = c("lightgrey", "#EF3B2C"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Signature/EP/JUN in EP2.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(data, reduction = "umap_wf6", features = "JUN", pt.size = 0.4, min.cutoff = "q1", cols = c("lightgrey", "#EF3B2C"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Signature/EP/EDNRB in EP2.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(data, reduction = "umap_wf6", features = "EDNRB", pt.size = 0.4, min.cutoff = "q1", cols = c("lightgrey", "#EF3B2C"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Signature/EP/DCT in EP2.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(data, reduction = "umap_wf6", features = "DCT", pt.size = 0.4, min.cutoff = "q1", cols = c("lightgrey", "#EF3B2C"))
dev.off()

tiff(filename = "../../../Figures/In-vivo/scRNAseq/Signature/EP/TYRP1 in EP2.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(data, reduction = "umap_wf6", features = "TYRP1", pt.size = 0.4, min.cutoff = "q1", cols = c("lightgrey", "#EF3B2C"))
dev.off()

################################################ EnrichR ################################################
library(ggplot2)
library(readxl)

g1 <- read_excel(path = "~/Documents/Wistar/Haiyin/MeRLin/Data/Source data/EnrichR.xlsx", sheet = "Group 1")
g2 <- read_excel(path = "~/Documents/Wistar/Haiyin/MeRLin/Data/Source data/EnrichR.xlsx", sheet = "Group 2")
g3 <- read_excel(path = "~/Documents/Wistar/Haiyin/MeRLin/Data/Source data/EnrichR.xlsx", sheet = "Group 3")
g4 <- read_excel(path = "~/Documents/Wistar/Haiyin/MeRLin/Data/Source data/EnrichR.xlsx", sheet = "Group 4")
#g5 <- read_excel(path = "~/Documents/Wistar/Haiyin/MeRLin/Data/Source data/EnrichR.xlsx", sheet = "Group 5")

g1$Group <- "Stress-like signature"
g2$Group <- "Lipid metabolism"
g3$Group <- "PI3K signaling"
g4$Group <- "ECM remodeling"

df <- rbind(g1, g2, g3, g4)
df$Group <- factor(x = df$Group, levels = c("Stress-like signature", "Lipid metabolism", "PI3K signaling", "ECM remodeling"))

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/scRNAseq/Signature/EnrichR.tiff", width = 1600, height = 2050, res = 300)
ggplot(data = df, mapping = aes(x = reorder(Name, -`Adjusted P-value`), y = -log10(`Adjusted P-value`), fill = Group)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Stress-like signature" = "#E64B35FF",
                               "Lipid metabolism" = "#3C5488FF", 
                               "PI3K signaling" = "#00A087FF", 
                               "ECM remodeling" = "#7E6148FF")) +
  facet_wrap(~ Group, scales = "free", ncol = 1) +
  xlab(label = "") +
  ylab(label = "") +
  theme_classic() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 10))
dev.off()                


