################################################# Step1 ####################################################

library(Seurat)

data <- readRDS(file = "EP_filtered.rds")
data <- subset(data, cells = names(which(!is.na(data$Barcode))))
write.csv(t(as.matrix(data@assays$RNA$counts)),file = "/Volumes/herlynm/linux/ychen/Haiyin/barcode_v3/scRNAseq/SECNIC/sce_exp_EP.csv")

################################################# Step2 ####################################################

import os, sys
os.getcwd()
os.listdir(os.getcwd()) 

import loompy as lp
import numpy as np
import scanpy as sc
x = sc.read_csv("sce_exp.csv"); 
row_attrs = {"Gene": np.array(x.var_names),};
col_attrs = {"CellID": np.array(x.obs_names)};
lp.create("sce.loom", x.X.transpose(), row_attrs, col_attrs)

################################################# Step3 ####################################################

pyscenic grn --num_workers 16 --sparse --output sce.adj.tsv --method grnboost2 sce.loom hs_hgnc_tfs.txt
pyscenic ctx --num_workers 16 sce.adj.tsv hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather --annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname sce.loom --output sce.regulon.csv --all_modules --mask_dropouts --mode "dask_multiprocessing" --min_genes 10
pyscenic aucell --num_workers 16 sce.loom sce.regulon.csv --output sce_SCENIC.loom

################################################# Step4 ####################################################

setwd("~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/scRNAseq/SECNIC/")
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(ggplot2)
library(pheatmap)

data <- readRDS("~/Documents/Wistar/Haiyin/MeRLin/Data/Processed data/in-vivo/EP_Clonocluster.rds")

sce_SCENIC <- open_loom("/Volumes/herlynm/linux/ychen/Haiyin/barcode_v3/scRNAseq/SECNIC/sce_EP_SCENIC.loom")
regulons_incidMat <- get_regulons(sce_SCENIC, column.attr.name="Regulons")

regulons <- regulonsToGeneLists(regulons_incidMat)

regulonAUC <- get_regulons_AUC(sce_SCENIC, column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(sce_SCENIC)

div <- readxl::read_excel(path = "../../../../Data/Processed data/in-vivo/div_all_15bp.xlsx")
data@meta.data$bc_type <- div$bc_type[match(x = data$Barcode, table = div$barcode)]
Type <- as.data.frame(subset(data@meta.data, select = "bc_type"))
selectedResolution <- "bc_type"

rss <- calcRSS(AUC = getAUC(regulonAUC), cellAnnotation = data@meta.data$bc_type)
rss <- na.omit(rss)
rssPlot <- plotRSS(rss, zThreshold = 0.8, cluster_columns = FALSE, order_rows = TRUE, thr = 0.01,
                   col.low = '#330066', col.mid = '#66CC66', col.high = '#FFCC33')
rssPlot

Group <- as.data.frame(subset(data@meta.data, select = "bc_group"))
selectedResolution <- "bc_group"

rss <- calcRSS(AUC = getAUC(regulonAUC), cellAnnotation = data@meta.data$bc_group)
rss <- na.omit(rss)
rss <- rss[,c(3,2,1,5,4)]
rss2 <- rss[c("ATF4(+)", "PRDX5(+)", "TAGLN2(+)",
              "ETV5(+)", "ETV4(+)", "FOSL1(+)",
              "LEF1(+)", "BNC1(+)", "SP6(+)",
              "ETS1(+)", "BRCA1(+)", "MTHFD1(+)",
              "JUN(+)", "FOSB(+)", "BPTF(+)"),]
rssPlot <- plotRSS(rss2, zThreshold = 0.8, cluster_columns = T, order_rows = T, thr = 0.1,
                   col.low = '#330066', col.mid = '#66CC66', col.high = '#FFCC33')

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/scRNAseq/SECNIC/EP/By barcode group.tiff", width = 960, height = 1410, res = 300)
rssPlot
dev.off()

regulon_AUC <- regulonAUC@NAMES
data@meta.data = cbind(data@meta.data, t(SummarizedExperiment::assay(regulonAUC[regulon_AUC,])))

for(i in colnames(data@meta.data)[10:65]){
  tiff(filename = paste0("UMAP_",i,".tif"), width = 572, height = 383, res = 100, compression = "lzw")
  p <- FeaturePlot(data, features = i, cols = c("lightgrey", "#FFD700", "#FF6B6B"))
  print(p)
  dev.off()
}

############################################## Addition ###################################################
library(Seurat)
library(readxl)

EP <- readRDS("~/Documents/Wistar/Haiyin/MeRLin/Data/Processed data/in-vivo/EP_Clonocluster.rds")
EP <- subset(EP, bc_group %in% c("Barcode group 1", "Barcode group 2", "Barcode group 3", "Barcode group 4", "Barcode group 5"))

FeaturePlot(EP, features = "ETS1")

Encode <- read_excel("~/Documents/Wistar/Haiyin/MeRLin/Data/Source data/ENCODE Transcription Factor Targets.xlsx")
targets <- intersect(Encode$Target[which(Encode$ETS1 == 1)], rownames(EP))

Signature <- read_excel(path = "~/Documents/Wistar/Haiyin/MeRLin/Data/Source data/MeRLin_signatures.xlsx", sheet = "Identified signatures")
ECM <- as.character(na.omit(Signature$ECM_remodeling))

Signature <- read_excel(path = "~/Documents/Wistar/Haiyin/MeRLin/Data/Source data/MeRLin_signatures.xlsx")
ECM <- as.character(na.omit(Signature$Group_4_signature))

G1 <- subset(EP, bc_group == "Barcode group 1")
G2 <- subset(EP, bc_group == "Barcode group 2")
G3 <- subset(EP, bc_group == "Barcode group 3")
G4 <- subset(EP, bc_group == "Barcode group 4")
G5 <- subset(EP, bc_group == "Barcode group 5")
NonG4 <- subset(EP, bc_group != "Barcode group 4")

mat <- GetAssayData(EP, assay = "RNA", layer = "scale.data")
ets1_expr <- mat["ETS1", ]
cors <- apply(mat, 1, function(x) cor(x, ets1_expr, method = "pearson"))
View(as.data.frame(cors[ECM]))


mat <- GetAssayData(NonG4, assay = "RNA", layer = "scale.data")
ets1_expr <- mat["ETS1", ]
cors <- apply(mat, 1, function(x) cor(x, ets1_expr, method = "pearson"))
View(as.data.frame(cors[ECM]))

