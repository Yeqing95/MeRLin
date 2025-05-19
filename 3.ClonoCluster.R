################################################# Step1 ####################################################

setwd("~/Documents/Wistar/Haiyin/MeRLin/Data/Processed data/in-vivo/")
library(Seurat)
library(dplyr)
`%ni%` <- Negate(`%in%`)

data <- readRDS(file = "EP_filtered.rds")
data <- subset(data, cells = names(which(!is.na(data$Barcode))))

data <- NormalizeData(object = data)
data <- ScaleData(object = data)

cm <- t(as.matrix(data@assays$RNA$scale.data)) %>% data.table::as.data.table(keep.rownames = TRUE)
data.table::fwrite(cm, "CloneCluster/EP_data.tsv", sep = "\t")

bt <- data.frame(rn = colnames(data), Barcode = data$Barcode)
data.table::fwrite(bt, "CloneCluster/EP_barcode.tsv", sep = "\t")

################################################# Step2 ####################################################

setwd("~/Documents/Wistar/Haiyin/MeRLin/Source data/Processed data/in-vivo/ClonoCluster/")
library(magrittr)
library(ClonoCluster)
library(ggplot2)

cm <- data.table::fread("D0_data.tsv")
cm %<>% dt2m
pca <- irlba_wrap(cm, npc = 10)

bt <- data.table::fread("D0_barcode.tsv")

clust <- clonocluster(pca, bt, alpha = seq(0, 1, by = 0.1), beta = 0.1, res = 1)

wfs <- seq(0, 10, by = 1)
umaps <- lapply(wfs, function(s){
  uws <- engage_warp(pca, bt, s)
  return(uws)
}) %>% data.table::rbindlist()

umaps <- merge(umaps, bt, by = "rn")
umaps[, Barcode :=
        ifelse(rn %>% unique %>% length > 1, Barcode, "Singlet"),
      by = "Barcode"]

tiff(filename = "../../../../Figures/In-vivo/scRNAseq/ClonoCluster/D0_diff_warps.tiff", width = 1200, height = 740, compression = "lzw", res = 150)
ggplot(umaps, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = 0.3, alpha = 0.5) +
  facet_wrap(~warp) +
  theme_void()
dev.off()

umap_wf6 <- engage_warp(pca, bt, s = 6) %>% select(c(rn, UMAP_1, UMAP_2)) %>% tibble::column_to_rownames(var = "rn")
data[["umap_wf6"]] <- CreateDimReducObject(embeddings = as.matrix(umap_wf6), key = "UMAP_", assay = DefaultAssay(data))

div <- readxl::read_excel(path = "../div_all_15bp.xlsx")
data@meta.data$bc_type <- div$bc_type[match(x = data$Barcode, table = div$barcode)]
data$bc_type[which(data$bc_type == "Multi_fate")] <- "Multi-fate"

tiff(filename = "../../../../Figures/In-vivo/scRNAseq/ClonoCluster/Sensitive barcodes in D0.tiff", width = 1400, height = 1000, compression = "lzw", res = 300)
highlight_cells <- WhichCells(data, expression = bc_type == "Sensitive")
DimPlot(data, reduction = "umap_wf6", cells.highlight = highlight_cells) + 
  scale_color_manual(values = c("grey", "#00BA38")) + NoLegend() + ggtitle("Sensitive barcodes")
dev.off()

tiff(filename = "../../../../Figures/In-vivo/scRNAseq/ClonoCluster/Persister barcodes in D0.tiff", width = 1400, height = 1000, compression = "lzw", res = 300)
highlight_cells <- WhichCells(data, expression = bc_type == "Persister")
DimPlot(data, reduction = "umap_wf6", cells.highlight = highlight_cells) + 
  scale_color_manual(values = c("grey", "#F8766D")) + NoLegend() + ggtitle("Persister barcodes")
dev.off()

tiff(filename = "../../../../Figures/In-vivo/scRNAseq/ClonoCluster/Multi-fate barcodes in D0.tiff", width = 1400, height = 1000, compression = "lzw", res = 300)
highlight_cells <- WhichCells(data, expression = bc_type == "Multi-fate")
DimPlot(data, reduction = "umap_wf6", cells.highlight = highlight_cells) + 
  scale_color_manual(values = c("grey", "#619CFF")) + NoLegend() + ggtitle("Multi-fate barcodes")
dev.off()

tiff(filename = "../../../../Figures/In-vivo/scRNAseq/ClonoCluster/Barcode types in D0.tiff", width = 1600, height = 1000, compression = "lzw", res = 300)
DimPlot(data, reduction = "umap_wf6", group.by = "bc_type") +
  scale_color_manual(values = c("Persister" = "#F8766D",
                                "Sensitive" = "#00BA38", 
                                "Multi-fate" = "#619CFF")) +
  ggtitle("Barcode types")
dev.off()

#data <- FindNeighbors(object = data, reduction = "umap_wf6", dims = 1:2)
#data <- FindClusters(object = data, resolution = 0.1)
#DimPlot(object = data, reduction = "umap_wf6", group.by = "seurat_clusters")

#df <- table(data$Barcode, data$seurat_clusters)
#max_clusters <- apply(df, 1, function(x) {
#  max_indices <- which(x == max(x))
#  if (length(max_indices) > 1) {
#    return(NA)
#  } else {
#    return(names(x)[max_indices])
#  }
#})

#highlight_cells <- WhichCells(data, expression = Barcode %in% names(which(max_clusters == "4")))
#DimPlot(data, reduction = "umap_wf6", cells.highlight = highlight_cells) + 
#  scale_color_manual(values = c("grey", "#F8766D")) + NoLegend()

#group1 <- names(which(max_clusters == "0")) #stress
#group2 <- names(which(max_clusters == "1")) #ncsc_bc1 (early)
#group3 <- names(which(max_clusters == "2")) #ncsc_bc2 (late)
#group4 <- names(which(max_clusters == "3")) #unknown
#group5 <- names(which(max_clusters == "4")) #dying

#write.table(x = group1, file = "group1_bc.txt", quote = F, sep = "\t", row.names = F, col.names = F)
#write.table(x = group2, file = "group2_bc.txt", quote = F, sep = "\t", row.names = F, col.names = F)
#write.table(x = group3, file = "group3_bc.txt", quote = F, sep = "\t", row.names = F, col.names = F)
#write.table(x = group4, file = "group4_bc.txt", quote = F, sep = "\t", row.names = F, col.names = F)
#write.table(x = group5, file = "group5_bc.txt", quote = F, sep = "\t", row.names = F, col.names = F)

group1_bc <- read.table(file = "group1_bc.txt")
group2_bc <- read.table(file = "group2_bc.txt")
group3_bc <- read.table(file = "group3_bc.txt")
group4_bc <- read.table(file = "group4_bc.txt")
group5_bc <- read.table(file = "group5_bc.txt")

data$bc_group <- "NA"
data$bc_group[which(data$Barcode %in% group1_bc$V1)] <- "Barcode group 1"
data$bc_group[which(data$Barcode %in% group2_bc$V1)] <- "Barcode group 2"
data$bc_group[which(data$Barcode %in% group3_bc$V1)] <- "Barcode group 3"
data$bc_group[which(data$Barcode %in% group4_bc$V1)] <- "Barcode group 4"
data$bc_group[which(data$Barcode %in% group5_bc$V1)] <- "Barcode group 5"

tiff(filename = "../../../../Figures/In-vivo/scRNAseq/ClonoCluster/Barcode groups in D0.tiff", width = 1700, height = 1000, compression = "lzw", res = 300)
DimPlot(data, reduction = "umap_wf6", group.by = "bc_group") + 
  scale_color_manual(values = c("Barcode group 1" = "#E64B35FF",
                                "Barcode group 2" = "#3C5488FF", 
                                "Barcode group 3" = "#00A087FF", 
                                "Barcode group 4" = "#7E6148FF", 
                                "Barcode group 5" = "grey")) + 
  ggtitle(label = "Barcode groups")
dev.off()

saveRDS(object = data, file = "../D0_Clonocluster.rds")

#highlight_cells <- WhichCells(data, expression = Barcode %in% stress_bc)
#tiff(filename = "../../../../Figures/In-vivo/scRNAseq/CloneCluster/Stress-like barcodes in D21.tiff", width = 1400, height = 1000, compression = "lzw", res = 300)
#DimPlot(data, reduction = "umap_wf6", cells.highlight = highlight_cells) + 
#  scale_color_manual(values = c("grey", "#E64B35FF")) + NoLegend()
#dev.off()

################################################# Step3 ####################################################
setwd("~/Documents/Wistar/Haiyin/MeRLin/Data/Processed data/in-vivo/")
library(Seurat)
library(ggplot2)

D0 <- readRDS("D0_Clonocluster.rds")
D21 <- readRDS("D21_Clonocluster.rds")
EP <- readRDS("EP_Clonocluster.rds")
EP <- subset(EP, bc_group %in% c("Barcode group 1", "Barcode group 2", "Barcode group 3", "Barcode group 4", "Barcode group 5"))

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/scRNAseq/ClonoCluster/D0/MKI67 in D0.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(D0, features = "MKI67", pt.size = 0.4, cols = c("lightgrey", "#EF3B2C"))
dev.off()

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/scRNAseq/ClonoCluster/D21/MKI67 in D21.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(D21, features = "MKI67", pt.size = 0.4, cols = c("lightgrey", "#EF3B2C"))
dev.off()

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/scRNAseq/ClonoCluster/EP/MKI67 in EP.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(EP, features = "MKI67", pt.size = 0.4, cols = c("lightgrey", "#EF3B2C"))
dev.off()

EP <- CellCycleScoring(EP, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
table(EP$bc_group, EP$Phase)

df <- as.data.frame(table(EP$bc_group, EP$Phase))
colnames(df) <- c("Group", "Phase", "Frequency")
df$Group <- factor(x = df$Group, levels = c("Barcode group 1","Barcode group 2","Barcode group 3","Barcode group 4","Barcode group 5"))
df$Phase <- factor(x = df$Phase, levels = c("G1","S","G2M"))

tiff(filename = "../../../Figures/In-vivo/scRNAseq/ClonoCluster/Cell cycle phase in EP (by bc_group).tiff", width = 1200, height = 800, compression = "lzw", res = 200)
ggplot(data = df, mapping = aes(x = Group, y = Frequency, fill = Phase)) +
  geom_bar(position="fill", stat="identity", width = 0.7) +
  scale_fill_manual(values = c("G1" = "#3F3E40", "S" = "#CCCDCE", "G2M" = "#7C7D80")) +
  xlab(label = "") +
  ylab(label = "Proportion") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

DimPlot(EP, group.by = "Phase")

D21 <- CellCycleScoring(D21, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
DimPlot(D21, group.by = "Phase")

D0 <- CellCycleScoring(D0, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
DimPlot(D0, group.by = "Phase")

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/scRNAseq/ClonoCluster/EP/Genes/SPP1 in EP.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(EP, features = "SPP1", pt.size = 0.4, cols = cols = c("lightgrey", "#EF3B2C"))
dev.off()

mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(96)

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/scRNAseq/Signature/EP/Heatmap.tiff", width = 1900, height = 2400, res = 300)
DoHeatmap(EP, features = c("ALDOA","BNIP3","FAM162A","P4HA1","PDK1","PGK1","SLC2A1",
                           "THAP5","FASN","APOC1","LXN","RDH5",
                           "ADCK2","MCAM", "HEY1", "SPP1","CCND1",
                           "FN1","BCL2","AKT3","FGFR1","TNC",
                           "MET","COL6A2","VCL","ECM1","HOOK1",
                           "EDNRB","DCT","JUN","MITF","TYRP1"),
          group.by = "bc_group", assay = "RNA", slot = "scale.data") +
  scale_fill_gradientn(colours = rev(mapal))
dev.off()

