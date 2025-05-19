setwd("~/Documents/Wistar/Haiyin/MeRLin/Data/Processed data/in-vitro/")
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
options(future.globals.maxSize = 12000 * 1024^2)

# Load data
CTRL <- readRDS("CTRL_filtered.rds")
#CTRL <- subset(CTRL, cells = names(which(!is.na(CTRL$Barcode))))

DT <- readRDS("DT_filtered.rds")
#DT <- subset(DT, cells = names(which(!is.na(DT$Barcode))))

Data <-  merge(x = CTRL, y = DT, add.cell.ids = c("CTRL", "DT"))
rm(CTRL, DT)

Data$Condition[which(Data$Condition == "DT")] <- "BRAFi/MEKi"

# Normalization
Data <- JoinLayers(object = Data, assay = "RNA")
Data <- NormalizeData(Data, verbose = FALSE)
Data <- FindVariableFeatures(Data, verbose = FALSE)
Data <- ScaleData(Data, verbose = FALSE)
Data <- RunPCA(Data, npcs = 20, verbose = FALSE)
ElbowPlot(Data)

Data <- RunUMAP(Data, dims = 1:10, verbose = FALSE)
Data$Condition <- factor(x = Data$Condition, levels = c("BRAFi/MEKi", "CTRL"))

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

Data <- CellCycleScoring(Data, s.features = s.genes, g2m.features = g2m.genes)
DimPlot(Data, group.by = "Phase")

df <- as.data.frame(table(Data$Condition, Data$Phase))
colnames(df) <- c("Condition", "Phase", "Frequency")
df$Condition <- factor(x = df$Condition, levels = c("CTRL", "BRAFi/MEKi"))
df$Phase <- factor(x = df$Phase, levels = c("G1","S","G2M"))

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vitro/scRNAseq/Cell cycle phase.tiff", width = 1150, height = 960, compression = "lzw", res = 300)
ggplot(data = df, mapping = aes(x = Condition, y = Frequency, fill = Phase)) +
  geom_bar(position="fill", stat="identity", width = 0.6) +
  scale_fill_manual(values = c("G1" = "#3F3E40", "S" = "#CCCDCE", "G2M" = "#7C7D80")) +
  xlab(label = "") +
  ylab(label = "Proportion") +
  theme_classic()
dev.off()

Data <- FindNeighbors(Data)
Data <- FindClusters(Data, resolution = 0.3)
DimPlot(Data, group.by = "seurat_clusters")

Data$Cell_State <- "Non-cycling cells"
Data$Cell_State[which(Data$seurat_clusters %in% c(4,6))] <- "Cycling cells"
DimPlot(Data, group.by = "Cell_State")

#selected_cells <- CellSelector(DimPlot(Data, group.by = "Cell_State"))
#Data$Cell_State[selected_cells] <- "Non-cycling cells"
#Data$Cell_State[selected_cells] <- "Cycling cells"

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vitro/scRNAseq/UMAP_State.tiff", width = 1600, height = 1000, res = 300)
DimPlot(object = Data, group.by = "Cell_State")
dev.off()

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vitro/scRNAseq/UMAP_Condition.tiff", width = 1650, height = 1000, res = 300)
DimPlot(object = Data, group.by = "Condition")
dev.off()

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vitro/scRNAseq/Sensitive clone.tiff", width = 1400, height = 1000, res = 300)
highlight_cells <- WhichCells(Data, expression = Barcode == "GCAGATGGTGAACGT")
DimPlot(Data, reduction = "umap", cells.highlight = highlight_cells) + 
  scale_color_manual(values = c("grey", "#33C0C8")) + NoLegend()
dev.off()

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vitro/scRNAseq/Resistent clone.tiff", width = 1400, height = 1000, res = 300)
highlight_cells <- WhichCells(Data, expression = Barcode == "GTTCGTGGTGATGTT")
DimPlot(Data, reduction = "umap", cells.highlight = highlight_cells) + 
  scale_color_manual(values = c("grey", "#CE5F90")) + NoLegend()
dev.off()


tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vitro/scRNAseq/JUN.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(Data, features = "JUN", cols = c("lightgrey", "#FF6B6B"))
dev.off()

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vitro/scRNAseq/MKI67.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(Data, features = "MKI67", cols = c("lightgrey", "#FF6B6B"))
dev.off()

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vitro/scRNAseq/TYRP1.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(Data, features = "TYRP1", cols = c("lightgrey", "#FF6B6B"))
dev.off()

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vitro/scRNAseq/DCT.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(Data, features = "DCT", cols = c("lightgrey", "#FF6B6B"))
dev.off()

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vitro/scRNAseq/EDNRB.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(Data, features = "EDNRB", cols = c("lightgrey", "#FF6B6B"))
dev.off()

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vitro/scRNAseq/CCND1.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(Data, features = "CCND1", cols = c("lightgrey", "#FF6B6B"))
dev.off()

df <- Data@meta.data %>%
  select(Barcode, Condition) %>%
  filter(!is.na(Barcode)) %>%
  group_by(Barcode, Condition) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = Condition, values_from = n, values_fill = 0) %>%
  mutate(`BRAFi/MEKi` = `BRAFi/MEKi` / 5031, CTRL = CTRL / 7375) %>%
  mutate(Change = `BRAFi/MEKi` - CTRL)

df_long <- df %>%
  pivot_longer(cols = c("CTRL", "BRAFi/MEKi"), names_to = "Condition", values_to = "Percentage")

df_long <- subset(df_long, !Barcode %in% c("CATCCAGTACAACTA","CTAGTAGCTGGAGCC"))
df_long$Condition <- factor(x = df_long$Condition, levels = c("CTRL", "BRAFi/MEKi"))

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vitro/scRNAseq/Barcode Abundance Change (Percentage).tiff", width = 1300, height = 1200, res = 300)
ggplot(df_long, aes(x = Condition, y = Percentage, group = Barcode)) +
  geom_line(aes(color = Change), linewidth = 1) +
  geom_point(size = 1.5) +
  scale_y_continuous(labels = scales::label_percent()) +
  scale_color_gradient2(low = "blue", mid = "gray", high = "red", midpoint = 0) +
  theme_classic(base_size = 12) +
  labs(title = "Barcode Abundance Change",
       x = "Condition",
       y = "Cell Percentage",
       color = "Change")
dev.off()

# addition 2025-04-29
hist(df_long$Change,40)
sensitive_bc <- unique(df_long$Barcode[which(df_long$Change <= -0.001)])
resistant_bc <- unique(df_long$Barcode[which(df_long$Change >= 0.001)])

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vitro/scRNAseq/Sensitive clones.tiff", width = 1400, height = 1000, res = 300)
highlight_cells <- WhichCells(Data, expression = Barcode %in% sensitive_bc)
DimPlot(Data, reduction = "umap", cells.highlight = highlight_cells, pt.size = 0.1) + 
  scale_color_manual(values = c("grey", "#33C0C8")) + NoLegend()
dev.off()

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vitro/scRNAseq/Resistent clones.tiff", width = 1400, height = 1000, res = 300)
highlight_cells <- WhichCells(Data, expression = Barcode %in% resistant_bc)
DimPlot(Data, reduction = "umap", cells.highlight = highlight_cells, pt.size = 0.1) + 
  scale_color_manual(values = c("grey", "#CE5F90")) + NoLegend()
dev.off()

Data$bc_type <- NA
Data$bc_type[which(Data$Barcode %in% sensitive_bc)] <- "Sensitive clones"
Data$bc_type[which(Data$Barcode %in% resistant_bc)] <- "Resistant clones"

Data_sub <- subset(Data, bc_type %in% c("Sensitive clones", "Resistant clones"))

df <- as.data.frame(table(Data_sub$Condition, Data_sub$Cell_State, Data_sub$bc_type))
colnames(df) <- c("Condition", "State", "Type", "Frequency")
df$Condition <- factor(x = df$Condition, levels = c("CTRL", "BRAFi/MEKi"))
df$State <- factor(x = df$State, levels = c("Cycling cells", "Non-cycling cells"))

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vitro/scRNAseq/Cell state.tiff", width = 1350, height = 960, compression = "lzw", res = 300)
ggplot(data = df, mapping = aes(x = Type, y = Frequency, fill = State)) +
  geom_bar(position="fill", stat="identity", width = 0.7) +
  scale_fill_manual(values = c("Non-cycling cells" = "#3F3E40", "Cycling cells" = "#CCCDCE")) +
  scale_y_continuous(transform = "log10", breaks = c(0,1,2,3,4,5,6,7,8,9,10)) +
  facet_grid(~Condition) +
  xlab(label = "") +
  ylab(label = "log10(Proportion)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
dev.off()

# AUCell
library(AUCell)
library(msigdbr)

Data_sub$group <- NA
df <- as.data.frame(table(subset(Data_sub, Condition == "BRAFi/MEKi")$Barcode))

Data_sub$group[which(Data_sub$bc_type == "Resistant clones" & Data_sub$Barcode %in% df$Var1[which(df$Freq <= 15)])] <- "Resistant Small-sized clones"
Data_sub$group[which(Data_sub$bc_type == "Resistant clones" & Data_sub$Barcode %in% df$Var1[which(df$Freq > 15)])] <- "Resistant Large-sized clones"
Data_sub$group[which(Data_sub$bc_type == "Sensitive clones")] <- "Sensitive clones"

DimPlot(Data_sub, group.by = "group")

counts <- GetAssayData(object = Data_sub, assay = "RNA", layer = "data")

hallmark <- msigdbr(species = "Homo sapiens", collection = "H")

gene_sets <- split(x = hallmark$gene_symbol, f = hallmark$gs_name)
OP <- gene_sets$HALLMARK_OXIDATIVE_PHOSPHORYLATION
mT <- gene_sets$HALLMARK_MTORC1_SIGNALING
Uf <- gene_sets$HALLMARK_UNFOLDED_PROTEIN_RESPONSE

emt <- gene_sets$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
Gly <- gene_sets$HALLMARK_GLYCOLYSIS
CH <- gene_sets$HALLMARK_CHOLESTEROL_HOMEOSTASIS

TNF <- gene_sets$HALLMARK_TNFA_SIGNALING_VIA_NFKB
TGF <- gene_sets$HALLMARK_TGF_BETA_SIGNALING
Hypo <- gene_sets$HALLMARK_HYPOXIA
Hed <- gene_sets$HALLMARK_HEDGEHOG_SIGNALING

IAR <- gene_sets$HALLMARK_INTERFERON_ALPHA_RESPONSE

Data_sub$`Oxidative Phosphorylation` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = OP)))
Data_sub$`mTORC1 Signaling` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = mT)))
Data_sub$`Unfolded Protein Response` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = Uf)))

Data_sub$`Epithelial Mesenchymal Transition` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = emt)))
Data_sub$`Glycolysis` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = Gly)))
Data_sub$`Cholesterol Homeostasis` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = CH)))

Data_sub$`TNF-alpha Signaling via NF-kB` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = TNF)))
Data_sub$`TGF-beta Signaling` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = TGF)))
Data_sub$`Hypoxia` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = Hypo)))
Data_sub$`Hedgehog Signaling` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = Hed)))

Data_sub$`Interferon Alpha Response` <- as.numeric(getAUC(AUCell_run(exprMat = counts, geneSets = IAR)))

df <- FetchData(Data_sub, vars = c("Condition","group","Unfolded Protein Response","Oxidative Phosphorylation","TGF-beta Signaling","Hedgehog Signaling")) %>%
  pivot_longer(cols = -c(Condition, group), names_to = "Pathway", values_to = "Score")

df$x_axis <- interaction(df$Condition, df$group, sep = ": ")
df$x_axis <- factor(df$x_axis, levels = c(
  "CTRL: Resistant Large-sized clones",
  "CTRL: Resistant Small-sized clones",
  "CTRL: Sensitive clones",
  "BRAFi/MEKi: Resistant Large-sized clones",
  "BRAFi/MEKi: Resistant Small-sized clones",
  "BRAFi/MEKi: Sensitive clones"
))

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vitro/scRNAseq/jitter plot.tiff", width = 2600, height = 800, compression = "lzw", res = 200)
ggplot(df, aes(x = x_axis, y = Score, colour = group)) +
  geom_jitter(width = 0.2, size = 0.5) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5,
    color = "black", linewidth = 0.6, aes(ymin = ..y.., ymax = ..y..)) +
  facet_wrap(~Pathway, scales = "free", ncol = 4) +
  theme_classic() +
  xlab(label = "") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.placement = "outside", 
        strip.background = element_blank()) +
  theme(strip.text.x = element_text(size = 10),
        panel.spacing = unit(1, "lines"))
dev.off()

matrix <- as.matrix(GetAssayData(Data_sub, layer = "data"))
Data_sub <- AddMetaData(Data_sub, metadata = matrix["BCL2L1",], col.name = "BCL2L1")
Data_sub <- AddMetaData(Data_sub, metadata = matrix["MAP3K3",], col.name = "MAP3K3")
Data_sub <- AddMetaData(Data_sub, metadata = matrix["PIK3CB",], col.name = "PIK3CB")
Data_sub <- AddMetaData(Data_sub, metadata = matrix["STAT2",], col.name = "STAT2")

df <- FetchData(Data_sub, vars = c("Condition","group","BCL2L1","MAP3K3","PIK3CB","STAT2")) %>%
  pivot_longer(cols = -c(Condition, group), names_to = "Pathway", values_to = "Score")

df$x_axis <- interaction(df$Condition, df$group, sep = ": ")
df$x_axis <- factor(df$x_axis, levels = c(
  "CTRL: Resistant Large-sized clones",
  "CTRL: Resistant Small-sized clones",
  "CTRL: Sensitive clones",
  "BRAFi/MEKi: Resistant Large-sized clones",
  "BRAFi/MEKi: Resistant Small-sized clones",
  "BRAFi/MEKi: Sensitive clones"
))

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vitro/scRNAseq/boxplot2.tiff", width = 1700, height = 700, compression = "lzw", res = 200)
ggplot(subset(df, Score > 0), aes(x = x_axis, y = Score, colour = group)) +
  geom_boxplot() +
  facet_wrap(~Pathway, scales = "free", ncol = 6) +
  theme_classic() +
  ylab(label = "Expression") +
  xlab(label = "") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.placement = "outside", 
        strip.background = element_blank()) +
  theme(strip.text.x = element_text(size = 10),
        panel.spacing = unit(1, "lines"))
dev.off()



Data_sub <- subset(Data_sub, Condition == "DT")
res <- FindAllMarkers(Data_sub, group.by = "group", only.pos = T, logfc.threshold = 0.5, min.pct = 0.5)
View(subset(res, p_val_adj < 0.05))


# Ctrl:Sensitive vs Resistant
Data_nc <- subset(Data, Cell_State == "Non-cycling cells")
Data_nc_ctrl <- subset(Data_nc, Condition == "CTRL" & bc_type %in% c("Sensitive clones", "Resistant clones"))

de1 <- FindAllMarkers(Data_nc_ctrl, group.by = "bc_type", only.pos = T)
WriteXLS::WriteXLS(x = de1, ExcelFileName = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vitro/scRNAseq/DEGs/CTRL_Sensitive_vs_Resistant.xlsx")

de1_sig <- subset(de1, p_val_adj <= 0.05 & avg_log2FC > 0.8 & pct.1 > 0.4)
VlnPlot(Data_nc_ctrl, features = c("FSTL5"), group.by = "bc_type")

# DT:Sensitive vs Resistant
Data_nc_dt <- subset(Data_nc, Condition == "DT" & bc_type %in% c("Sensitive clones", "Resistant clones"))

de2 <- FindAllMarkers(Data_nc_dt, group.by = "bc_type", only.pos = T)
WriteXLS::WriteXLS(x = de2, ExcelFileName = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vitro/scRNAseq/DEGs/DT_Sensitive_vs_Resistant.xlsx")

de2_sig <- subset(de2, p_val_adj <= 0.05 & avg_log2FC > 0.8 & pct.1 > 0.4)
VlnPlot(Data_nc_dt, features = "ZNF627", group.by = "bc_type")

# Resistant:DT vs CTRL
Data_nc_r <- subset(Data_nc, bc_type == "Resistant clones")

de3 <- FindAllMarkers(Data_nc_r, group.by = "Condition", only.pos = T)
WriteXLS::WriteXLS(x = de3, ExcelFileName = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vitro/scRNAseq/DEGs/Resistant_DT_vs_CTRL.xlsx")

de3_sig <- subset(de3, p_val_adj <= 0.05 & avg_log2FC > 1 & pct.1 > 0.5)
VlnPlot(Data_nc_r, features = c("ASPA","ETV1"), group.by = "Condition")

# Sensitive:DT vs CTRL
Data_nc_s <- subset(Data_nc, bc_type == "Sensitive clones")

de4 <- FindAllMarkers(Data_nc_s, group.by = "Condition", only.pos = T)
WriteXLS::WriteXLS(x = de4, ExcelFileName = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vitro/scRNAseq/DEGs/Sensitive_DT_vs_CTRL.xlsx")

de4_sig <- subset(de4, p_val_adj <= 0.05 & avg_log2FC > 1 & pct.1 > 0.5)
VlnPlot(Data_nc_s, features = "ZNF627", group.by = "bc_type")

# All:DT vs CTRL
de <- FindAllMarkers(Data_nc, group.by = "Condition", only.pos = T)
WriteXLS::WriteXLS(x = de4, ExcelFileName = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vitro/scRNAseq/DEGs/All_DT_vs_CTRL.xlsx")

de_sig <- subset(de, p_val_adj <= 0.05 & avg_log2FC > 1 & pct.1 > 0.5)

library(VennDiagram)

set1 <- de_sig$gene[which(de_sig$cluster == "CTRL")]
set2 <- de3_sig$gene[which(de_sig$cluster == "CTRL")]
set3 <- de4_sig$gene[which(de_sig$cluster == "CTRL")]

library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

# Chart
venn.diagram(
  x = list(set1, set2, set3),
  category.names = c("All", "Resistant", "Sensitive"),
  filename = 'CTRL.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

# MR
library(readxl)
MR <- read_excel(path = "../../Source data/MeRLin_signatures.xlsx", sheet = "MR & MR_lit")
MR_lit <- as.character(na.omit(MR$MR_lit))
MR <- as.character(na.omit(MR$MR))

plot_every_gene <- function(object, gene_list){
  for(i in seq_along(gene_list)){
    gene <- gene_list[i]
    p <- FeaturePlot(object, reduction = "umap", pt.size = 0.4, min.cutoff = "q1", features = gene, cols = c("lightgrey", "#EF3B2C"))
    
    figure_name <- paste0("~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vitro/scRNAseq/Genes/MR_lit/", gene, ".tiff")
    tiff(filename = figure_name, width = 700, height = 500, res = 150)
    print(p)
    dev.off()
  }
}

plot_every_gene(object = Data, gene_list = MR_lit)



df <- data.frame(All = de1$avg_log2FC, RC = de2$avg_log2FC)
df$fc_diff <- df$RC - df$All
ggplot(data = df, mapping = aes(x = All, y = RC, color = fc_diff)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "gray", high = "red", midpoint = 0) +
  geom_abline(intercept = 1, slope = 1, linetype = "dashed") +
  geom_abline(intercept = -1, slope = 1, linetype = "dashed") +
  xlim(c(0,10)) +
  ylim(c(0,10)) +
  theme_bw()

