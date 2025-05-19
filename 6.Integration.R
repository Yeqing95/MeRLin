setwd("~/Documents/Wistar/Haiyin/MeRLin/Data/Processed data/in-vivo/")
library(Seurat)
library(readxl)
library(ggplot2)
`%ni%` <- Negate(`%in%`)
options(future.globals.maxSize = 24000 * 1024^2)

# load data
D0 <- readRDS(file = "D0_filtered.rds")
D21 <- readRDS(file = "D21_filtered.rds")
EP <- readRDS(file = "EP_filtered.rds")

D0$Condition <- "Day 0"
D21$Condition <- "Day 21"
EP$Condition <- "Endpoint"

# integration
Data_combined <-  merge(x = D0, y = c(D21, EP), add.cell.ids = c("D0", "D21", "EP"), project = "Data_combined")
rm(D0, D21, EP)

Data_combined <- JoinLayers(Data_combined, assay = "RNA")
Data_combined <- subset(Data_combined, cells = names(which(!is.na(Data_combined$Barcode))))

div <- read_excel(path = "../../Processed data/in-vivo/div_all_15bp.xlsx")
Data_combined@meta.data$bc_type <- div$bc_type[match(x = Data_combined$Barcode, table = div$barcode)]
Data_combined$bc_type[which(Data_combined$bc_type == "Multi_fate")] <- "Multi-fate"

g1 <- read.table("ClonoCluster/group1_bc.txt")
g2 <- read.table("ClonoCluster/group2_bc.txt")
g3 <- read.table("ClonoCluster/group3_bc.txt")
g4 <- read.table("ClonoCluster/group4_bc.txt")
g5 <- read.table("ClonoCluster/group5_bc.txt")

Data_combined@meta.data$bc_group <- "Other"
Data_combined@meta.data$bc_group[which(Data_combined@meta.data$Barcode %in% g1$V1)] <- "Barcode group 1"
Data_combined@meta.data$bc_group[which(Data_combined@meta.data$Barcode %in% g2$V1)] <- "Barcode group 2"
Data_combined@meta.data$bc_group[which(Data_combined@meta.data$Barcode %in% g3$V1)] <- "Barcode group 3"
Data_combined@meta.data$bc_group[which(Data_combined@meta.data$Barcode %in% g4$V1)] <- "Barcode group 4"
Data_combined@meta.data$bc_group[which(Data_combined@meta.data$Barcode %in% g5$V1)] <- "Barcode group 5"

VlnPlot(object = Data_combined, layer = "counts", features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), group.by = "Condition", ncol = 3, pt.size = 0.1)

# normalization
Data_combined <- SCTransform(object = Data_combined, vst.flavor = "v2")
Data_combined <- RunPCA(Data_combined, npcs = 30)
ElbowPlot(Data_combined, ndims = 30)
Data_combined <- RunUMAP(Data_combined, reduction = "pca", n.components = 2, dims = 1:20)

# clustering
Data_combined <- FindNeighbors(Data_combined)
Data_combined <- FindClusters(Data_combined, resolution = 0.1)
DimPlot(Data_combined, reduction = "umap", group.by = "seurat_clusters")

# cell cycle
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
Data_combined <- CellCycleScoring(Data_combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

df <- as.data.frame(table(Data_combined$Condition, Data_combined$Phase))
colnames(df) <- c("Condition", "Phase", "Frequency")
df$Condition <- factor(x = df$Condition, levels = c("Day 0","Day 21","Endpoint"))
df$Phase <- factor(x = df$Phase, levels = c("G1","S","G2M"))

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/scRNAseq/Integration/Cell cycle phase.tiff", width = 1200, height = 900, compression = "lzw", res = 300)
ggplot(data = df, mapping = aes(x = Condition, y = Frequency, fill = Phase)) +
  geom_bar(position="fill", stat="identity", width = 0.6) +
  scale_fill_manual(values = c("G1" = "#3F3E40", "S" = "#CCCDCE", "G2M" = "#7C7D80")) +
  xlab(label = "Sample") +
  ylab(label = "Proportion") +
  theme_classic()
dev.off()

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/scRNAseq/Integration/2D UMAP Condition.tiff", width = 1400, height = 1000, res = 300)
DimPlot(Data_combined, reduction = "umap", group.by = "Condition") +
  scale_color_manual(values = c("Day 0" = "#F8766D",
                                "Day 21" = "#00BA38", 
                                "Endpoint" = "#619CFF"))
dev.off()

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/scRNAseq/Integration/2D UMAP Barcode type.tiff", width = 1500, height = 1000, res = 300)
DimPlot(Data_combined, reduction = "umap", group.by = "bc_type") +
  scale_color_manual(values = c("Persister" = "#F8766D",
                                "Sensitive" = "#00BA38", 
                                "Multi-fate" = "#619CFF"))
dev.off()

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/scRNAseq/Integration/MKI67.tiff", width = 1400, height = 1000, res = 300)
FeaturePlot(Data_combined, features = "MKI67", pt.size = 0.4, min.cutoff = "q1", cols = c("lightgrey", "#EF3B2C"))
dev.off()

EP <- subset(Data_combined, Condition == "EP")
EP <- subset(EP, bc_group %in% c("Barcode group 1", "Barcode group 2", "Barcode group 3", "Barcode group 4", "Barcode group 5"))

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/scRNAseq/Integration/2D UMAP Barcode type in EP.tiff", width = 1500, height = 1000, res = 300)
DimPlot(EP, reduction = "umap", group.by = "bc_type") +
  scale_color_manual(values = c("Persister" = "#F8766D",
                                "Sensitive" = "#00BA38", 
                                "Multi-fate" = "#619CFF"))
dev.off()

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/scRNAseq/Integration/2D UMAP Barcode group in EP.tiff", width = 1550, height = 1000, res = 300)
DimPlot(EP, reduction = "umap", group.by = "bc_group") + 
  scale_color_manual(values = c("Barcode group 1" = "#E64B35FF",
                                "Barcode group 2" = "#3C5488FF", 
                                "Barcode group 3" = "#00A087FF", 
                                "Barcode group 4" = "#7E6148FF", 
                                "Barcode group 5" = "grey",
                                "Other" = "white"))
dev.off()

# escape rate
table(EP$seurat_clusters, EP$bc_group)
df <- data.frame(Group = factor(x = c("Barcode group 1", "Barcode group 2", "Barcode group 3", "Barcode group 4", "Barcode group 5")),
                 Rate = c(870/999, 431/643, 234/256, 264/378, 46/267))

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/scRNAseq/Integration/Escape rate.tiff", width = 1450, height = 1500, compression = "lzw", res = 300)
ggplot(data = df, mapping = aes(x = Group, y = Rate, fill = Group)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c("Barcode group 1" = "#E64B35FF",
                               "Barcode group 2" = "#3C5488FF", 
                               "Barcode group 3" = "#00A087FF", 
                               "Barcode group 4" = "#7E6148FF", 
                               "Barcode group 5" = "grey")) +
  geom_text(aes(label = sprintf("%.1f%%", Rate * 100)), vjust = -0.5, size = 4) +
  xlab(label = "") +
  ylab(label = "Escape rate") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

table(EP$seurat_clusters, EP$bc_type)
df <- data.frame(Group = factor(x = c("Persister", "Sensitive", "Multi-fate"), levels = c("Persister", "Sensitive", "Multi-fate")),
                 Rate = c(1756/2242, 11/74, 78/227))

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/scRNAseq/Integration/Escape rate 2.tiff", width = 1450, height = 1500, compression = "lzw", res = 300)
ggplot(data = df, mapping = aes(x = Group, y = Rate, fill = Group)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c("Persister" = "#F8766D",
                               "Sensitive" = "#00BA38", 
                               "Multi-fate" = "#619CFF")) +
  #geom_text(aes(label = sprintf("%.1f%%", Rate * 100)), vjust = -0.5, size = 4) +
  xlab(label = "") +
  ylab(label = "Escape rate") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()


library(plotly)
Data_combined <- RunUMAP(Data_combined, reduction = "pca", n.components = 3, dims = 1:20)
Data_combined@meta.data = cbind(Data_combined@meta.data, Embeddings(object = Data_combined, reduction = "umap"))
plot.data <- FetchData(object = Data_combined, vars = c("umap_1", "umap_2", "umap_3", "Condition"))
plot.data$label <- paste(rownames(plot.data))

plot_ly(data = plot.data, 
        x = ~umap_1, y = ~umap_2, z = ~umap_3, 
        color = ~Condition, size = I(5),
        colors = c("#F8766D", "#00BA38", "#619CFF"),
        type = "scatter3d", 
        mode = "markers", 
        marker = list(size = 3, width=2),
        text=~label,
        hoverinfo="text")

# APA
library(MAAPER)
options(future.globals.maxSize = 24000 * 1024^2)

pas_annotation = readRDS(file = "human.PAS.hg38.rds")
gtf = "/Volumes/herlynm/linux/ychen/refSeq/10x/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz"
bam_c1 = "/Volumes/herlynm/linux/ychen/Haiyin/barcode_v3/scRNAseq/D0/outs/possorted_genome_bam.bam"
#bam_c1 = "/Volumes/herlynm/linux/ychen/Haiyin/barcode_v3/scRNAseq/CTRL/outs/possorted_genome_bam.bam"
#bam_c2 = "/Volumes/herlynm/linux/ychen/Haiyin/barcode_v3/scRNAseq/D21_DT/outs/possorted_genome_bam.bam"
bam_c2 = "/Volumes/herlynm/linux/ychen/Haiyin/barcode_v3/scRNAseq/EP_DT/outs/possorted_genome_bam.bam"
#bam_c2 = "/Volumes/herlynm/linux/ychen/Haiyin/barcode_v3/scRNAseq/DT/outs/possorted_genome_bam.bam"

maaper(gtf, # full path of the GTF file
       pas_annotation, # PAS annotation
       output_dir = "./CTRL_vs_DT", # output directory
       bam_c1, bam_c2, # full path of the BAM files
       read_len = 280, # read length
       ncores = 1  # number of cores used for parallel computation 
)

setwd("~/Documents/Wistar/Haiyin/MeRLin/Data/Processed data/in-vivo/APA/")
library(ggplot2)

APA1 <- read.table("./D0_vs_D21/gene.txt", header =T)
APA2 <- read.table("./D0_vs_D57/gene.txt", header =T)
APA3 <- read.table("./D0_vs_EP/gene.txt", header =T)
APA4 <- read.table("./CTRL_vs_DT/gene.txt", header =T)

APA1$REDu.pval.adj <- p.adjust(APA1$REDu.pval, method = "BH")
APA1$REDi.pval.adj <- p.adjust(APA1$REDi.pval, method = "BH")
APA2$REDu.pval.adj <- p.adjust(APA2$REDu.pval, method = "BH")
APA2$REDi.pval.adj <- p.adjust(APA2$REDi.pval, method = "BH")
APA3$REDu.pval.adj <- p.adjust(APA3$REDu.pval, method = "BH")
APA3$REDi.pval.adj <- p.adjust(APA3$REDi.pval, method = "BH")
APA4$REDu.pval.adj <- p.adjust(APA4$REDu.pval, method = "BH")
APA4$REDi.pval.adj <- p.adjust(APA4$REDi.pval, method = "BH")

APA1$REDu.regu <- "No"
APA1$REDu.regu[which(APA1$REDu.pval.adj < 0.05 & APA1$REDu > log(1.2))] <- "Lengthened"
APA1$REDu.regu[which(APA1$REDu.pval.adj < 0.05 & APA1$REDu < -log(1.2))] <- "Shortened"
APA2$REDu.regu <- "No"
APA2$REDu.regu[which(APA2$REDu.pval.adj < 0.05 & APA2$REDu > log(1.2))] <- "Lengthened"
APA2$REDu.regu[which(APA2$REDu.pval.adj < 0.05 & APA2$REDu < -log(1.2))] <- "Shortened"
APA3$REDu.regu <- "No"
APA3$REDu.regu[which(APA3$REDu.pval.adj < 0.05 & APA3$REDu > log(1.2))] <- "Lengthened"
APA3$REDu.regu[which(APA3$REDu.pval.adj < 0.05 & APA3$REDu < -log(1.2))] <- "Shortened"
APA4$REDu.regu <- "No"
APA4$REDu.regu[which(APA4$REDu.pval.adj < 0.05 & APA4$REDu > log(1.2))] <- "Lengthened"
APA4$REDu.regu[which(APA4$REDu.pval.adj < 0.05 & APA4$REDu < -log(1.2))] <- "Shortened"

APA1$REDi.regu <- "No"
APA1$REDi.regu[which(APA1$REDi.pval.adj < 0.05 & APA1$REDi > log(1.2))] <- "Lengthened"
APA1$REDi.regu[which(APA1$REDi.pval.adj < 0.05 & APA1$REDi < -log(1.2))] <- "Shortened"
APA2$REDi.regu <- "No"
APA2$REDi.regu[which(APA2$REDi.pval.adj < 0.05 & APA2$REDi > log(1.2))] <- "Lengthened"
APA2$REDi.regu[which(APA2$REDi.pval.adj < 0.05 & APA2$REDi < -log(1.2))] <- "Shortened"
APA3$REDi.regu <- "No"
APA3$REDi.regu[which(APA3$REDi.pval.adj < 0.05 & APA3$REDi > log(1.2))] <- "Lengthened"
APA3$REDi.regu[which(APA3$REDi.pval.adj < 0.05 & APA3$REDi < -log(1.2))] <- "Shortened"
APA4$REDi.regu <- "No"
APA4$REDi.regu[which(APA3$REDi.pval.adj < 0.05 & APA3$REDi > log(1.2))] <- "Lengthened"
APA4$REDi.regu[which(APA3$REDi.pval.adj < 0.05 & APA3$REDi < -log(1.2))] <- "Shortened"

APA1$Comparation <- "D0_vs_D21"
APA2$Comparation <- "D0_vs_D57"
APA3$Comparation <- "D0_vs_EP"
APA4$Comparation <- "CTRL_vs_DT"

APA1$Experiment <- "In-vivo"
APA2$Experiment <- "In-vivo"
APA3$Experiment <- "In-vivo"
APA4$Experiment <- "In-vitro"

common_genes <- Reduce(intersect, list(APA1$gene, APA2$gene, APA3$gene, APA4$gene))

df <- rbind(subset(APA1, gene %in% common_genes), 
            subset(APA2, gene %in% common_genes), 
            subset(APA3, gene %in% common_genes),
            subset(APA4, gene %in% common_genes))

df.REDu <- as.data.frame(table(df$Comparation, df$REDu.regu))
colnames(df.REDu) <- c("Comparation", "Type", "Counts")
df.REDu <- subset(df.REDu, Type != "No")
df.REDu$Counts <- ifelse(df.REDu$Type == "Shortened", -df.REDu$Counts, df.REDu$Counts)
df.REDu$Comparation <- factor(x = df.REDu$Comparation, levels = c("D0_vs_EP", "D0_vs_D57", "D0_vs_D21", "CTRL_vs_DT"))
df.REDu$Experiment <- c("In-vitro", "In-vivo", "In-vivo", "In-vivo","In-vitro", "In-vivo", "In-vivo", "In-vivo")

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/scRNAseq/Integration/APA (utr).tiff", width = 1500, height = 900, compression = "lzw", res = 300)
ggplot(data = df.REDu, mapping = aes(x = Counts, y = Comparation, fill = Type)) +
  geom_bar(stat="identity", width = 0.7) +
  facet_grid(rows = vars(Experiment), scales = "free_y", space = "free_y") +
  theme_classic()
dev.off()

df.REDi <- as.data.frame(table(df$Comparation, df$REDi.regu))
colnames(df.REDi) <- c("Comparation", "Type", "Counts")
df.REDi <- subset(df.REDi, Type != "No")
df.REDi$Counts <- ifelse(df.REDi$Type == "Shortened", -df.REDi$Counts, df.REDi$Counts)
df.REDi$Comparation <- factor(x = df.REDi$Comparation, levels = c("D0_vs_EP", "D0_vs_D57", "D0_vs_D21", "CTRL_vs_DT"))
df.REDi$Experiment <- c("In-vitro", "In-vivo", "In-vivo", "In-vivo","In-vitro", "In-vivo", "In-vivo", "In-vivo")

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/scRNAseq/Integration/APA (intron).tiff", width = 1500, height = 900, compression = "lzw", res = 300)
ggplot(data = df.REDi, mapping = aes(x = Counts, y = Comparation, fill = Type)) +
  geom_bar(stat="identity", width = 0.7) +
  facet_grid(rows = vars(Experiment), scales = "free_y", space = "free_y") +
  theme_classic()
dev.off()






