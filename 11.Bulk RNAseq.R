######################################### Convert raw to tpm ###########################################
setwd("~/Documents/Wistar/Haiyin/MeRLin/Data/Original data/In-vivo/Bulk RNAseq/")
library(stringr)
library(dplyr)

list <- list.files(path = ".", pattern = "*tab", full.names = T)
raw_matrix <- matrix(data = NA, nrow = 28277, ncol = length(list))
colnames(raw_matrix) <- str_remove_all(string = list, pattern = "^./|.ReadsPerGene.out.tab")

for(i in 1:length(list)){
  temp <- read.table(file = list[i], header = F)
  raw_matrix[,i] <- temp$V2[5:28281]
}
rownames(raw_matrix) <- temp$V1[5:28281]

ribosomal_genes <- grep(pattern = "^RNA18S|^RNA28S|^RNA45S|^RNA5-8S|^RNA5S", x = rownames(raw_matrix), value = T)
#View(raw_matrix[which(rownames(raw_matrix) %in% ribosomal_genes),])
raw_matrix <- raw_matrix[which(!rownames(raw_matrix) %in% ribosomal_genes),]
raw_matrix <- raw_matrix[order(rownames(raw_matrix)),]
WriteXLS::WriteXLS(x = as.data.frame(raw_matrix), ExcelFileName = "../../../Processed data/in-vivo/Bulk RNAseq/raw_counts.xlsx", row.names = T, col.names = T)

hg38 <- read.table("~/Desktop/hg38.refGene.gtf", sep = "\t", header = F)
hg38 <- hg38[which(hg38$V3 == "transcript"),]
length <- hg38$V5 - hg38$V4 + 1
gene <- str_remove_all(string = str_split_fixed(string = hg38$V9, pattern = ";", n = 2)[,1], pattern = "gene_id ")
df <- data.frame(gene = gene, length = length)
df <- df %>% group_by(gene) %>% summarize(length = max(length))
df <- df[which(!df$gene %in% ribosomal_genes),]
df <- df[order(df$gene),]
df <- df[which(df$gene != "BAGE5"),]

tpm <- function(counts){
  tpm = 1000000*(counts/df$length)/sum(counts/df$length)
  return(tpm)
}

tpm_matrix <- apply(raw_matrix, 2, tpm)
WriteXLS::WriteXLS(x = as.data.frame(tpm_matrix), ExcelFileName = "../../../Processed data/in-vivo/Bulk RNAseq/tpm_counts.xlsx", row.names = T, col.names = T)

############################################ 5 group DEGS #############################################
setwd("~/Documents/Wistar/Haiyin/MeRLin/Data/Processed data/In-vivo/Bulk RNAseq/")
library(readxl)
library(tidyr)
library(stringr)
library(ggplot2)

tpm_counts <- read_xlsx("tpm_counts.xlsx", col_types = rep(x = c("text", "numeric"), times = c(1,17)))
tpm_vivo <- tpm_counts[,c(1, 4:17)]

DEGs <- read_excel(path = "~/Documents/Wistar/Haiyin/MeRLin/Data/Source data/MeRLin_signatures.xlsx")
Group_1_signature <- as.character(na.omit(DEGs$Group_1_signature))
Group_2_signature <- as.character(na.omit(DEGs$Group_2_signature))
Group_3_signature <- as.character(na.omit(DEGs$Group_3_signature))
Group_4_signature <- as.character(na.omit(DEGs$Group_4_signature))
Group_5_signature <- as.character(na.omit(DEGs$Group_5_signature))

plot_gene <- function(tpm, gene_list){
  for(i in seq_along(gene_list)){
    tpm_long <- pivot_longer(tpm[which(tpm$Symbol == gene_list[i]),], cols = colnames(tpm)[-1], names_to = "Sample", values_to = "Expr")
    tpm_long$Sample <- factor(x = tpm_long$Sample , levels = colnames(tpm)[-1])
    tpm_long$Group <- str_split_fixed(string = tpm_long$Sample, pattern = "_", n = 4)[,3]
    
    p <- ggplot(data = tpm_long, mapping = aes(x = Group, y = Expr, fill = Group)) +
      geom_boxplot(width = 0.7) +
      xlab(label = "") +
      ylab(label = "TPM experssion") +
      ggtitle(label = gene_list[i]) +
      theme_bw()
    
    tiff(filename = paste0("~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/RNAseq/Gene/",gene_list[i],".tiff"), width = 450, height = 330, res = 100)
    print(p)
    dev.off()
  }
}

plot_gene(tpm_vivo, Group_1_signature)
plot_gene(tpm_vivo, Group_2_signature)
plot_gene(tpm_vivo, Group_3_signature)
plot_gene(tpm_vivo, Group_4_signature)
plot_gene(tpm_vivo, Group_5_signature)

######################################### PCA and similarity ###########################################
setwd("~/Documents/Wistar/Haiyin/MeRLin/Data/Processed data/In-vivo/Bulk RNAseq/")
library(readxl)
library(tibble)
library(stringr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(pheatmap)

#load data
tpm_counts <- read_xlsx("tpm_counts.xlsx", col_types = rep(x = c("text", "numeric"), times = c(1,17)))
tpm_vitro <- tpm_counts[,c(1:3)]
tpm_vivo <- tpm_counts[,c(1, 4:17)]

#barplot
g_list <- c("MITF", "BRAF", "JUN", "SLC2A1", "SERPINE2", "VCL")
tpm_long <- pivot_longer(tpm_vitro[which(tpm_vitro$Symbol %in% g_list),], cols = colnames(tpm_vitro)[-1], names_to = "Sample", values_to = "Expr")
tpm_long$Symbol <- factor(x = tpm_long$Symbol , levels = c("MITF", "BRAF", "JUN", "SLC2A1", "SERPINE2", "VCL"))
tpm_long$Sample <- factor(x = tpm_long$Sample , levels = colnames(tpm_vitro)[-1])
tpm_long$Group <- str_split_fixed(string = tpm_long$Sample, pattern = "_", n = 2)[,1]


ggplot(data = tpm_long, mapping = aes(x = Sample, y = Expr, fill = Group)) +
  geom_bar(stat = "identity", width = 0.7) +
  xlab(label = "") +
  ylab(label = "") +
  facet_wrap(facets = ~Symbol, nrow = 2, scales = "free") +
  theme_bw()


g_list <- c("TYRP1", "JUN", "BRAF", "MKI67", "SLC2A1", "MCAM", "TIMP1", "VCL")
tpm_long <- pivot_longer(tpm_vivo[which(tpm_vitro$Symbol %in% g_list),], cols = colnames(tpm_vivo)[-1], names_to = "Sample", values_to = "Expr")
tpm_long$Symbol <- factor(x = tpm_long$Symbol , levels = c("TYRP1", "JUN", "BRAF", "MKI67", "SLC2A1", "MCAM", "TIMP1", "VCL"))
tpm_long$Sample <- factor(x = tpm_long$Sample , levels = colnames(tpm_vivo)[-1])
tpm_long$Group <- str_split_fixed(string = tpm_long$Sample, pattern = "_", n = 4)[,3]

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin//Figures/In-vivo/RNAseq/Gene/Selected genes.tiff", width = 2400, height = 1120, res = 200)
ggplot(data = tpm_long, mapping = aes(x = Group, y = Expr, fill = Group)) +
  geom_boxplot() +
  xlab(label = "") +
  ylab(label = "") +
  facet_wrap(facets = ~Symbol, nrow = 2, scales = "free") +
  theme_bw()
dev.off()

#reshape matrix
tpm_counts <- column_to_rownames(tpm_counts, var = "Symbol")
tpm_counts <- tpm_counts[-which(rowSums(tpm_counts) == 0),]

#PCA
pca <- prcomp(t(tpm_counts), scale=TRUE)
pca.var <- pca$sdev^2  
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)  
pca.data <- data.frame(Sample = rownames(pca$x),
                       PC1 = pca$x[,1],
                       PC2 = pca$x[,2],
                       Experiment = rep(x = c("In-vitro", "In-vivo"), times = c(2,15)),
                       Group = rep(x = c("CTRL", "DT", "D0", "D21", "D57", "EP"), times = c(1,1,3,3,3,6)))

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin//Figures/In-vivo/RNAseq/PCA.tiff", width = 2400, height = 1540, res = 200)
ggplot(data = pca.data, aes(x = PC1, y = PC2, color = Group, shape = Experiment)) +
  geom_point(size = 3) +
  geom_label_repel(aes(label = Sample), box.padding = 0.3) +
  theme_bw() +
  xlab(paste("PC1(",pca.var.per[1],"%"," variance)", sep="")) +
  ylab(paste("PC2(",pca.var.per[2],"%"," variance)", sep="")) +
  theme(legend.position = "none", plot.title=element_text(size = 14, hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 12), axis.title.x = element_text(size = 12)) +
  labs(title = "Sample relations based on all genes")
dev.off()

# similarity
library(proxy)

tpm_counts_select <- tpm_counts[c(Stress, NC, Lipid, PI3K, ECM, ME),]

cos_sim_matrix <- dist(t(tpm_counts_select), method = "cosine")
cos_sim_matrix <- as.matrix(1 - cos_sim_matrix)

pheatmap(mat = cos_sim_matrix, cluster_rows = F, cluster_cols = F,
         clustering_distance_rows = 'correlation', 
         clustering_distance_cols = 'correlation')


############################################## Signature ###############################################
setwd("~/Documents/Wistar/Haiyin/MeRLin/Data/Processed data/In-vivo/Bulk RNAseq/")
library(readxl)
library(tibble)
library(GSVA)
library(ggplot2)
library(pheatmap)
library(reshape2)
library(dplyr)
library(stringr)

# load tpm data
tpm <- read_xlsx("tpm_counts.xlsx", col_types = rep(x = c("text", "numeric"), times = c(1,17)))
tpm <- column_to_rownames(tpm, var = "Symbol")[,3:17]

# Marine signature
Pathway <- read_excel(path = "~/Documents/Wistar/Haiyin/MeRLin/Data/Source data/Marine_signatures.xlsx", sheet = "Sheet2")
Stem <- toupper(as.character(na.omit(Pathway$`Stem-like`)))
NCSC <- toupper(as.character(na.omit(Pathway$`Neural Crest-like`)))
Mesenchymal <- toupper(as.character(na.omit(Pathway$`Mesenchymal-like`)))
RNA <- toupper(as.character(na.omit(Pathway$`RNA processing`)))
AP <- toupper(as.character(na.omit(Pathway$`Antigen Presentation`)))
Stress <- toupper(as.character(na.omit(Pathway$`Stress-like`)))
Melanocytic <- toupper(as.character(na.omit(Pathway$Melanocytic)))

ssGSEA_Scores <- gsva(param = ssgseaParam(exprData = as.matrix(tpm), 
                                          geneSets = list("Stem-like" = Stem,
                                                          "Neural Crest-like" = NCSC,
                                                          "Mesenchymal-like" = Mesenchymal,
                                                          "RNA processing" = RNA,
                                                          "Antigen Presentation" = AP,
                                                          "Stress-like" = Stress,
                                                          "Melanocytic" = Melanocytic)))

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/RNAseq/Signature/Marine signature (Heatmap).tiff", width = 3300, height = 1650, res = 300)
pheatmap(mat = ssGSEA_Scores, fontsize = 14,
         color = colorRampPalette(c("#2600D1FF", "white", "#D60C00FF"))(256), 
         scale = "row", cluster_rows = T, cluster_cols = F, 
         show_colnames = T, angle_col = 315)
dev.off()

ssGSEA_Scores <- rownames_to_column(as.data.frame(ssGSEA_Scores), var = "Signature")
ssGSEA_Scores_long <- melt(ssGSEA_Scores, id.vars = "Signature", variable.name = "Sample", value.name = "Score") %>%
  mutate(Group = str_extract(Sample, "D0|D21|D57|EP"))

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/RNAseq/Signature/Marine signature (Boxplot).tiff", width = 3300, height = 2200, res = 300)
ggplot(ssGSEA_Scores_long, mapping = aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot() +
  facet_wrap(.~ Signature, scales = "free") +
  theme_classic()
dev.off()

# MeRLin signature
Signature <- read_excel(path = "~/Documents/Wistar/Haiyin/MeRLin/Data/Source data/MeRLin_signatures.xlsx", sheet = "Identified signatures")
Stress <- as.character(na.omit(Signature$Stress_like_signature))
NC <- as.character(na.omit(Signature$NC_like_signature))
Lipid <- as.character(na.omit(Signature$Lipid_metabolism))
PI3K <- as.character(na.omit(Signature$PI3K_signaling))
ECM <- as.character(na.omit(Signature$ECM_remodeling))
ME <- as.character(na.omit(Signature$Melanocytic_markers))

ssGSEA_Scores <- gsva(param = ssgseaParam(exprData = as.matrix(tpm), 
                                          geneSets = list("Stress-like signature" = Stress,
                                                          "Neural Crest-like" = NC,
                                                          "Lipid metabolism" = Lipid,
                                                          "PI3K signaling" = PI3K,
                                                          "ECM remodeling" = ECM,
                                                          "Melanocytic markers" = ME)))

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/RNAseq/Signature/MeRLin signature (Heatmap).tiff", width = 3300, height = 1650, res = 300)
pheatmap(mat = ssGSEA_Scores, fontsize = 14,
         color = colorRampPalette(c("#2600D1FF", "white", "#D60C00FF"))(256), 
         scale = "row", cluster_rows = T, cluster_cols = F, 
         show_colnames = T, angle_col = 315)
dev.off()

ssGSEA_Scores <- rownames_to_column(as.data.frame(ssGSEA_Scores), var = "Signature")
ssGSEA_Scores_long <- melt(ssGSEA_Scores, id.vars = "Signature", variable.name = "Sample", value.name = "Score") %>%
  mutate(Group = str_extract(Sample, "D0|D21|D57|EP"))
ssGSEA_Scores_long$Signature <- factor(ssGSEA_Scores_long$Signature, levels = c("Melanocytic markers", "Neural Crest-like", "Lipid metabolism",
                                                                                "PI3K signaling", "ECM remodeling", "Stress-like signature"))
ssGSEA_Scores_long$Group[which(ssGSEA_Scores_long$Group == "D0")] <- "Day 0"
ssGSEA_Scores_long$Group[which(ssGSEA_Scores_long$Group == "D21")] <- "Day 21"
ssGSEA_Scores_long$Group[which(ssGSEA_Scores_long$Group == "D57")] <- "Day 57"
ssGSEA_Scores_long$Group[which(ssGSEA_Scores_long$Group == "EP")] <- "Endpoint"

ssGSEA_Scores_long <- subset(ssGSEA_Scores_long, Signature %in% c("Stress-like signature", "Lipid metabolism", "PI3K signaling", "ECM remodeling"))
ssGSEA_Scores_long$Signature <- factor(x = ssGSEA_Scores_long$Signature, 
                                       levels = c("Stress-like signature", "Lipid metabolism", "PI3K signaling", "ECM remodeling"))

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/RNAseq/Signature/MeRLin signature (Boxplot).tiff", width = 4260, height = 1080, res = 300)
ggplot(ssGSEA_Scores_long, mapping = aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot() +
  facet_wrap(.~ Signature, scales = "free", nrow = 1) +
  theme_classic()
dev.off()

