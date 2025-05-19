################################################# Step1 ####################################################
setwd("~/Documents/Wistar/Haiyin/MeRLin/Data/Processed data/in-vivo/")
library(Seurat)

data <- readRDS(file ="EP_Clonocluster.rds")
data <- subset(data, bc_group %in% c("Barcode group 1", "Barcode group 2", "Barcode group 3", "Barcode group 4", "Barcode group 5"))

M <- as.matrix(GetAssayData(object = data, assay = "RNA", layer = "count"))
saveRDS(object = M, file = "inferCNV/EP_Matrix.rds")

G <- data.frame(Group = data$bc_group, row.names = colnames(data))
saveRDS(object = G, file = "inferCNV/EP_Group.rds")

################################################# Step2 ####################################################
setwd("~/Documents/Wistar/Haiyin/MeRLin/Data/Processed data/in-vivo/inferCNV/")
library(infercnv)

M <- readRDS("EP_Matrix.rds")
G <- readRDS("EP_Group.rds")
P <- readRDS("GRCh38_gene_pos_gene_id.rds")

infercnv_obj = CreateInfercnvObject(raw_counts_matrix = M, 
                                    annotations_file = G, 
                                    gene_order_file = P, 
                                    ref_group_names = NULL)

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff = 0.1,
                             out_dir = ".",
                             cluster_by_groups = FALSE,
                             denoise = TRUE,
                             HMM = TRUE,
                             analysis_mode = "subclusters", 
                             resume_mode = TRUE,
                             num_threads = 4)

################################################# Step3 ####################################################
setwd("~/Documents/Wistar/Haiyin/MeRLin/Data/Processed data/in-vivo/inferCNV/")
library(infercnv)
library(ggplot2)
library(tibble)
library(Seurat)
library(readxl)

infercnv_obj <- readRDS("run.final.infercnv_obj")

cnv_matrix <- round(infercnv_obj@expr.data, 5) - 1.00469
scores <- apply(infercnv_obj@expr.data, 2, function(x){sum(x < 0.9 | x > 1.1)/length(x)})

scores <- as.data.frame(scores)
scores <- rownames_to_column(scores, var = "cell")

G <- readRDS("EP_Group.rds")
G <- rownames_to_column(G, var = "cell")

df <- merge(x = G, y = scores, by = "cell")
df$Group <- factor(x = df$Group, levels = c("Barcode group 5", "Barcode group 4", "Barcode group 3", "Barcode group 2", "Barcode group 1"))
ggplot(data = df, mapping = aes(x = Group, y = scores, fill = Group)) +
  geom_boxplot(outliers = F) +
  scale_fill_manual(values = c("Barcode group 1" = "#E64B35FF",
                                "Barcode group 2" = "#3C5488FF", 
                                "Barcode group 3" = "#00A087FF", 
                                "Barcode group 4" = "#7E6148FF", 
                                "Barcode group 5" = "grey")) +
  xlab(label = "") +
  ylab(label = "CNV Scores") +
  coord_flip() +
  theme_classic() +
  theme(axis.ticks.length = unit(0.2, "cm"),
        axis.ticks.y = element_blank())

df2 <- merge(x = G, y = rownames_to_column(as.data.frame(t(cnv_matrix)), var = "cell"), by = "cell")
df2$Group <- factor(x = df2$Group, levels = c("Barcode group 5", "Barcode group 4", "Barcode group 3", "Barcode group 2", "Barcode group 1"))
ggplot(data = df2, mapping = aes(x = Group, y = CCND1, fill = Group)) +
  geom_boxplot(outliers = F) +
  scale_fill_manual(values = c("Barcode group 1" = "#E64B35FF",
                               "Barcode group 2" = "#3C5488FF", 
                               "Barcode group 3" = "#00A087FF", 
                               "Barcode group 4" = "#7E6148FF", 
                               "Barcode group 5" = "grey")) +
  xlab(label = "") +
  ylab(label = "") +
  coord_flip() +
  theme_classic() +
  theme(theme(legend.position = "none"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks.y = element_blank())

# corr
Signature <- read_excel(path = "~/Documents/Wistar/Haiyin/MeRLin/Data/Source data/MeRLin_signatures.xlsx")
Group_1_signature <- as.character(na.omit(Signature$Group_1_signature))
Group_2_signature <- as.character(na.omit(Signature$Group_2_signature))
Group_3_signature <- as.character(na.omit(Signature$Group_3_signature))
Group_4_signature <- as.character(na.omit(Signature$Group_4_signature))
Group_5_signature <- as.character(na.omit(Signature$Group_5_signature))

data <- readRDS(file = "../EP_Clonocluster.rds")
data <- subset(data, bc_group %in% c("Barcode group 1", "Barcode group 2", "Barcode group 3", "Barcode group 4", "Barcode group 5"))
expr_matrix <- as.matrix(GetAssayData(object = data, assay = "RNA", layer = "data"))

cal_corr <- function(gene_list, cnv_matrix, expr_matrix, group){
  res_list <- lapply(gene_list, function(gene){
    if(gene %in% rownames(expr_matrix) & gene %in% rownames(cnv_matrix)){
      cnv <- cnv_matrix[gene, names(which(data$bc_group == group))]
      expr <- expr_matrix[gene, names(which(data$bc_group == group))]
      temp <- cor.test(cnv, expr, method = "pearson")
      return(data.frame(gene = gene,
                        corr = as.numeric(temp$estimate),
                        p = temp$p.value))
    }
  })
  res_df <- do.call(rbind, res_list)
  return(res_df)
}

corr_1 <- cal_corr(Group_1_signature, cnv_matrix, expr_matrix, "Barcode group 1")
corr_2 <- cal_corr(Group_2_signature, cnv_matrix, expr_matrix, "Barcode group 2")
corr_3 <- cal_corr(Group_3_signature, cnv_matrix, expr_matrix, "Barcode group 3")
corr_4 <- cal_corr(Group_4_signature, cnv_matrix, expr_matrix, "Barcode group 4")
corr_5 <- cal_corr(Group_5_signature, cnv_matrix, expr_matrix, "Barcode group 5")

corr_1$group <- "Barcode group 1"
corr_2$group <- "Barcode group 2"
corr_3$group <- "Barcode group 3"
corr_4$group <- "Barcode group 4"
corr_5$group <- "Barcode group 5"

corr_all <- rbind(corr_1, corr_2, corr_3, corr_4, corr_5)
corr_all$group <- factor(x = corr_all$group, levels = c("Barcode group 5", "Barcode group 4", "Barcode group 3", "Barcode group 2", "Barcode group 1"))

ggplot(corr_all, aes(x = group, y = corr, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 0.7, alpha = 0.6) +
  scale_fill_manual(values = c("Barcode group 1" = "#E64B35FF",
                               "Barcode group 2" = "#3C5488FF",
                               "Barcode group 3" = "#00A087FF",
                               "Barcode group 4" = "#7E6148FF",
                               "Barcode group 5" = "grey")) +
  ylab("Pearson Correlation (CNV vs Expression)") +
  xlab("") +
  coord_flip() +
  theme_classic() +
  theme(legend.position = "none",
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks.y = element_blank())

library(pheatmap)
anno <- column_to_rownames(G, var = "cell")
anno_colors <- list(Group = c("Barcode group 1" = "#E64B35FF",
                              "Barcode group 2" = "#3C5488FF",
                              "Barcode group 3" = "#00A087FF",
                              "Barcode group 4" = "#7E6148FF",
                              "Barcode group 5" = "grey"),
                    Signature = c("Group 1 signature" = "#E64B35AA",
                                  "Group 2 signature" = "#3C5488AA",
                                  "Group 3 signature" = "#00A087AA",
                                  "Group 4 signature" = "#7E6148AA",
                                  "Group 5 signature" = "lightgrey"))

anno2 <- data.frame(gene = unique(c(Group_5_signature, Group_4_signature, Group_3_signature, Group_2_signature, Group_1_signature)),
                    Signature = rep(x = c("Group 5 signature", "Group 4 signature", "Group 3 signature", "Group 2 signature", "Group 1 signature"), 
                                    times = c(length(Group_5_signature), length(Group_4_signature)-1, length(Group_3_signature)-1, length(Group_2_signature)-10, length(Group_1_signature))))
anno2 <- subset(anno2, gene %in% rownames(cnv_matrix))
rownames(anno2) <- NULL
anno2 <- column_to_rownames(anno2, var = "gene")

cnv_matrix_heatmap <- cnv_matrix[rownames(anno2),]

n_colors <- 101
min_val <- -0.5
max_val <- 5
left <- seq(min_val, 0, length.out = floor(n_colors / 2) + 1)
right <- seq(0, max_val, length.out = floor(n_colors / 2) + 1)[-1]
breaks <- c(left, right)

pheatmap(mat = cnv_matrix_heatmap[rownames(anno2)[order(anno2$Signature)],rownames(anno)[order(anno$Group)]],
         fontsize = 6,
         color = colorRampPalette(c("#9797f8", "white", "#fe2f2f"))(100),
         breaks = breaks,
         scale = "column",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_col = anno,
         annotation_row = anno2,
         annotation_colors = anno_colors)

