############################################## In-vivo #################################################
setwd("~/Documents/Wistar/Haiyin/MeRLin/Data/Processed data/")
library(dplyr)
library(tibble)
library(ggplot2)
library(ggpattern)
library(ggsci)
`%ni%` <- Negate(`%in%`)

# load filtered data
CTRL <-  readRDS(file = "in-vitro/CTRL_filtered.rds")
DT <-  readRDS(file = "in-vitro/DT_filtered.rds")
EP <-  readRDS(file = "in-vivo/EP_filtered.rds")
D0 <-  readRDS(file = "in-vivo/D0_filtered.rds")
D21 <-  readRDS(file = "in-vivo/D21_filtered.rds")
EP <-  readRDS(file = "in-vivo/EP_filtered.rds")

# extract barcode information
CTRL_bc <- CTRL@meta.data %>% select(Barcode) %>% rownames_to_column(var = "Cell") %>% filter(!is.na(Barcode))
DT_bc <- DT@meta.data %>% select(Barcode) %>% rownames_to_column(var = "Cell") %>% filter(!is.na(Barcode))
D0_bc <- D0@meta.data %>% select(Barcode) %>% rownames_to_column(var = "Cell") %>% filter(!is.na(Barcode))
D21_bc <- D21@meta.data %>% select(Barcode) %>% rownames_to_column(var = "Cell") %>% filter(!is.na(Barcode))
EP_bc <- EP@meta.data %>% select(Barcode) %>% rownames_to_column(var = "Cell") %>% filter(!is.na(Barcode))

DT_bc <- DT_bc[which(DT_bc$Barcode %in% CTRL_bc$Barcode),]

df <- data.frame(Sample = factor(x = c("Control", "BRAFi/MEKi", "Day 0", "Day 21", "Endpoint"), 
                                 levels = c("Control", "BRAFi/MEKi", "Day 0", "Day 21", "Endpoint")),
                 Experiment = c("In-vitro", "In-vitro", "In-vivo", "In-vivo", "In-vivo"),
                 Total_cells = c(ncol(CTRL), ncol(DT), ncol(D0), ncol(D21), ncol(EP)),
                 Total_barcodes = c(length(CTRL_bc$Barcode), length(DT_bc$Barcode), length(D0_bc$Barcode), length(D21_bc$Barcode), length(EP_bc$Barcode)),
                 Unique_barcodes = c(length(unique(CTRL_bc$Barcode)), length(unique(DT_bc$Barcode)), length(unique(D0_bc$Barcode)), length(unique(D21_bc$Barcode)), length(unique(EP_bc$Barcode))))

tiff(filename = "../../Figures/In-vivo/scRNAseq/Barcode/Total cells.tiff", width = 1400, height = 1200, compression = "lzw", res = 300)
ggplot(data = df, mapping = aes(x = Sample, y = Total_cells, pattern = Sample)) +
  geom_bar_pattern(stat = "identity", colour = "black", fill = "white", position = "dodge", width = 0.8, show.legend = FALSE)+
  scale_pattern_manual(values = c("Control" = "none", "BRAFi/MEKi" = "stripe", "Day 0" = "none", "Day 21" = "stripe", "Endpoint" = "weave")) +
  facet_grid(. ~ Experiment, scales = "free_x", space = "free_x") +
  ggtitle(label = "Total cells") +
  xlab(label = "") +
  ylab(label = "") +
  theme_bw()
dev.off()

tiff(filename = "../../Figures/In-vivo/scRNAseq/Barcode/Total barcodes.tiff", width = 1400, height = 1200, compression = "lzw", res = 300)
ggplot(data = df, mapping = aes(x = Sample, y = Total_barcodes, pattern = Sample)) +
  geom_bar_pattern(stat = "identity", colour = "black", fill = "white", position = "dodge", width = 0.8, show.legend = FALSE)+
  scale_pattern_manual(values = c("Control" = "none", "BRAFi/MEKi" = "stripe", "Day 0" = "none", "Day 21" = "stripe", "Endpoint" = "weave")) +
  facet_grid(. ~ Experiment, scales = "free_x", space = "free_x") +
  ggtitle(label = "Total barcodes") +
  xlab(label = "") +
  ylab(label = "") +
  theme_bw()
dev.off()

tiff(filename = "../../Figures/In-vivo/scRNAseq/Barcode/Unique barcodes.tiff", width = 1400, height = 1200, compression = "lzw", res = 300)
ggplot(data = df, mapping = aes(x = Sample, y = Unique_barcodes, pattern = Sample)) +
  geom_bar_pattern(stat = "identity", colour = "black", fill = "white", position = "dodge", width = 0.8, show.legend = FALSE)+
  scale_pattern_manual(values = c("Control" = "none", "BRAFi/MEKi" = "stripe", "Day 0" = "none", "Day 21" = "stripe", "Endpoint" = "weave")) +
  facet_grid(. ~ Experiment, scales = "free_x", space = "free_x") +
  ggtitle(label = "Unique barcodes") +
  xlab(label = "") +
  ylab(label = "") +
  theme_bw()
dev.off()

tiff(filename = "../../Figures/In-vivo/scRNAseq/Barcode/Unique barcodes (In-vitro).tiff", width = 800, height = 1000, compression = "lzw", res = 300)
ggplot(data = subset(df, Experiment == "In-vitro"), mapping = aes(x = Sample, y = Unique_barcodes, pattern = Sample)) +
  geom_bar_pattern(stat = "identity", colour = "black", fill = "white", position = "dodge", width = 0.7, show.legend = FALSE)+
  scale_pattern_manual(values = c("Control" = "none", "BRAFi/MEKi" = "stripe")) +
  ggtitle(label = "Unique barcodes") +
  xlab(label = "") +
  ylab(label = "") +
  theme_bw()
dev.off()

# calculate barcode diversity
library(vegan)
CTRL_bc_div <- table(CTRL_bc$Barcode)
DT_bc_div <- table(DT_bc$Barcode)
D0_bc_div <- table(D0_bc$Barcode)
D21_bc_div <- table(D21_bc$Barcode)
EP_bc_div <- table(EP_bc$Barcode)

df <- data.frame(Sample = factor(x = c("Control", "BRAFi/MEKi", "Day 0", "Day 21", "Endpoint"), 
                                 levels = c("Control", "BRAFi/MEKi", "Day 0", "Day 21", "Endpoint")),
                 Experiment = c("In-vitro", "In-vitro", "In-vivo", "In-vivo", "In-vivo"),
                 Shannon = c(diversity(CTRL_bc_div, index = "shannon"), 
                             diversity(DT_bc_div, index = "shannon"), 
                             diversity(D0_bc_div, index = "shannon"), 
                             diversity(D21_bc_div, index = "shannon"), 
                             diversity(EP_bc_div, index = "shannon")), 
                 Pielou = c(diversity(CTRL_bc_div, index = "shannon")/log(length(CTRL_bc_div)), 
                            diversity(DT_bc_div, index = "shannon")/log(length(DT_bc_div)), 
                            diversity(D0_bc_div, index = "shannon")/log(length(D0_bc_div)), 
                            diversity(D21_bc_div, index = "shannon")/log(length(D21_bc_div)), 
                            diversity(EP_bc_div, index = "shannon")/log(length(EP_bc_div))),
                 Simpson = c(diversity(CTRL_bc_div, index = "invsimpson"), 
                             diversity(DT_bc_div, index = "invsimpson"), 
                             diversity(D0_bc_div, index = "invsimpson"), 
                             diversity(D21_bc_div, index = "invsimpson"), 
                             diversity(EP_bc_div, index = "invsimpson")), 
                 Margalef = c((length(CTRL_bc_div)-1)/log(sum(CTRL_bc_div)),
                              (length(DT_bc_div)-1)/log(sum(DT_bc_div)),
                              (length(D0_bc_div)-1)/log(sum(D0_bc_div)),
                              (length(D21_bc_div)-1)/log(sum(D21_bc_div)),
                              (length(EP_bc_div)-1)/log(sum(EP_bc_div))),
                 group = "BC")

tiff(filename = "../../Figures/In-vivo/scRNAseq/Barcode/Shannon–Wiener Diversity Index.tiff", width = 1500, height = 1200, compression = "lzw", res = 300)
ggplot(data = df, mapping = aes(x = Sample, y = Shannon, group = group)) +
  geom_point(size = 3.5, show.legend = FALSE) +
  geom_line(linewidth = 1.5) +
  ggtitle(label = "Shannon–Wiener Diversity Index") +
  xlab(label = "") +
  ylab(label = "") +
  theme_bw() +
  facet_grid(. ~ Experiment, scales = "free_x", space = "free_x") +
  theme(axis.text = element_text(size = 12),
        title = element_text(size = 14, face = "bold"))
dev.off()

tiff(filename = "../../Figures/In-vivo/scRNAseq/Barcode/Shannon–Wiener Diversity Index (In-vivo).tiff", width = 1500, height = 1200, compression = "lzw", res = 300)
ggplot(data = subset(df, Experiment == "In-vivo"), mapping = aes(x = Sample, y = Shannon, group = group)) +
  geom_point(size = 3.5, show.legend = FALSE) +
  geom_line(linewidth = 1.5) +
  ggtitle(label = "Shannon–Wiener Diversity Index") +
  xlab(label = "") +
  ylab(label = "") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        title = element_text(size = 14, face = "bold"))
dev.off()

tiff(filename = "../../Figures/In-vivo/scRNAseq/Barcode/Pielou Evenness Index.tiff", width = 1500, height = 1200, compression = "lzw", res = 300)
ggplot(data = df, mapping = aes(x = Sample, y = Pielou, group = group)) +
  geom_point(size = 3.5, show.legend = FALSE) +
  geom_line(size = 1.5) +
  ggtitle(label = "Pielou Evenness Index") +
  xlab(label = "") +
  ylab(label = "") +
  theme_bw() +
  facet_grid(. ~ Experiment, scales = "free_x", space = "free_x") +
  theme(axis.text = element_text(size = 12),
        title = element_text(size = 14, face = "bold"))
dev.off()

tiff(filename = "../../Figures/In-vivo/scRNAseq/Barcode/Inversed Simpson Diversity Index.tiff", width = 1500, height = 1200, compression = "lzw", res = 300)
ggplot(data = df, mapping = aes(x = Sample, y = Simpson, group = group)) +
  geom_point(size = 3.5, show.legend = FALSE) +
  geom_line(size = 1.5) +
  ggtitle(label = "Inversed Simpson Diversity Index") +
  xlab(label = "") +
  ylab(label = "") +
  theme_bw() +
  facet_grid(. ~ Experiment, scales = "free_x", space = "free_x") +
  theme(axis.text = element_text(size = 12),
        title = element_text(size = 14, face = "bold"))
dev.off()

tiff(filename = "../../Figures/In-vivo/scRNAseq/Barcode/Margalef Richness Index.tiff", width = 1500, height = 1200, compression = "lzw", res = 300)
ggplot(data = df, mapping = aes(x = Sample, y = Margalef, group = group)) +
  geom_point(size = 3.5, show.legend = FALSE) +
  geom_line(size = 1.5) +
  ggtitle(label = "Margalef Richness Index") +
  xlab(label = "") +
  ylab(label = "") +
  theme_bw() +
  facet_grid(. ~ Experiment, scales = "free_x", space = "free_x") +
  theme(axis.text = element_text(size = 12),
        title = element_text(size = 14, face = "bold"))
dev.off()

# barcode proportion in EP
data <- data.frame(
  group = factor(x = c("GATGATCAACAAGGA", "GATCTACCTCCAGAT", "GTTGAACGACCACAA", "GCAGGTGTTCAACGA", "CCTGGTCTAGGAGGT",
                       "Top6-Top10 barcodes", "Top11-Top20 barcodes", "Minor population", "Singletons"), 
                 levels = c("GATGATCAACAAGGA", "GATCTACCTCCAGAT", "GTTGAACGACCACAA", "GCAGGTGTTCAACGA", "CCTGGTCTAGGAGGT",
                            "Top6-Top10 barcodes", "Top11-Top20 barcodes", "Minor population", "Singletons")),
  value = c(291, 200, 192, 155, 147, 531, 369, 533, 147))

data <- data %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(data$value) *100) %>%
  mutate(ypos = cumsum(prop) - 0.5*prop)

my_colors <- c(pal_frontiers("default", alpha = 0.8)(8), "#808180CC")

tiff(filename = "../../Figures/In-vivo/scRNAseq/Barcode/Barcode proportion in EP.tiff", width = 2600, height = 1600, compression = "lzw", res = 300)
ggplot(data, aes(x = "", y = prop, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  scale_fill_manual(values = my_colors) +
  coord_polar("y", start = 0) +
  theme_void() + 
  geom_text(aes(y = ypos, label = paste0(round(prop,1), "%")), color = "white",size = 8, family = "Arial", fontface = "bold", position = position_nudge(x = 0.2)) +
  theme(text = element_text(family = "Arial", size = 20),
        legend.box.margin=margin(-15,-15,-15,-15),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"))
dev.off()

# merge barcode information
D0_div_df <- as.data.frame(D0_bc_div);colnames(D0_div_df) <- c("barcode", "D0")
D21_div_df <- as.data.frame(D21_bc_div);colnames(D21_div_df) <- c("barcode", "D21")
EP_div_df <- as.data.frame(EP_bc_div);colnames(EP_div_df) <- c("barcode", "EP")

merge_all <- function(x, y){merge(x, y, all = TRUE)}
div_all <- Reduce(merge_all, list(D0_div_df, D21_div_df, EP_div_df))
div_all[is.na(div_all)] = 0

D0_singleton <- div_all$barcode[which(div_all$D0 == 1 & div_all$D21 == 0 & div_all$EP == 0)];print(length(D0_singleton))
EP_singleton <- div_all$barcode[which(div_all$D0 == 0 & div_all$D21 <= 1 & div_all$EP <= 1)];print(length(EP_singleton))
even_barcode <- div_all$barcode[which(div_all$D0 == div_all$D21 & div_all$D21 == div_all$EP)];print(length(even_barcode))

div_all <- subset(div_all, barcode %ni% c(D0_singleton, EP_singleton, even_barcode))
div_all$D0_adj <- div_all$D0*nrow(D0_bc)/ncol(D0)
div_all$D21_adj <- div_all$D21*nrow(D21_bc)/ncol(D21)
div_all$EP_adj <- div_all$EP*nrow(EP_bc)/ncol(EP)

# test
library(minpack.lm)
library(scales)
library(cowplot)

## Exponential death model
rss <- c()
plot_list <- list()
plot_count <- 0

for(i in seq_along(div_all$barcode)) {
  data <- data.frame(
    x = c(0, 21, 94), 
    y = rescale(x = c(div_all$D0[i], div_all$D21[i], div_all$EP[i]), to = c(0, 1))
  )
  model <- nlsLM(y ~ exp(-k * x), data = data,
                 start = c(k = 0.01),
                 lower = c(k = 0),
                 upper = c(k = Inf),
                 control = nls.lm.control(maxiter = 50))
  rss[i] <- sum(residuals(model)^2)
  
  if(sum(residuals(model)^2) <= 0.1) {
    plot_data <- data.frame(x = seq(min(data$x), max(data$x), length.out = 100))
    plot_data$y <- predict(model, newdata = plot_data)
    
    p <- ggplot(data, aes(x = x, y = y)) + 
      geom_point(size = 3, color = "blue") + 
      geom_line(data = plot_data, aes(x = x, y = y), color = "red") +
      labs(x = "Time (days)", y = "Survival rate", title = div_all$barcode[i])
    
    plot_list[[length(plot_list) + 1]] <- p
    plot_count <- plot_count + 1
    
    if (plot_count == 30) {
      pp <- plot_grid(plotlist = plot_list, ncol = 5, align = 'v')
      print(pp)
      plot_list <- list()
      plot_count <- 0
    }
  }
}

if (length(plot_list) > 0) {
  plot_grid(plotlist = plot_list, ncol = 5, align = 'v')
}

sensitive_bc <- div_all$barcode[which(rss < 0.1)]

rss2 <- c()
plot_list <- list()
plot_count <- 0

for(i in seq_along(div_all$barcode)) {
  data <- data.frame(
    x = c(0, 21, 94), 
    y = rescale(x = c(div_all$D0[i], div_all$D21[i], div_all$EP[i]), to = c(0, 1))
  )
  model <- nlsLM(y ~ P0 * exp(r * x), data = data,
                 start = c(P0 = 0.1, r = 0.01),
                 lower = c(P0 = 0, r = 0),
                 upper = c(P0 = Inf, r = Inf),
                 control = nls.lm.control(maxiter = 50))
  rss2[i] <- sum(residuals(model)^2)
  
  if(sum(residuals(model)^2) <= 0.1) {
    plot_data <- data.frame(x = seq(min(data$x), max(data$x), length.out = 100))
    plot_data$y <- predict(model, newdata = plot_data)
    
    p <- ggplot(data, aes(x = x, y = y)) + 
      geom_point(size = 3, color = "blue") + 
      geom_line(data = plot_data, aes(x = x, y = y), color = "red") +
      labs(x = "Time (days)", y = "Survival rate", title = div_all$barcode[i])
    
    plot_list[[length(plot_list) + 1]] <- p
    plot_count <- plot_count + 1
    
    if (plot_count == 25) {
      pp <- plot_grid(plotlist = plot_list, ncol = 5, align = 'v')
      print(pp)
      plot_list <- list()
      plot_count <- 0
    }
  }
}

persister_bc <- div_all$barcode[which(rss2 < 0.1)]

multi_fate_bc <- setdiff(div_all$barcode, c(D0_singleton, EP_singleton, even_barcode, sensitive_bc, persister_bc))

# plot
div_all <- Reduce(merge_all, list(D0_div_df, D21_div_df, EP_div_df))
div_all[is.na(div_all)] = 0

div_all$`Barcode type` <- NA
div_all$`Barcode type`[which(div_all$barcode %in% D0_singleton)] = "Multi-fate"
div_all$`Barcode type`[which(div_all$barcode %in% EP_singleton)] = "Multi-fate"
div_all$`Barcode type`[which(div_all$barcode %in% even_barcode)] = "Multi-fate"
div_all$`Barcode type`[which(div_all$barcode %in% sensitive_bc)] = "Sensitive"
div_all$`Barcode type`[which(div_all$barcode %in% persister_bc)] = "Persister"
div_all$`Barcode type`[which(div_all$barcode %in% multi_fate_bc)] = "Multi_fate"

#GATGCTGGACTTGCT	73	47	93	multi_fate_bc	persister_bc
div_all$`Barcode type`[which(div_all$barcode == "GATGCTGGACTTGCT")] <- "Persister"

div_all_long <- reshape2::melt(div_all, id.vars = c("barcode", "Barcode type"), variable.name = "timepoint", value.name = "value") %>%
  mutate(timepoint = recode(timepoint, "D0" = "Day 0")) %>%
  mutate(timepoint = recode(timepoint, "D21" = "Day 21")) %>%
  mutate(timepoint = recode(timepoint, "EP" = "Endpoint"))

div_all_long$`Barcode type` <- factor(div_all_long$`Barcode type`, levels = c("Persister", "Sensitive", "Multi-fate"))

tiff(filename = "../../Figures/In-vivo/scRNAseq/Barcode/Each Barcode Across Timepoints.tiff", width = 2000, height = 1000, compression = "lzw", res = 300)
ggplot(div_all_long, aes(x = timepoint, y = value, group = barcode, color = `Barcode type`)) +
  geom_line() +
  facet_wrap(~`Barcode type`, ncol = 3) +
  labs(title = "Each Barcode Across Timepoints",
       x = "", y = "Cell Count") +
  theme_bw()
dev.off()

WriteXLS::WriteXLS(x = div_all, ExcelFileName = "div_all_15bp.xlsx")

div_summary <- div_all_long %>%
  group_by(timepoint, `Barcode type`) %>%
  summarise(total = sum(value), .groups = 'drop')
div_summary$`Barcode type` <- factor(x = div_summary$`Barcode type`, levels = c("Sensitive", "Multi-fate", "Persister"))

D0_total <- sum(741, 2324, 990);D21_total <- sum(703, 459, 892);EP_total <- sum(2246, 80, 239);
div_summary$total <- c(741/D0_total, 2324/D0_total, 990/D0_total, 703/D21_total, 459/D21_total, 892/D21_total, 2246/EP_total, 80/EP_total, 239/EP_total)

ggplot(div_summary, aes(x = timepoint, y = total, fill = `Barcode type`, group = `Barcode type`)) +
  geom_area() +
  scale_fill_manual(values = c("Persister" = "#F8766D", "Sensitive" = "#00BA38", "Multi-fate" = "#619CFF")) +
  labs(title = "Stacked Area Chart of Barcode by Type", 
       x = "", y = "Cell (%)") +
  theme_bw() +
  theme(legend.title = element_text(size = 12), legend.position = "right")

############################################## In-vitro ################################################
setwd("~/Documents/Wistar/Haiyin/MeRLin/Data/Original data/In-vivo/Bulk RNAseq/")
library(readxl)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggpattern)
library(ggsci)
library(vegan)
`%ni%` <- Negate(`%in%`)

data <- read_excel("~/Documents/Wistar/Haiyin/Barcode_v3/RNAseq/barcode/div_all_15bp_v2.xlsx")
sample_names <- colnames(data)[colnames(data) != "barcode"]
data[is.na(data)] <- 0

diversity_list <- lapply(sample_names, function(sample) {
  sample_table <- table(rep(data$barcode, times = data[[sample]]))
  
  shannon <- diversity(sample_table, index = "shannon")
  simpson <- diversity(sample_table, index = "invsimpson")
  pielou <- shannon / log(length(sample_table))
  margalef <- (length(sample_table) - 1) / log(sum(sample_table))
  
  return(data.frame(Sample = sample,
                    Shannon = shannon,
                    Pielou = pielou,
                    Simpson = simpson,
                    Margalef = margalef))
})

df <- do.call(rbind, diversity_list)

df$Group <- case_when(
  grepl("^D0_", df$Sample) ~ "Day 0",
  grepl("^D21_", df$Sample) ~ "Day 21",
  grepl("^D57_", df$Sample) ~ "Day 57",
  grepl("^EP_", df$Sample) ~ "Endpoint",
  TRUE ~ "Other"
)

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/RNAseq/Barcode/Shannon Index.tiff", width = 1500, height = 1050, res = 300)
ggplot(df, aes(x = Group, y = Shannon, fill = Group)) +
  geom_boxplot() +
  theme_classic(base_size = 14) +
  labs(title = "Shannon Diversity by Group",
       x = NULL, y = "Shannon Index") +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "right")
dev.off()

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/RNAseq/Barcode/Pielou Index.tiff", width = 1500, height = 1050, res = 300)
ggplot(df, aes(x = Group, y = Pielou, fill = Group)) +
  geom_boxplot() +
  theme_classic(base_size = 14) +
  labs(title = "Pielou Evenness by Group",
       x = NULL, y = "Pielou Index") +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "right")
dev.off()

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/RNAseq/Barcode/Inversed Simpson Index.tiff", width = 1500, height = 1050, res = 300)
ggplot(df, aes(x = Group, y = Simpson, fill = Group)) +
  geom_boxplot() +
  theme_classic(base_size = 14) +
  labs(title = "Inversed Simpson Diversity by Group",
       x = NULL, y = "Inversed Simpson Index") +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "right")
dev.off()

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/RNAseq/Barcode/Margalef Index.tiff", width = 1500, height = 1050, res = 300)
ggplot(df, aes(x = Group, y = Margalef, fill = Group)) +
  geom_boxplot() +
  theme_classic(base_size = 12) +
  labs(title = "Margalef Richness by Group",
       x = NULL, y = "Margalef Index") +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "right")
dev.off()

# test
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(minpack.lm)

data <- read_excel("~/Documents/Wistar/Haiyin/Barcode_v3/RNAseq/barcode/div_all_15bp_v2.xlsx")
data[is.na(data)] <- 0

tpm <- readxl::read_xlsx(path = "~/Documents/Wistar/Haiyin/MeRLin/Data/Processed data/in-vivo/Bulk RNAseq/tpm_counts.xlsx")

df <- data.frame(Sample = c("D0_126", "D0_130", "D0_134", "D21_139", "D21_140", "D21_142", 
                            "D57_137", "D57_141", "D57_149", "EP_144", "EP_146", "EP_150", "EP_127", "EP_129", "EP_138"),
                 Timepoint = c("Day 0", "Day 0", "Day 0", "Day 21", "Day 21", "Day 21", "Day 57", "Day 57", "Day 57", "Endpoint", "Endpoint", "Endpoint", "Endpoint", "Endpoint", "Endpoint"),
                 Total_barcodes = c(sum(data$D0_126), sum(data$D0_130), sum(data$D0_134), sum(data$D21_139), sum(data$D21_140), sum(data$D21_142),
                                    sum(data$D57_137), sum(data$D57_141), sum(data$D57_149), sum(data$EP_144), sum(data$EP_146), sum(data$EP_150), sum(data$EP_127), sum(data$EP_129), sum(data$EP_138)),
                 Non_zero_genes = c(sum(tpm$`WM4237-1_BCv3_D0_126` != 0), sum(tpm$`WM4237-1_BCv3_D0_130` != 0), sum(tpm$`WM4237-1_BCv3_D0_134` != 0), 
                                    sum(tpm$`WM4237-1_BCv3_D21_139` != 0), sum(tpm$`WM4237-1_BCv3_D21_140` != 0), sum(tpm$`WM4237-1_BCv3_D21_142` != 0), 
                                    sum(tpm$`WM4237-1_BCv3_D57_137` != 0), sum(tpm$`WM4237-1_BCv3_D57_141` != 0), sum(tpm$`WM4237-1_BCv3_D57_149` != 0), 
                                    sum(tpm$`WM4237-1_BCv3_EP_144` != 0), sum(tpm$`WM4237-1_BCv3_EP_146` != 0), sum(tpm$`WM4237-1_BCv3_EP_150` != 0),
                                    sum(tpm$`WM4237-1_BCv3_EP_127` != 0), sum(tpm$`WM4237-1_BCv3_EP_129` != 0), sum(tpm$`WM4237-1_BCv3_EP_138` != 0)),
                 Active_genes = c(sum(tpm$`WM4237-1_BCv3_D0_126` > 35), sum(tpm$`WM4237-1_BCv3_D0_130` > 35), sum(tpm$`WM4237-1_BCv3_D0_134` > 35), 
                                  sum(tpm$`WM4237-1_BCv3_D21_139` > 35), sum(tpm$`WM4237-1_BCv3_D21_140` > 35), sum(tpm$`WM4237-1_BCv3_D21_142` > 35), 
                                  sum(tpm$`WM4237-1_BCv3_D57_137` > 35), sum(tpm$`WM4237-1_BCv3_D57_141` > 35), sum(tpm$`WM4237-1_BCv3_D57_149` > 35), 
                                  sum(tpm$`WM4237-1_BCv3_EP_144` > 35), sum(tpm$`WM4237-1_BCv3_EP_146` > 35), sum(tpm$`WM4237-1_BCv3_EP_150` > 35),
                                  sum(tpm$`WM4237-1_BCv3_EP_127` > 35), sum(tpm$`WM4237-1_BCv3_EP_129` > 35), sum(tpm$`WM4237-1_BCv3_EP_138` > 35)))

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/RNAseq/Non_zero_genes.tiff", width = 950, height = 700, res = 200)
ggplot(data = df, mapping = aes(x = Total_barcodes, y = Non_zero_genes,)) +
  geom_point(mapping = aes(color = Timepoint)) +
  geom_smooth(method = "lm", se = F) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  theme_bw()
dev.off()

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/RNAseq/Active_genes.tiff", width = 950, height = 700, res = 200)
ggplot(data = df, mapping = aes(x = Total_barcodes, y = Active_genes)) +
  geom_point(mapping = aes(color = Timepoint)) +
  geom_smooth(method = "lm", se = F) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  theme_bw()
dev.off()

expr <-  readxl::read_xlsx(path = "~/Documents/Wistar/Haiyin/MeRLin/Data/Processed data/in-vivo/Bulk RNAseq/tpm_counts.xlsx")
exprMat <- as.matrix(expr[,c(2:10,14:16)])
rownames(exprMat) <- expr$Symbol

nonzero_per_sample <- apply(exprMat, 2, function(x) sum(x != 0))
total_nonzero_genes <- sum(rowSums(exprMat != 0) > 0)
ratio <- nonzero_per_sample / total_nonzero_genes


data$D0_126 <- data$D0_126*1*1.3494141
data$D0_130 <- data$D0_130*1*1.5082067
data$D0_134 <- data$D0_134*1*1.4058427
data$D21_139 <- data$D21_139*0.42*0.8337355
data$D21_140 <- data$D21_140*0.36*0.8550690
data$D21_142 <- data$D21_142*0.4*0.8794547
data$D57_137 <- data$D57_137*0.43*0.9978267
data$D57_141 <- data$D57_141*0.49*0.9992620
data$D57_149 <- data$D57_149*0.48*1.0495041
data$EP_144 <- data$EP_144*1.19*0.8660453
data$EP_146 <- data$EP_146*0.71*1.1975728
data$EP_150 <- data$EP_150*1.57*0.8746964
data$EP_127 <- data$EP_127*0.89*
data$EP_129 <- data$EP_129*1.69*
data$EP_138 <- data$EP_138*0.71*

sample_names <- colnames(data)[colnames(data) != "barcode"]

time_map <- case_when(
  grepl("^D0_", sample_names) ~ "Day 0",
  grepl("^D21_", sample_names) ~ "Day 21",
  grepl("^D57_", sample_names) ~ "Day 57",
  grepl("^EP_", sample_names) ~ "Endpoint"
)
names(time_map) <- sample_names

long_df <- data %>%
  pivot_longer(-barcode, names_to = "Sample", values_to = "Abundance") %>%
  mutate(Time_point = time_map[Sample])

avg_df <- long_df %>%
  group_by(barcode, Time_point) %>%
  summarise(mean_abundance = mean(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Time_point, values_from = mean_abundance, values_fill = 0)

colnames(avg_df)[which(colnames(avg_df) %in% c("Day 0", "Day 21", "Day 57", "Endpoint"))] <- c("Day 0", "Day 21", "Day 57", "Endpoint")

rss1 <- rss2 <- rep(Inf, nrow(avg_df))
sensitive_bc <- persister_bc <- c()

for (i in 1:nrow(avg_df)) {
  vec <- as.numeric(avg_df[i, c("Day 0", "Day 21", "Day 57", "Endpoint")])
  x <- c(0, 21, 57, 94)
  y <- rescale(vec, to = c(0, 1))
  data_model <- data.frame(x = x, y = y)
  
  try({
    model1 <- nlsLM(y ~ exp(-k * x), data = data_model,
                    start = c(k = 0.01), lower = c(k = 0), upper = c(k = Inf),
                    control = nls.lm.control(maxiter = 50))
    rss1[i] <- sum(residuals(model1)^2)
  }, silent = TRUE)
  
  try({
    model2 <- nlsLM(y ~ P0 * exp(r * x), data = data_model,
                    start = c(P0 = 0.1, r = 0.01), lower = c(P0 = 0, r = 0),
                    upper = c(P0 = Inf, r = Inf),
                    control = nls.lm.control(maxiter = 50))
    rss2[i] <- sum(residuals(model2)^2)
  }, silent = TRUE)
}

avg_df$Type <- "Multi-fate"
avg_df$Type[rss1 < 0.1] <- "Sensitive"
avg_df$Type[rss2 < 0.1] <- "Persister"
avg_df$Type[which(c(avg_df$D0 + avg_df$D21 + avg_df$D57 + avg_df$EP) <= 1)] <- "Multi-fate"

avg_df$Type[which(avg_df$barcode %in% c("GATGATCAACAAGGA",
                                        "GATGGTGTTGAAGGA",
                                        "GATGATCATGTTCTA",
                                        "CCTGTACATGGTGGA",
                                        "GAACAACTAGAAGAT",
                                        "GATGGTCAAGCAGGA",
                                        "CTAGTAGCTGGAGCC",
                                        "GTAGGTGCTCATGAT",
                                        "GTTCTAGATGGTGTT",
                                        "CATCCAGTACAACTA",
                                        "CATGATCTAGGTGGT",
                                        "GAACAACAAGGAGCA",
                                        "GTTGGACTAGATGCT",
                                        "GCAGCAGTAGTAGTA",
                                        "GTTCTAGTAGATGGT"))] <- "Persister"

avg_df$Type <- factor(avg_df$Type, levels = c("Sensitive", "Persister", "Multi-fate"))

plot_df <- avg_df %>%
  select(barcode, Type, `Day 0`, `Day 21`, `Day 57`, Endpoint) %>%
  pivot_longer(cols = c("Day 0", "Day 21", "Day 57", "Endpoint"), names_to = "Time_point", values_to = "Abundance")

plot_df$Type <- factor(plot_df$Type, levels = c("Persister", "Sensitive", "Multi-fate"))

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/RNAseq/Barcode/Each Barcode Across Timepoints.tiff", width = 2250, height = 1000, compression = "lzw", res = 300)
ggplot(plot_df, aes(x = Time_point, y = Abundance, group = barcode, color = Type)) +
  geom_line() +
  scale_fill_manual(values = c("Persister" = "#F8766D", "Sensitive" = "#00BA38", "Multi-fate" = "#619CFF")) +
  facet_wrap(~Type, ncol = 3) +
  labs(title = "Each Barcode Across Timepoints",
       x = "", y = "Barcode Count") +
  theme_bw()
dev.off()

div_summary <- plot_df %>%
  group_by(Time_point, Type) %>%
  summarise(total = sum(Abundance), .groups = 'drop')
div_summary$Type <- factor(x = div_summary$Type, levels = c("Sensitive", "Multi-fate", "Persister"))

div_summary$total[which(div_summary$Time_point == "Day 0")] <- div_summary$total[which(div_summary$Time_point == "Day 0")]/12
div_summary$total[which(div_summary$Time_point == "Day 21")] <- div_summary$total[which(div_summary$Time_point == "Day 21")]/9
div_summary$total[which(div_summary$Time_point == "Day 57")] <- div_summary$total[which(div_summary$Time_point == "Day 57")]/6
div_summary$total[which(div_summary$Time_point == "Endpoint")] <- div_summary$total[which(div_summary$Time_point == "Endpoint")]/3

#div_summary$total[which(div_summary$Time_point == "D0")] <- div_summary$total[which(div_summary$Time_point == "D0")]/mean(1.3494141, 1.5082067, 1.4058427)
#div_summary$total[which(div_summary$Time_point == "D21")] <- div_summary$total[which(div_summary$Time_point == "D21")]/mean(0.8337355, 0.8550690, 0.8794547)
#div_summary$total[which(div_summary$Time_point == "D57")] <- div_summary$total[which(div_summary$Time_point == "D57")]/mean(0.9978267, 0.9992620, 1.0495041)
#div_summary$total[which(div_summary$Time_point == "EP")] <- div_summary$total[which(div_summary$Time_point == "EP")]/mean(0.8660453, 1.1975728, 0.8746964)

tiff(filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/RNAseq/Barcode/Stacked Barcode Across Timepoints.tiff", width = 1600, height = 1000, compression = "lzw", res = 300)
ggplot(div_summary, aes(x = Time_point, y = total, fill = Type, group = Type)) +
  geom_area() +
  scale_fill_manual(values = c("Persister" = "#F8766D", "Sensitive" = "#00BA38", "Multi-fate" = "#619CFF")) +
  labs(title = "Stacked Area Chart of Barcode by Type", 
       x = "", y = "Pseudo Barcode Count") +
  theme_bw() +
  theme(legend.title = element_text(size = 12), legend.position = "right")
dev.off()


