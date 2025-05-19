################################################# Step1 ####################################################

setwd("/Volumes/herlynm/linux/ychen/Geneformer/WM4237/h5ad/")
library(Seurat)
library(SeuratDisk)

data <- readRDS(file = "~/Documents/Wistar/Haiyin/MeRLin/Data/Processed data/in-vivo/EP_Clonocluster.rds")
data <- subset(data, bc_group %in% c("Barcode group 1", "Barcode group 2", "Barcode group 3", "Barcode group 4", "Barcode group 5"))

group <- as.data.frame(data$bc_group)
colnames(group) <- "bc_group"

options(Seurat.object.assay.version = "v3")
data0 <- Read10X(data.dir = "/Volumes/herlynm/linux/ychen/Haiyin/barcode_v3/scRNAseq/EP_DT/outs/filtered_feature_bc_matrix/", gene.column = 1)
data0 <- CreateSeuratObject(data0)
data <- Read10X(data.dir = "/Volumes/herlynm/linux/ychen/Haiyin/barcode_v3/scRNAseq/EP_DT/outs/filtered_feature_bc_matrix/", gene.column = 2)
data <- CreateSeuratObject(data)
data <- subset(x = data, cells = rownames(group))
data <- subset(x = data, features = setdiff(rownames(data), "barcode"))
data0 <- subset(x = data0, features = setdiff(rownames(data0), "barcode"))
data <- AddMetaData(data, metadata = group)
data[["RNA"]]@meta.features$ensembl_id <- rownames(data0)
data@meta.data$orig.ident <- NULL
colnames(data@meta.data)[colnames(data@meta.data) == "nCount_RNA"] <- "n_counts"
colnames(data@meta.data)[colnames(data@meta.data) == "nFeature_RNA"] <- "n_genes_by_counts"

SaveH5Seurat(object = data, filename = "EP.h5seurat", overwrite = T)
Convert(source = "EP.h5seurat", dest = "h5ad", assay = "RNA", overwrite = T)

################################################# Step2 ####################################################
from geneformer import TranscriptomeTokenizer

tk = TranscriptomeTokenizer({"group": "group"}, nproc=4)
tk.tokenize_data(data_directory = "/wistar/herlynm/ychen/Geneformer/WM4237/h5ad/", 
                 output_directory = "/wistar/herlynm/ychen/Geneformer/WM4237/token/", 
                 output_prefix = "EP",
                 file_format = "h5ad")

################################################# Step4 ####################################################
from geneformer import InSilicoPerturber
from geneformer import EmbExtractor

import logging
logging.basicConfig(level=logging.DEBUG)

cell_states_to_model = {"state_key": "group",
  "start_state": "Persister",
  "goal_state": "Sensitive"}

filter_data_dict={"group":["Persister","Sensitive"]}

embex = EmbExtractor(
  model_type = "CellClassifier",
  filter_data = filter_data_dict,
  emb_layer = -1,
  summary_stat = "exact_mean",
  forward_batch_size = 128,
  nproc = 4,
  token_dictionary_file = "/wistar/herlynm/ychen/Geneformer/Geneformer/geneformer/token_dictionary_gc95M.pkl")

state_embs_dict = embex.get_state_embs(
  cell_states_to_model = cell_states_to_model,
  model_directory = "/wistar/herlynm/ychen/Geneformer/Geneformer/gf-12L-95M-i4096_CLcancer",
  input_data_file = "/wistar/herlynm/ychen/Geneformer/WM4237/token/EP.dataset",
  output_directory = "/wistar/herlynm/ychen/Geneformer/WM4237/embedding/",
  output_prefix = "EP")

isp = InSilicoPerturber(
  perturb_type = "delete",
  perturb_rank_shift = None,
  genes_to_perturb = "all",
  cell_states_to_model = cell_states_to_model,
  state_embs_dict = state_embs_dict,
  combos = 0,
  anchor_gene = None,
  model_type = "CellClassifier",
  emb_mode = "cls",
  cell_emb_style = "mean_pool",
  filter_data = filter_data_dict,
  forward_batch_size = 128,
  max_ncells = None,
  emb_layer = -1,
  nproc = 1)

isp.perturb_data(
  model_directory = "/wistar/herlynm/ychen/Geneformer/Geneformer/gf-12L-95M-i4096_CLcancer",
  input_data_file = "/wistar/herlynm/ychen/Geneformer/WM4237/token/EP.dataset",
  output_directory = "/wistar/herlynm/ychen/Geneformer/WM4237/perturbation/",
  output_prefix = "EP")

################################################# Step3 ####################################################
from geneformer import InSilicoPerturberStats

cell_states_to_model = {"state_key": "bc_group",
  "start_state": "Persister",
  "goal_state": "Sensitive"}

ispstats = InSilicoPerturberStats(mode = "goal_state_shift",
                                  combos = 0,
                                  anchor_gene = None,
                                  token_dictionary_file = "/wistar/herlynm/ychen/Geneformer/Geneformer/geneformer/token_dictionary_gc95M.pkl",
                                  cell_states_to_model = cell_states_to_model)

ispstats.get_stats("/wistar/herlynm/ychen/Geneformer/WM4237/perturbation/",
                   None,
                   "/wistar/herlynm/ychen/Geneformer/WM4237/stats/",
                   "EP")

################################################# Step4 ####################################################
library(VennDiagram)
library(RColorBrewer)

data <- read_excel("~/Documents/Wistar/Haiyin/MeRLin/Data/Source data/EP_DEGs_persisters_vs_group5.xlsx", "Sheet1")

set1 <- na.omit(data$DEGs)
set2 <- na.omit(data$Geneformer)
set3 <- na.omit(data$Druggable)

Reduce(f = intersect, x = list(set1, set2, set3))

myCol <- brewer.pal(3, "Pastel2")

venn.diagram(
  x = list(set1, set2, set3),
  category.names = c("DEGs" , "Geneformer" , "Druggable"),
  filename = "~/Documents/Wistar/Haiyin/MeRLin/Figures/In-vivo/scRNAseq/Geneformer/Venn.tiff",
  output=TRUE,
  
  # Output features
  imagetype="tiff" ,
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

