setwd("~/Documents/Wistar/Haiyin/Barcode_v3/RNAseq/sam/")
library(data.table)
library(stringr)
library(Biostrings)

convert_quality_scores <- function(quality_string, offset = 33){
  quality_scores <- as.integer(charToRaw(quality_string)) - offset
  return(mean(quality_scores))
}

extraction_tool <- function(sam, seq = "CAGATCTTAGCCACTTTTTA"){
  df <- fread(file = sam, skip = 3)
  df <- subset(df, agrepl(pattern = seq, x = df$V10, max.distance = 1))
  match_positions <- unlist(lapply(df$V10, function(seq2){
    matches <- matchPattern(pattern = seq, subject = DNAString(seq2), max.mismatch = 1, with.indels = T)
    start(matches)
  }))

  df$barcode <- substr(x = df$V10, start = 1, stop = match_positions + nchar(seq) - 1)
  df$barcode_quality_score <- unlist(lapply(substr(x = df$V11, start = 1, stop = match_positions + 14), convert_quality_scores, offset = 33))
  df <- subset(df, barcode_quality_score == 40)
  
  df$len <- nchar(df$barcode)
  df <- df[order(df$len, decreasing = T),]
  
  View(as.data.frame(df$barcode))
}

extraction_tool(sam = "WM4237-1_EP_150_Unmapped.sam", seq = "GGTGGTGTTCGTGTTCAGATCTTAGCCACTTTTTA")














