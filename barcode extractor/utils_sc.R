library(GenomicAlignments)
library(stringr)
library(stringdist)
library(Biostrings)
is.not.na <- Negate(is.na)

#create strong-random-weak regular expression for barcode correction
generate_pattern <- function(len){
  pattern <- paste0(rep(x = "[CG][ATCG][AT]", times = floor(len/3)), collapse = "")
  if(len %% 3 == 1){
    pattern <- paste0("[AT]", pattern)
  }else if(len %% 3 == 2){
    pattern <- paste0("[ATCG][AT]", pattern)
  }else{
    pattern <- pattern
  }
  return(paste0("^", pattern, "$"))
}

#calculate hamming distances between a given string to another set of strings
calculate_hamming_distances <- function(target_string, strings){
  if (all(nchar(target_string) == nchar(strings))){
    distances <- sapply(strings, function(x) stringdist(x, target_string, method = "hamming"))
    return(distances)
  }else{
    stop("Strings must be of equal length.")
  }
}

#allowed mismatch ratio
allowed_mismatch <- function(string, ratio = 0.033){
  return(ceiling(nchar(string)*ratio))
}

extract_barcode <- function(bam, len = 12, known_seq = "CAGATCTTAGCCACTTTTTA"){
  #load bam
  param <- ScanBamParam(what = "seq", tag = c("AS", "nM", "RG", "TX", "CR", "CB", "UR", "UB"))
  index <- str_replace(string = bam, pattern = "bam$", replacement = "bam.bai")
  df <- as.data.frame(readGAlignments(file = bam, index = index, param = param, use.names = TRUE))
  
  #extract + strand reads and remove cells with empty CB tag
  df <- subset(df, strand == "+" & is.not.na(CB))
  
  #match the known sequence on the 3 prime end
  df <- subset(df, agrepl(pattern = known_seq, x = df$seq, max.distance = allowed_mismatch(string = known_seq)))
  
  #find the end position of barcodes
  df$match_positions <- unlist(lapply(df$seq, function(seq){
    matches <- matchPattern(pattern = known_seq, subject = DNAString(seq), max.mismatch = 1, with.indels = T)
    start(matches)
  }))
  
  #remove short barcode
  df <- subset(df, match_positions > len)

  #remove UMI duplicates
  df <- df[order(df$UB, df$match_positions, decreasing = T),]
  df <- df[-which(duplicated(paste0(df$CB, df$UB))),]
  
  #truncate barcodes to specified length and calculate quality scores
  df$barcode <- substr(x = df$seq, start = df$match_positions - len, stop = df$match_positions - 1)

  #remove empty vector
  empty_vector <- "TCCAGCGGACCTTCCTTCCCGCGGCCTGCTGCCGGCTCTGCGGCCTCTTCCGCGTCTTCGCCTTCGCCCTCAGACGAGTCGGATCTCCCTTTGGGCCGCCTCCCCGCCTGGTACTTTTAAGACCAATGACTTACAAGGCAGCTGTGGTAC"
  truncated_empty_vector <- substring(text = empty_vector, first = nchar(empty_vector) - len + 1)
  empty_vector_index <- agrep(pattern = truncated_empty_vector, x = df$barcode, max.distance = allowed_mismatch(string = truncated_empty_vector))
  df <- df[-empty_vector_index,]
  print(paste("Number of empty vector", length(empty_vector_index), sep = ": "))
  
  #correct barcode
  pattern <- generate_pattern(len)
  whitelist <- unique(grep(pattern, df$barcode, value = T))
  blacklist <- setdiff(df$barcode, whitelist)
  
  for(i in blacklist){
    h_dist <- calculate_hamming_distances(i, whitelist)
    if(min(h_dist) == 1){
      index <- which(h_dist == 1)
      if(length(index) == 1){
        df$barcode[which(df$barcode == i)] <- whitelist[index]
      }else{
        df <- df[-which(df$barcode == i),]
      }
    }else{
      if(length(which(df$barcode == i)) > 10){
        df <- df
      }else{
        df <- df[-which(df$barcode == i),]
      }
    }
  }
  
  #remove cells with multiple barcodes
  cell <- unique(df$CB)
  ind <- c()
  
  for(i in cell){
    temp <- df[which(df$CB == i),]
    len_bc <- length(unique(temp$barcode))
    if(len_bc >= 3){
      freq <- table(temp$barcode)/nrow(temp)
      if(any(freq >= 0.3)){
        ind <- c(ind, which(df$CB == i & df$barcode %in% names(which(freq < 0.3))))
      }else{
        ind <- c(ind, which(df$CB == i))
      }
    }else if(len_bc == 2){
      freq <- table(temp$barcode)/nrow(temp)
      if(any(freq >= 0.5)){
        ind <- c(ind, which(df$CB == i & df$barcode %in% names(which(freq < 0.5))))
      }else{
        ind <- c(ind, which(df$CB == i))
      }
    }else{}
  }
  
  df <- df[-ind,]
  
  #remove CB duplicates
  df <- df[-which(duplicated(df$CB)),]
  
  #organize df
  df <- df[,c("CB", "barcode")]
  
  return(df)
}

