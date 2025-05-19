library(data.table)
library(stringr)
library(stringdist)
library(Biostrings)

#convert fastq quality sequence to scores 
convert_quality_scores <- function(quality_string, offset = 33){
  quality_scores <- as.integer(charToRaw(quality_string)) - offset
  return(mean(quality_scores))
}

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

extract_barcode <- function(sam, len = 15, known_seq = "CAGATCTTAGCCACTTTTTA"){
  #loading sam
  df <- fread(file = sam, skip = 3)
  
  #match the known sequence on the 3 prime end (用‘[]’直接索引计算速度更快,但是用subset语法更加直观，容易理解和维护)
  df <- subset(df, agrepl(pattern = known_seq, x = df$V10, max.distance = allowed_mismatch(string = known_seq)))
  
  #find the end position of barcodes
  match_positions <- unlist(lapply(df$V10, function(seq){
    matches <- matchPattern(pattern = known_seq, subject = DNAString(seq), max.mismatch = 1, with.indels = T)
    start(matches)
  }))
  
  #remove short barcode
  valid_positions <- match_positions > len
  df <- df[valid_positions,]
  match_positions <- match_positions[valid_positions]
  
  #truncate barcodes to specified length and calculate quality scores
  df$barcode <- substr(x = df$V10, start = match_positions - len, stop = match_positions - 1)
  df$barcode_quality_score <- unlist(lapply(substr(x = df$V11, start = match_positions - len, stop = match_positions - 1), convert_quality_scores, offset = 33))
  
  #remove low quality reads
  df <- subset(df, barcode_quality_score >= 20)
  
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
  return(as.vector(df$barcode))
}

