
source("R/biomart.R")
source("R/genomicRanges.R")


# the function that to filter biomart and return genome information.
get_intervals <- function(x, assembly, ID.type, URL) {
  if (length(x) == 0 || !is.character(ID.type)) {
    return(NULL)
  }
  
  func_map <- list(
    ensembl_gene_name = filterGeneName,
    ensembl_transcript_id = filterTranscriptID,
    ensembl_transcript_id_version = filterTranscriptID_V,
    ensembl_gene_id_version = filterGeneID_V,
    ensembl_gene_id = filterGeneID
  )
  
  func <- func_map[[ID.type]]
  if (!is.null(func)) {
    return(func(x, assembly, URL))
  }
  
  return(NULL)
}


# the function is called as b is returned overlap positions between genes and repeats
get_overlaps <- function(g, r, strand, distance, repeat_type) {
  options(warn = -1)

  # Giriş doğrulama
  if (!is.data.frame(g) || !is.data.frame(r) || !is.numeric(distance)) {
    stop("Invalid input: 'g' and 'r' must be data frames, 'distance' must be numeric.")
  }
  
  # repeat_type filtresi
  if (!is.null(repeat_type)) {
    r <- r[r$repeat_type == repeat_type, ]
  } else {
    stop("Please choose one of the repeat options.")
  }

  # Mesafeye göre upstream/downstream işlemi
  if (distance > 0) {
    g <- getUpstream(g, distance, FALSE)
  } else if (distance < 0) {
    g <- getDownstreams(g, distance, FALSE)
  }

  # GRange nesnesine dönüştürme
  g <- makeGRangeObj(g)
  r <- makeGRangeObj(r)

  # Örtüşmeleri bulma
  toFindOverlaps(r, g, strand)
}


rm_format <- function(filepath) {
  dt <- biomartr::read_rm(filepath)
  last <- as.data.frame(stringr::str_split_fixed(dt$matching_class, "/", 2))
  
  dt <- data.frame(
    chr = dt$qry_id,
    start = dt$qry_start,
    end = dt$qry_end,
    strand = dt$matching_repeat,
    repeat_name = dt$repeat_id,
    repeat_type = last$V1,
    repeat_family = last$V2
  )
  
  dt$strand <- ifelse(dt$strand == "C", "-", as.character(dt$strand))
  return(dt)
}

count_repeats <- function(bamlist, namelist, ranges) {
  bamfiles <- paste(bamlist, collapse = " ")
  bamFile <- Rsamtools::BamFile(bamlist[1])
  
  if (!stringr::str_detect(GenomeInfoDb::seqnames(GenomeInfoDb::seqinfo(bamFile)), "chr")) {
    seqlevels(ranges) <- gsub("chr", " ", GenomeInfoDb::seqlevels(ranges))
  }

  df <- as.data.frame(ranges)[c(1, 13, 14, 6, 4, 5, 2, 3, 7, 10, 11, 12)]
  df <- as.data.frame(apply(df, 2, function(x) gsub("\\s+", "", x))) 
  df$repeat_family <- sub("^$", ".", df$repeat_family)
  
  write.table(df, "overlapped.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  system(paste("bedtools multicov -s -f 1 -D -bams", bamfiles, "-bed overlapped.bed > counts.txt"))
  
  counts <- read.csv("counts.txt", sep = "\t", header = FALSE)
  colnames(counts) <- c(colnames(df), namelist)
  return(counts)
}

summarize_repeat_counts <- function(counts, namelist) {
  if (is.null(counts) || is.null(namelist)) return(NULL)

  col_indexes <- which(colnames(counts) %in% namelist)
  aggregate(counts[, col_indexes], 
            by = list(geneName = counts$geneName, repeatClass = counts$repeat_type, repeatName = counts$repeat_name), 
            FUN = sum)
}

apply_lm <- function(gene.annotation, gene.counts, repeat.counts, covariates = NULL, prefix) {
  df1 <- gene.annotation[, 5:6]
  y <- data.frame(geneID = row.names(gene.counts), gene.counts)
  df <- merge(df1, y, by = "geneID")[, -1]
  
  e1 <- repeat.counts[, 4:ncol(repeat.counts)]
  e1$geneName <- seq.int(nrow(e1))
  count.matrix <- rbind(df, e1)
  
  dge <- edgeR::DGEList(count.matrix[, -1])
  keep <- edgeR::filterByExpr(dge, min.total.count = 10)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)
  v <- limma::voom(dge)
  
  vall <- count.matrix[row.names(v$E), ]
  new_voom <- cbind(data.frame(geneName = vall$geneName), v$E)
  
  v_ids_for_repeats <- setdiff(vall[, 1], df$geneName)
  col_indexes <- which(vall[, 1] %in% v_ids_for_repeats)
  temp <- repeat.counts[v_ids_for_repeats, ]
  temp$geneName <- paste(temp$geneName, temp$repeatClass, temp$repeatName, sep = ":")
  
  prepLMdata_repeats <- temp
  prepLMdata_repeats[, 4:ncol(prepLMdata_repeats)] <- new_voom[col_indexes, -1]

  v_ids_for_genes <- match(temp$geneName, new_voom$geneName)
  prepLMdata_genes <- unique(na.omit(new_voom[v_ids_for_genes, ]))

  writingResultOfVoom(vall, prefix)

  lm_list <- list()
  for (r in unique(prepLMdata_repeats$geneName)) {
    hit1 <- prepLMdata_genes$geneName == r
    s1 <- prepLMdata_genes[hit1, -1]
    
    hit2 <- prepLMdata_repeats$geneName == r
    s2 <- prepLMdata_repeats[hit2, ]
    
    if (nrow(s1) == 1 && nrow(s2) > 0) {
      output <- data.frame(matrix(0, nrow = 1, ncol = ncol(s2) - 3))
      colnames(output) <- c(unique(s2$geneName), s2$repeatName)
      
      output[1:length(s1), 1] <- as.numeric(s1)
      for (var in 1:(ncol(output) - 1)) {
        output[1:length(s2[var, 4:ncol(s2)]), var + 1] <- as.numeric(s2[var, 4:ncol(s2)])
      }
      
      if (!is.null(covariates) && nrow(covariates) == ncol(count.matrix) - 1) {
        output <- cbind(output, covariates)
      }
      
      formula <- paste(colnames(output)[1], "~ .", sep = "")
      lm.out <- lm(formula = formula, data = output)
      lm_list <- append(lm_list, list(lm.out))
    }
  }
  
  writingResultOfLM(lm_list, ncol(covariates), prefix)
  return(lm_list)
}

writingResultOfLM <- function(lm_list, ncov, prefix) {
  results <- data.frame(matrix(ncol = 6, nrow = length(lm_list)))
  colnames(results) <- c("GeneName", "RepeatName", "r.squared", "adjusted-r.squared", "model-p.value", "individual-p.vals")

  for (i in seq_along(lm_list)) {
    lm_model <- lm_list[[i]]
    results$GeneName[i] <- colnames(lm_model$model)[1]
    results$RepeatName[i] <- paste(colnames(lm_model$model)[2:(ncol(lm_model$model) - ncov)], collapse = " ")
    results$r.squared[i] <- summary(lm_model)$r.squared
    results$adjusted.r.squared[i] <- summary(lm_model)$adj.r.squared
    results$model.p.value[i] <- lmp(lm_model)

    p_vals <- summary(lm_model)$coefficients[, 4]
    non_na_vals <- names(p_vals[-1])
    missing_vals <- setdiff(names(lm_model$coefficients[-1]), non_na_vals)
    
    results$individual.p.vals[i] <- paste(c(
      paste(non_na_vals, p_vals[-1], sep = " : "),
      paste(missing_vals, ": NA", sep = "")
    ), collapse = " // ")
  }

  write.table(results, paste0(prefix, "-lm-results.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}

lmp <- function(model) {
  if (!inherits(model, "lm")) stop("Not an 'lm' object")
  f <- summary(model)$fstatistic
  p <- pf(f[1], f[2], f[3], lower.tail = FALSE)
  return(as.numeric(p))
}

writingResultOfVoom <- function(v, prefix) {
  write.table(v, paste0(prefix, "-cpm-values.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}
