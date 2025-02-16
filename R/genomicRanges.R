getUpstream <- function(df, length, isWithGeneBody) {
  dfnew <- df

  indexP <- df$strand == "+"
  indexN <- df$strand == "-"

  dfnew$start[indexP] <- df$start[indexP] - length
  dfnew$end[indexN] <- df$end[indexN] + length

  if (!isWithGeneBody) {
    dfnew$end[indexP] <- df$start[indexP]
    dfnew$start[indexN] <- df$end[indexN]
  }

  return(dfnew)
}

getDownstream <- function(df, length, isWithGeneBody) {
  dfnew <- df

  indexP <- df$strand == "+"
  indexN <- df$strand == "-"

  dfnew$start[indexN] <- df$start[indexN] - length
  dfnew$end[indexP] <- df$end[indexP] + length

  if (!isWithGeneBody) {
    dfnew$end[indexN] <- df$start[indexN]
    dfnew$start[indexP] <- df$end[indexP]
  }

  return(dfnew)
}

getDownAndUpStream <- function(df, len_up, len_down) {
  dfnew <- df

  indexP <- df$strand == "+"
  indexN <- df$strand == "-"

  dfnew$start[indexP] <- df$start[indexP] - len_up
  dfnew$end[indexP] <- df$end[indexP] + len_down

  dfnew$start[indexN] <- df$start[indexN] - len_down
  dfnew$end[indexN] <- df$end[indexN] + len_up

  return(dfnew)
}

makeGRangeObj <- function(df) {
  gr <- GenomicRanges::GRanges(df$chr, IRanges(df$start, df$end), strand = df$strand)
  mcols(gr) <- df[, 5:ncol(df)]
  mcols(gr) <- cbind(mcols(gr), df[, 2:3])
  return(gr)
}

makeGrObj_Unstrand <- function(df) {
  gr <- GenomicRanges::GRanges(df$chr, IRanges(df$start, df$end), strand = "*")
  mcols(gr) <- df[, 5:ncol(df)]
  mcols(gr) <- cbind(mcols(gr), df[, 2:3])
  return(gr)
}

toFindOverlaps <- function(gr_repeats, gr_genome, strand) {
  names(mcols(gr_genome))[3:4] <- c("gStart", "gEnd")
  names(mcols(gr_repeats))[4:5] <- c("rStart", "rEnd")

  ignore_strand <- strand != "same"
  m <- GenomicRanges::findOverlaps(gr_genome, gr_repeats, ignore.strand = ignore_strand)

  gr_genome.matched <- gr_genome[queryHits(m)]
  mcols(gr_genome.matched) <- cbind(mcols(gr_genome.matched), mcols(gr_repeats[subjectHits(m)]))

  return(gr_genome.matched)
}
