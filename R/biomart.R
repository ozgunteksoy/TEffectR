get_ensembl_data <- function(attributes, filter_type, values, assembly, URL) {
  if (length(values) == 0) return(NULL)
  
  ensembl <- returnEnsembl(assembly, URL)
  if (is.null(ensembl)) return(NULL)
  
  data <- biomaRt::getBM(attributes = attributes, filters = filter_type, values = values, mart = ensembl, verbose = FALSE)
  if (nrow(data) == 0) return(NULL)

  names(data) <- c("geneID", "geneName", "chr", "start", "end", "strand")
  data <- data[c("chr", "start", "end", "strand", "geneID", "geneName")]

  data$strand[data$strand == "-1"] <- "-"
  data$strand[data$strand == "1"] <- "+"
  
  no_chr_prefix <- !stringr::str_detect(data$chr, "CHR")
  data$chr[no_chr_prefix] <- paste0("chr", data$chr[no_chr_prefix])

  return(data)
}

filterTranscriptID <- function(transcript_ids, assembly, URL) {
  get_ensembl_data(
    attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", 
                   "start_position", "end_position", "strand"),
    filter_type = "ensembl_transcript_id",
    values = transcript_ids,
    assembly = assembly,
    URL = URL
  )
}

filterTranscriptID_V <- function(transcript_ids, assembly, URL) {
  get_ensembl_data(
    attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", 
                   "start_position", "end_position", "strand"),
    filter_type = "ensembl_transcript_id_version",
    values = transcript_ids,
    assembly = assembly,
    URL = URL
  )
}

filterGeneName <- function(genes, assembly, URL) {
  if (length(genes) == 0) return(NULL)
  
  ensembl <- returnEnsembl(assembly, URL)
  if (is.null(ensembl)) return(NULL)
  
  gene_ids <- biomaRt::getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), 
                             filters = "hgnc_symbol", values = genes, mart = ensembl, verbose = FALSE)

  if (nrow(gene_ids) == 0) return(NULL)
  
  get_ensembl_data(
    attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", 
                   "start_position", "end_position", "strand"),
    filter_type = "ensembl_gene_id",
    values = gene_ids$ensembl_gene_id,
    assembly = assembly,
    URL = URL
  )
}

filterGeneID_V <- function(gene_ids, assembly, URL) {
  get_ensembl_data(
    attributes = c("ensembl_gene_id_version", "external_gene_name", "chromosome_name", 
                   "start_position", "end_position", "strand"),
    filter_type = "ensembl_gene_id_version",
    values = gene_ids,
    assembly = assembly,
    URL = URL
  )
}

filterGeneID <- function(gene_ids, assembly, URL) {
  get_ensembl_data(
    attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", 
                   "start_position", "end_position", "strand"),
    filter_type = "ensembl_gene_id",
    values = gene_ids,
    assembly = assembly,
    URL = URL
  )
}

returnEnsembl <- function(assembly, URL) {
  assembly_map <- list(
    hg19 = list(dataset = "hsapiens_gene_ensembl", GRCh = 37),
    Grch37 = list(dataset = "hsapiens_gene_ensembl", GRCh = 37),
    hg38 = list(dataset = "hsapiens_gene_ensembl", GRCh = 38),
    Grch38 = list(dataset = "hsapiens_gene_ensembl", GRCh = 38),
    mm9 = list(dataset = "mmusculus_gene_ensembl"),
    Grcm37 = list(dataset = "mmusculus_gene_ensembl"),
    mm10 = list(dataset = "mmusculus_gene_ensembl"),
    Grcm38 = list(dataset = "mmusculus_gene_ensembl")
  )

  assembly_info <- assembly_map[[assembly]]
  if (is.null(assembly_info)) {
    message("Not found assembly")
    return(NULL)
  }

  return(biomaRt::useEnsembl(
    biomart = "ensembl", 
    dataset = assembly_info$dataset, 
    GRCh = assembly_info$GRCh %||% NULL, 
    host = URL, 
    verbose = TRUE
  ))
}
