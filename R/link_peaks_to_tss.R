#' @title Identifies ATAC-seq peaks that fall within a user-defined genomic window surrounding each gene.
#'
#' @description This function identifies ATAC-seq peaks that are located within a specified window around
#' the transcription start sites (TSS) of genes, and links those peaks to the respective genes.
#'
#' @param gtf_file A character string specifying the path to a GTF file containing gene annotations.
#' @param peaks_gr A \code{GRanges} object representing the ATAC-seq peaks with their genomic locations and associated metadata.
#' @param gene_list A character vector of gene names to filter for. Only peaks within the window around TSS of these genes will be considered. Default is \code{NA}, in which case all genes in the GTF will be used.
#' @param tss_window An integer specifying the window size around the TSS. Default is \code{100000} (±100kb).
#'
#' @return A data frame containing the following columns:
#' \item{peak}{The genomic coordinates of the peaks.}
#' \item{peak_id}{The unique ID of each peak.}
#' \item{gene_id}{The gene ID associated with the peak.}
#' \item{gene_name}{The gene name associated with the peak.}
#' \item{transcript_id}{The transcript ID associated with the peak.}
#' \item{gene_type}{The type of the gene (e.g., protein-coding).}
#' \item{distance}{The distance from the peak center to the TSS.}
#'
#' @examples
#' # Load bulk mouse dataset
#' data(multiome_human_mouse)  # This will load atacseq_data, rnaseq_data, peaks_gr
#'
#' # Inspect loaded data
#' head(mouse_atacseq_data)
#' head(mouse_rnaseq_data)
#' head(mouse_peaks_gr)
#'
#'
#' # Path to the GTF file in the package
#' gtf_file <- system.file("extdata", "gencode.vM25.annotation.sample.gtf", package = "ocrRBBR")
#'
#' # Example usage for linking peaks to TSS
#' linked_peaks <- link_peaks_to_tss(
#'   gtf_file = gtf_file,
#'   peaks_gr = mouse_peaks_gr,
#'   gene_list = c("Rag2", "Spi1"),
#'   tss_window = 100000
#' )
#'
#' # Filter results for a specific gene
#' linked_peaks_gene <- linked_peaks[linked_peaks$gene_name == "Rag2", ]
#' print(linked_peaks_gene)
#'
#' @export
link_peaks_to_tss <- function(gtf_file, peaks_gr, gene_list = NA, tss_window = NA){

  # window threshold default
  GenomeInfoDb::seqlevels(peaks_gr)
  if(is.na(tss_window)){
    tss_window <- 100000
  }

  # Import GTF
  gtf_data <- rtracklayer::import(gtf_file)

  # Filter transcript entries
  tx_gtf <- gtf_data[gtf_data$type == "transcript"]

  # Extract transcript information
  tx_info <- data.frame(
    seqnames      = as.character(GenomicRanges::seqnames(tx_gtf)),
    start         = GenomicRanges::start(tx_gtf),
    end           = GenomicRanges::end(tx_gtf),
    strand        = as.character(GenomicRanges::strand(tx_gtf)),
    gene_id       = gsub("\\.\\d+$", "", tx_gtf$gene_id),
    transcript_id = gsub("\\.\\d+$", "", tx_gtf$transcript_id),
    gene_name     = tx_gtf$gene_name,
    gene_type     = tx_gtf$gene_type,
    stringsAsFactors = FALSE
  )

  # Convert back to GRanges
  tx_gr <- GenomicRanges::makeGRangesFromDataFrame(
    tx_info,
    seqnames.field = "seqnames",
    start.field    = "start",
    end.field      = "end",
    strand.field   = "strand",
    keep.extra.columns = TRUE
  )

  # Extract TSS for each transcript
  tss_gr <- GenomicRanges::resize(tx_gr, width = 1, fix = "start")

  # Filter by gene list (optional)
  if (length(gene_list) > 0 && all(!is.na(gene_list))) {
    tss_gr <- tss_gr[tss_gr$gene_name %in% gene_list]
  }

  # Define ± window around TSS
  tss_window_gr <- GenomicRanges::promoters(
    tss_gr,
    upstream   = tss_window,
    downstream = tss_window
  )

  # Find peaks overlapping the TSS windows
  overlaps <- suppressWarnings(
    GenomicRanges::findOverlaps(peaks_gr, tss_window_gr)
  )

  # Extract overlapping entries
  peak_hits <- peaks_gr[S4Vectors::queryHits(overlaps)]
  tss_hits  <- tss_gr[S4Vectors::subjectHits(overlaps)]

  # Compute distance between peak center and TSS
  peak_center <- GenomicRanges::start(peak_hits) +
    (GenomicRanges::width(peak_hits) / 2)

  tss_position <- GenomicRanges::start(tss_hits)

  distance <- abs(peak_center - tss_position)

  # Build metadata table
  df <- data.frame(
    peak          = paste0(
      as.character(GenomicRanges::seqnames(peak_hits)),
      ":",
      GenomicRanges::start(peak_hits),
      "-",
      GenomicRanges::end(peak_hits)
    ),
    peak_id       = S4Vectors::mcols(peak_hits)$peakID,
    gene_id       = tss_hits$gene_id,
    gene_name     = tss_hits$gene_name,
    transcript_id = tss_hits$transcript_id,
    gene_type     = tss_hits$gene_type,
    distance      = distance,
    stringsAsFactors = FALSE
  )

  # Keep unique peak–gene pairs with minimum distance
  linked_peaks <- df |>
    dplyr::group_by(peak, gene_id) |>
    dplyr::summarise(
      gene_name     = dplyr::first(gene_name),
      gene_type     = dplyr::first(gene_type),
      transcript_id = dplyr::first(transcript_id),
      peak_id       = dplyr::first(peak_id),
      min_distance  = min(distance),
      .groups       = "drop"
    )


  return(linked_peaks)
}
