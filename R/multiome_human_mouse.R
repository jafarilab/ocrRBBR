#' multiome ATAC-seq and RNA-seq datasets from mouse and human
#'
#' Paired chromatin accessibility (ATAC-seq) and gene expression (RNA-seq)
#' datasets across diverse mouse and human hematopoietic and immune cell types.
#' These data enable integrative regulatory analyses linking cis-regulatory regions
#' (open chromatin regions, OCRs) to gene expression variability.
#'
#' @details
#' The dataset consists of multiple objects for mouse and human:
#'
#' \describe{
#'   \item{\code{mouse_atacseq_data}}{
#'     A numeric matrix of chromatin accessibility values for mouse.
#'     Rows correspond to ATAC-seq peaks (cis-regulatory regions) and columns
#'     correspond to 85 mouse cell types. Values represent quantile-normalized accessibility
#'     signals for each peak in each cell type.
#'   }
#'
#'   \item{\code{mouse_rnaseq_data}}{
#'     A numeric matrix of gene expression values for mouse.
#'     Rows correspond to genes and columns correspond to the same set of mouse
#'     cell types as in \code{mouse_atacseq_data}. Values represent quantile-normalized RNA-seq
#'     expression levels.
#'   }
#'
#'   \item{\code{mouse_peaks_gr}}{
#'     A \code{\link[GenomicRanges]{GRanges}} object containing genomic
#'     coordinates and annotations for all ATAC-seq peaks in mouse, including chromosome,
#'     genomic ranges, peak identifiers, evolutionary conservation scores
#'     (phastCons), and peak significance metrics.
#'   }
#'
#'   \item{\code{human_atacseq_data}}{
#'     A numeric matrix of chromatin accessibility values for human.
#'     Rows correspond to ATAC-seq peaks and columns to cells or samples.
#'     Values are normalized by the \emph{ReadsInTSS} metric per cell, where raw peak counts
#'     are divided by the number of Tn5 insertions in transcription start site (TSS) regions
#'     for each cell to correct for sequencing depth and accessibility differences.
#'   }
#'
#'   \item{\code{human_rnaseq_data}}{
#'     A numeric matrix of gene expression values for human.
#'     Rows correspond to genes and columns to cells or samples.
#'     Values are normalized using Seurat's \code{LogNormalize} method with a scale factor of 10,000:
#'     \code{NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 1e4)}.
#'   }
#'
#'   \item{\code{human_peaks_gr}}{
#'     A \code{\link[GenomicRanges]{GRanges}} object containing genomic
#'     coordinates and annotations for all ATAC-seq peaks in human, including chromosome,
#'     genomic ranges, and peak identifiers.
#'   }
#'   
#'   \item{\code{human_meta_data}}{
#'     A matrix that contains additional information per cell for human samples.
#'     This matrix must include the following columns:
#'     \describe{
#'       \item{nCount_RNA}{The total number of RNA molecules (unique molecular identifiers, UMIs) detected per cell, calculated as the sum of UMI counts across all genes.}
#'       \item{nFeature_RNA}{The number of genes detected per cell (genes with at least one UMI).}
#'       \item{Mitochondrial percentage}{The percentage of reads that map to mitochondrial genes, which can be used to assess the quality of the sample.}
#'     }
#'   }
#' }
#'
#' @format
#' A list-like dataset containing:
#' \describe{
#'   \item{\code{mouse_atacseq_data}}{numeric matrix (peaks x cell types)}
#'   \item{\code{mouse_rnaseq_data}}{numeric matrix (genes x cell types)}
#'   \item{\code{mouse_peaks_gr}}{\code{GRanges} object with peak annotations}
#'   \item{\code{human_atacseq_data}}{numeric matrix (peaks x cells)}
#'   \item{\code{human_rnaseq_data}}{numeric matrix (genes x cells)}
#'   \item{\code{human_peaks_gr}}{\code{GRanges} object with peak annotations}
#'   \item{\code{human_meta_data}}{\code{matrix} of metadata information, including RNA counts, feature counts, and mitochondrial percentage per cell}
#' }
#'
#' @source
#' Mouse data: GSE100738.
#' Human data: Derived from published human single-cell multiome ATAC-seq and RNA-seq experiments
#' across hematopoietic and immune cell populations.
#'
#' @usage
#' data(multiome_human_mouse)
#'
#' @examples
#' data(multiome_human_mouse)
#'
#' # Inspect mouse chromatin accessibility
#' head(mouse_atacseq_data)
#'
#' # Inspect mouse gene expression
#' head(mouse_rnaseq_data)
#'
#' # Mouse genomic coordinates of peaks
#' mouse_peaks_gr
#'
#' # Inspect human chromatin accessibility
#' head(human_atacseq_data)
#'
#' # Inspect human gene expression
#' head(human_rnaseq_data)
#'
#' # Human genomic coordinates of peaks
#' human_peaks_gr
#'
#' # Inspect human metadata
#' head(human_meta_data)