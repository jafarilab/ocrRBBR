#' @title Estimate Effective Sample Size (ESS) from Single-Cell RNA-seq Data
#'
#' @description This function estimates the effective sample size (ESS) of single-cell RNA-seq data
#' by accounting for correlation among cells within the same cell type.
#'
#' @param rnaseq_data A numeric matrix of RNA-seq expression values.
#'   Rows correspond to genes and columns correspond to cells.
#'   Expression values are assumed to be normalized (e.g., Seurat LogNormalize with scale.factor = 1e4).
#'
#' @param cell_type A data frame containing cell-type information.
#'   Row names must match the column names of `rnaseq_data`.
#'   Must include a column named `cell_type` specifying cell types.
#'
#' @return A numeric value representing the effective sample size (ESS),
#'   adjusted for within-cell-type correlation.
#'
#' @examples
#' \donttest{
#' # Load single-cell human dataset
#' data(multiome_human_mouse)  # loads atacseq_data, rnaseq_data, peaks_gr, cell_type, meta.data
#'
#' # Inspect loaded data
#' head(human_rnaseq_data)
#' head(human_cell_type)
#'
#' ess <- ESS(rnaseq_data = human_rnaseq_data, cell_type = human_cell_type)
#' print(ess)
#' }
#' @export
ESS <- function(rnaseq_data, cell_type) {

  # Cell name consistency ----
  if (!identical(colnames(rnaseq_data), rownames(cell_type))) {
    stop("ERROR: Column names of rnaseq_data must match row names of cell_type.")
  }

  # Annotate cell types ----
  cell_type <- factor(cell_type$cell_type)
  cell_types <- levels(cell_type)

  # Prepare output ----
  rho_per_type <- numeric(length(cell_types))
  names(rho_per_type) <- cell_types

  # Loop over cell types ----
  for (ct in cell_types) {
    message("Processing: ", ct)

    cells_ct <- which(cell_type == ct)

    if (length(cells_ct) < 10) {
      rho_per_type[ct] <- NA
      next
    }

    expr_ct <- rnaseq_data[, cells_ct, drop = FALSE]

    # Sample to reduce computation (optional)
    ncells <- ncol(expr_ct)
    if (ncells > 500) {
      # set.seed(1)
      sample_idx <- sample(seq_len(ncells), 500)
      expr_ct <- expr_ct[, sample_idx, drop = FALSE]
    }

    # Sample genes to limit size (optional)
    ngenes <- nrow(expr_ct)
    if (ngenes > 2000) {
      # set.seed(1)
      gene_idx <- sample(seq_len(ngenes), 2000)
      expr_ct <- expr_ct[gene_idx, , drop = FALSE]
    }

    # Convert to dense for correlation ----
    expr_ct_dense <- as.matrix(Matrix::t(expr_ct))  # cells as rows

    # # Compute correlation with suppressed warnings
    cor_ct <- tryCatch(
      suppressWarnings(cor(expr_ct_dense, use = "pairwise.complete.obs")),
      warning = function(w) NULL  # suppress the warning
    )

    # ---- Step 6: Compute correlation ----
    # cor_ct <- cor(expr_ct_dense, use = "pairwise.complete.obs")
    diag(cor_ct) <- NA

    # Estimate intra-cell-type correlation ----
    rho_ct <- mean(cor_ct, na.rm = TRUE)
    rho_per_type[ct] <- rho_ct
  }

  # Summarize ----
  rho_global <- mean(rho_per_type, na.rm = TRUE)
  N_cells <- ncol(rnaseq_data)
  ess_value <- round( N_cells / (1 + (N_cells - 1) * rho_global) )

  return(ess_value)
}
