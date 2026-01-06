#' @title Predicts OCR-driven Boolean rules for a gene
#'
#' @description This function predicts Boolean rule sets for a given gene using single-cell multi-omics datasets,
#' including RNA-seq gene expression and ATAC-seq peak signals in the gene's flanking regions,
#' across samples.
#'
#' @param rnaseq_data A numeric matrix of RNA-seq expression values. Rows correspond to genes, columns
#'   correspond to cells or samples. RNA-seq values are assumed to be normalized using Seurat's
#'   \code{LogNormalize} method with a scale factor of 10,000:
#'   \code{NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 1e4)}.
#' @param atacseq_data A numeric matrix of ATAC-seq signal intensities. Rows correspond to peaks,
#'   columns correspond to cells or samples. ATAC-seq counts are assumed to be normalized
#'   by the \emph{ReadsInTSS} metric on a per-cell basis, where raw peak counts are divided by
#'   the number of Tn5 insertions falling within transcription start site (TSS) regions for each cell.
#'   This normalization corrects for differences in sequencing depth and chromatin accessibility
#'   signal across cells.
#'
#'   ReadsInTSS values are typically obtained from ArchR and applied as column-wise scaling
#'   factors to the ATAC-seq count matrix prior to downstream analysis.
#'
#' @param gene_name A character string specifying the gene for which to infer Boolean rules.
#' @param peak_ids A vector of peak identifiers corresponding to rows in \code{atacseq_data} to be used as candidate regulatory regions for \code{gene_name}.
#' @param max_feature An integer specifying the maximum number of input features allowed in a Boolean rule. Default is 3.
#' @param slope The slope parameter for the sigmoid activation function. Default is 10.
#' @param num_cores Number of parallel workers to use for computation. Adjust according to your system. Default is NA.
#' @param ESS Effective sample size of the single-cell data after accounting for noise and cell-to-cell correlation.
#' @param meta_data A matrix that contains additional information per cell, such as:
#' - **nCount_RNA**: The total number of RNA molecules (unique molecular identifiers, UMIs) detected per cell, calculated as the sum of UMI counts across all genes.
#' - **nFeature_RNA**: The number of genes detected per cell (genes with at least one UMI).
#' - **Mitochondrial percentage**: The percentage of reads that map to mitochondrial genes, which can be used to assess the quality of the sample.
#'
#' This information is typically stored as columns in the `meta_data` object, which is associated with each cell in the dataset.
#'
#' @return A list containing predicted Boolean rules and associated metrics for the input gene.
#'
#' @examples
#' # Load single-cell human dataset
#' data(multiome_human_mouse)  # loads atacseq_data, rnaseq_data, peaks_gr, cell_type, meta_data
#'
#' # Inspect loaded data
#' head(human_atacseq_data)
#' head(human_rnaseq_data)
#' head(human_peaks_gr)
#' head(human_meta_data)
#'
#' # Example usage:
#' peak_ids <- c(83456, 83458, 83460, 83475, 83482)
#'
#' boolean_rules <- ocrRBBR_single_cell(human_rnaseq_data, human_atacseq_data,
#'        "CD74", peak_ids = peak_ids, max_feature = 3, slope = 6,
#'         num_cores = 1, ESS = 261, meta_data = human_meta_data)
#'
#' print(boolean_rules)
#' @export
ocrRBBR_single_cell <- function(rnaseq_data, atacseq_data, gene_name, peak_ids, max_feature = NA, slope = NA, num_cores = NA, ESS = NA, meta_data){
  message("Starting processing for gene: ", gene_name, " ...")

  ## ---- Check 1: gene_name exists in rnaseq_data ----
  if (!(gene_name %in% rownames(rnaseq_data))) {
    stop("ERROR: gene_name '", gene_name,
         "' does not exist in rownames(rnaseq_data).")
  }

  ## ---- Check 2: peak_ids exist in atacseq_data ----
  # Identify the missing peaks
  missing_peaks <- peak_ids[!peak_ids %in% rownames(atacseq_data)]

  if (length(missing_peaks) == length(peak_ids)) {
    stop("ERROR: None of the provided peak_ids exist in atacseq_data.")
  }

  # If there are any missing peaks, show a warning message
  if (length(missing_peaks) > 0) {
    warning("WARNING: The following peak_ids do not exist in atacseq_data: ",
            paste(missing_peaks, collapse = ", "))
  }

  # Check which peak_ids exist in atacseq_data
  peak_ids <- peak_ids[peak_ids %in% rownames(atacseq_data)]

  ## ---- Check 3: Column names match between RNA-seq & ATAC-seq ----
  if (!identical(colnames(rnaseq_data), colnames(atacseq_data))) {
    stop("ERROR: Column names of rnaseq_data and atacseq_data do NOT match.\n",
         "This indicates differences in cell identities or ordering.")
  }

  ## ---- Check 4: cell names match between RNA-seq & meta_data ----
  if (!identical(colnames(rnaseq_data), rownames(meta_data))) {
    stop("ERROR: Column names of rnaseq_data and meta_data do NOT match.\n",
         "This indicates differences in cell identities or ordering between RNA-seq data and metadata.")
  }

  ## If all checks passed:
  message("All input checks passed.")
  # -------------------------------
  # 1. Prepare RNA data (log10 -> scaled)
  # -------------------------------
  rnaseq_vec <- rnaseq_data[gene_name, ]
  rnaseq_vec <- rnaseq_vec - min(rnaseq_vec)  # shift to zero
  rnaseq_vec <- rnaseq_vec / quantile(rnaseq_vec, 0.975, na.rm = TRUE)  # scale

  # -------------------------------
  # 2. Prepare ATAC data (log10 -> z-score -> sigmoid)
  # -------------------------------
  atacseq_mat <- atacseq_data[as.character(peak_ids), , drop = FALSE]

  # Check for rows with zero variance and remove them
  var_check <- apply(atacseq_mat, 1, var)  # Calculate variance for each row (peak)

  # Identify rows with zero variance
  nonzero_var_rows <- which(var_check > 0)

  # Remove rows with zero variance
  atacseq_mat <- atacseq_mat[nonzero_var_rows, ]

  # If no rows remain, stop the analysis
  if (nrow(atacseq_mat) == 0) {
    stop("ERROR: All selected peaks have zero variance. Cannot proceed with analysis.")
  }

  atacseq_mat <- Matrix::t(atacseq_mat)

  # sigmoid-transform z-scores
  atacseq_mat <- apply(atacseq_mat, 2, function(x){
    z <- (x - mean(x)) / sd(x)
    1 / (1 + exp(-z))
  })

  atacseq_mat <- as.data.frame(atacseq_mat)

  # -------------------------------
  # 3. Combine ATAC + RNA
  # -------------------------------
  data <- cbind(atacseq_mat, RNA = as.numeric(rnaseq_vec))

  data[data >= 1] <- 0.9999
  data[data <= 0] <- 0.0001
  ##############################################################################
  # sigmoid-transform z-scores
  meta_data <- apply(meta_data, 2, function(x){
    z <- (x - mean(x)) / sd(x)
    1 / (1 + exp(-z))
  })
  ##############################################################################
  weight_threshold <- 0

  if (is.na(max_feature)) {
    max_feature <- 3
  }
  max_feature <- min(max_feature, ncol(data)-1)

  if (is.na(slope)) {
    slope <- 10
  }

  if (is.na(num_cores)) {
    num_cores <- parallel::detectCores()
  }

  if (is.na(ESS)){
    ESS <- nrow(data)
  }
  cat("training process started with ", num_cores, " computing cores\n")

  progress_percent <- 0
  pb <- txtProgressBar(min = 0, max = 100, style = 3, width = 20)
  sigmoid <- function(x) {
    1/(1 + exp(-slope * (x - 0.5)))
  }
  ToComputeLogic <- function(X, ORD) {
    if (is.vector(X)) {
      SAMP <- X
    }
    else {
      SAMP <- X[, c(ORD)]
    }
    ncol_SAMP <- length(c(ORD))
    if (length(c(ORD)) > 1) {
      SAMP_PARTITIONED <- matrix(rep(1, nrow(SAMP) * (2^ncol_SAMP)), nrow(SAMP), 2^ncol_SAMP)
    }
    else {
      SAMP <- matrix(SAMP, ncol = ncol_SAMP)
      SAMP_PARTITIONED <- matrix(rep(1, length(SAMP) * (2^ncol_SAMP)), length(SAMP), 2^ncol_SAMP)
    }
    new_logic_index <- 1
    for (gene_id in 1:ncol_SAMP) {
      SAMP_PARTITIONED[, new_logic_index] <- SAMP_PARTITIONED[, new_logic_index] * SAMP[, gene_id]
    }
    for (k in 1:ncol_SAMP) {
      D = combn(1:ncol_SAMP, k)
      for (i in 1:ncol(D)) {
        index_of_not <- t(D[, i])
        new_logic_index <- new_logic_index + 1
        for (j in 1:ncol(index_of_not)) {
          SAMP[, index_of_not[j]] <- (1 - SAMP[, index_of_not[j]])
        }
        for (gene_id in 1:ncol_SAMP) {
          SAMP_PARTITIONED[, new_logic_index] <- SAMP_PARTITIONED[, new_logic_index] * SAMP[, gene_id]
        }
        if (is.vector(X)) {
          SAMP <- X
        }
        else {
          SAMP <- X[, c(ORD)]
        }
      }
    }
    return(SAMP_PARTITIONED)
  }

  predictive_data <- data[, -ncol(data)]
  n <- ncol(data) - 1
  total_iterations <- 0
  for (i in 1:max_feature) {
    total_iterations <- total_iterations + ncol(combn(n, i))
  }
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  current_it <- 0
  LOGIC_VALUES <- list(`1` = list(), `2` = list())
  for (k in 1:max_feature) {
    C <- combn(1:n, k)
    current_it <- current_it + ncol(C)
    results <- foreach::foreach(j = 1:ncol(C), .combine = rbind) %dopar%
      {
        ORD <- t(C[, j])
        x <- cbind(ToComputeLogic(predictive_data, ORD), meta_data)
        y <- data[, ncol(data)]
        model <- glmnet::glmnet(x, y, alpha = 0)
        cv_model <- glmnet::cv.glmnet(x, y, alpha = 0)
        best_lambda <- cv_model$lambda.min
        best_model <- glmnet::glmnet(x, y, alpha = 0, lambda = best_lambda)
        y_predicted <- predict(model, s = best_lambda,
                               newx = x)
        y_predicted <- sigmoid(y_predicted)
        sst <- sum((y - mean(y))^2)
        sse <- sum((y_predicted - y)^2)
        rsq <- 1 - sse/sst
        p <- (2^k) + 1
        rsq_adj <- 1 - ((1 - rsq) * (nrow(data) - 1)/(nrow(data) - p - 1))
        BIC <- ESS * log(sse/ESS) + p * log(ESS)
        list(AGRE_OUT = c(ORD, coef(best_model)[1:length(coef(best_model))], rsq, rsq_adj, BIC), AGRE_OUT_pred = c(), R2 = rsq_adj)
      }
    progress_percent <- round(100 * current_it/total_iterations, 2)
    utils::setTxtProgressBar(pb, progress_percent)
    LOGIC_VALUES[[as.character(k)]] <- results
  }
  W2SYMBOL_ORIGINAL <- function(logic_significance, predicted_impact_set) {
    k <- length(predicted_impact_set)
    reg_char <- predicted_impact_set
    x <- list()
    logical_sets <- array(x, c(2^k, 2))
    logical_sets_index <- 1
    logical_sets[[1, 1]] <- 0
    logical_sets[[1, 2]] <- paste0("AND(", paste(reg_char, collapse = ","), ")")
    for (i in 1:k) {
      not_comb <- combn(1:k, i)
      for (j in 1:ncol(not_comb)) {
        reg_char_tmp <- reg_char
        logical_sets_index <- logical_sets_index + 1
        logical_sets[[logical_sets_index, 1]] <- not_comb[, j]
        for (l in 1:(length(not_comb[, j]))) {
          reg_char_tmp[not_comb[, j][l]] <- paste0("\u00AC", reg_char[not_comb[, j][l]])
        }
        logical_sets[[logical_sets_index, 2]] <- paste0("AND(", paste0(reg_char_tmp, collapse = ","), ")")
      }
    }
    LOGIC_VECTOR <- c()
    for (i in 1:(2^k)) {
      if (logic_significance[i] == 1) {
        LOGIC_VECTOR <- c(LOGIC_VECTOR, logical_sets[[i, 2]])
      }
    }
    LOGIC_VECTOR <- paste0("[", paste0(LOGIC_VECTOR, collapse = ","), "]")
    return(LOGIC_VECTOR)
  }
  W2SYMBOL <- function(logic_significance, predicted_impact_set, R) {
    if ((R == 1) || (R > 2)) {
      predicted_impact_set <- colnames(data)[c(predicted_impact_set)]
    }
    if (length(predicted_impact_set) >= 2) {
      LOGIC_VECTOR <- W2SYMBOL_ORIGINAL(logic_significance, predicted_impact_set)
    }
    else if (length(predicted_impact_set) == 1) {
      if (all(logic_significance == c(0, 0))) {
        LOGIC_VECTOR <- 0
      }
      else if (all(logic_significance == c(1, 0))) {
        LOGIC_VECTOR <- paste0("[", predicted_impact_set, "]")
      }
      else if (all(logic_significance == c(0, 1))) {
        LOGIC_VECTOR <- paste0("[", "\u00AC", predicted_impact_set, "]")
      }
      else if (all(logic_significance == c(1, 1))) {
        LOGIC_VECTOR <- 1
      }
    }
    return(LOGIC_VECTOR)
  }
  gate_info <- c()
  for (k in 1:max_feature) {
    if (is.vector(LOGIC_VALUES[[as.character(k)]])) {
      in_nodes <- LOGIC_VALUES[[as.character(k)]][[1]][1:k]
      coef <- LOGIC_VALUES[[as.character(k)]][[1]][(k + 2):(k + 2 + 2^k - 1)]
      logic_significance <- as.integer(coef > weight_threshold)
      LOGIC_VECTOR <- W2SYMBOL(logic_significance, in_nodes, R = 1)
      bic_val <- LOGIC_VALUES[[as.character(k)]][[1]][length(LOGIC_VALUES[[as.character(k)]][[1]])]
      r2_adj <- LOGIC_VALUES[[as.character(k)]][[1]][length(LOGIC_VALUES[[as.character(k)]][[1]]) - 1]
      Rule_num <- sum(logic_significance)
      Rule_coef <- c(paste(round(coef, 2), collapse = ":"), " ", " ")
      gate_info <- rbind(gate_info, c(LOGIC_VECTOR, r2_adj, bic_val, k, i, paste(colnames(data)[in_nodes], collapse = "."), Rule_num, Rule_coef))
    }
    else {
      results <- foreach(i = 1:nrow(LOGIC_VALUES[[as.character(k)]]),
                         .combine = rbind) %dopar% {
                           in_nodes <- LOGIC_VALUES[[as.character(k)]][[i, 1]][1:k]
                           coef <- LOGIC_VALUES[[as.character(k)]][[i, 1]][(k + 2):(k + 2 + 2^k - 1)]
                           logic_significance <- as.integer(coef > weight_threshold)
                           LOGIC_VECTOR <- W2SYMBOL(logic_significance, in_nodes, R = 1)
                           bic_val <- LOGIC_VALUES[[as.character(k)]][[i, 1]][length(LOGIC_VALUES[[as.character(k)]][[i, 1]])]
                           r2_adj <- LOGIC_VALUES[[as.character(k)]][[i, 1]][length(LOGIC_VALUES[[as.character(k)]][[i, 1]]) - 1]
                           Rule_num <- sum(logic_significance)
                           Rule_coef <- c(paste(round(coef, 2), collapse = ":"), " ", " ")
                           c(LOGIC_VECTOR, r2_adj, bic_val, k, i, paste(colnames(data)[in_nodes], collapse = "."), Rule_num, Rule_coef)
                         }
      gate_info <- rbind(gate_info, results)
    }
  }

  gate_info <- as.data.frame(gate_info)
  colnames(gate_info) <- c("Boolean_Rule", "R2", "BIC", "Input_Size",
                           "Index", "Features", "Active_Conjunctions", "Rule_Coefficients",
                           "Weights Layer1, Sub-Rule2", "Weights Layer2")

  gate_info <- gate_info[, 1:(ncol(gate_info) - 2)]
  gate_info <- gate_info[, c("Boolean_Rule", "R2", "BIC", "Rule_Coefficients")]

  gate_info$BIC <- round(as.numeric(gate_info$BIC), 2)
  gate_info$R2 <- round(as.numeric(gate_info$R2), 2)
  boolean_rules <- gate_info[order(gate_info$BIC), ]
  rownames(boolean_rules) <- NULL
  parallel::stopCluster(cl)
  foreach::registerDoSEQ()
  close(pb)
  return(boolean_rules)
}
