#' Compare Gene Expression Between Cell Types Using t-test
#'
#' Performs comprehensive cell type-specific differential expression analysis
#' including t-tests and Cohen's d effect size calculations for spatial transcriptomics data.
#'
#' @param normalized_count A matrix or data frame of normalized gene expression counts
#'   with genes as rows and samples as columns
#' @param meta_data A data frame containing sample metadata with columns for
#'   cell type and diagnosis information
#' @param target_cell_type Character string specifying the target cell type for analysis
#' @param diagnosis Character string specifying the diagnosis condition (e.g., "IPF")
#' @param celltype_col Character string specifying the column name for cell type
#'   in meta_data (default: "celltype")
#' @param diagnosis_col Character string specifying the column name for diagnosis
#'   in meta_data (default: "Diagnosis")
#' @param min_samples Minimum number of samples required for analysis (default: 3)
#' @param correction_method Character string specifying multiple testing correction method.
#'   Options: "BH" for Benjamini-Hochberg FDR correction, "Nominal" for no correction (default: "Nominal")
#'
#' @return A data frame containing:
#'   \item{index}{Gene index/row number}
#'   \item{gene_id}{Gene identifier (if available)}
#'   \item{d}{Cohen's d effect size}
#'   \item{stat.sig}{P-value from t-test (or FDR-adjusted p-value if correction_method = "BH")}
#'   \item{mean_target}{Mean expression in target cell type}
#'   \item{mean_non_target}{Mean expression in non-target cells}
#'   \item{n_target}{Number of target cell type samples}
#'   \item{n_non_target}{Number of non-target samples}
#'
#' @details
#' This function performs differential expression analysis comparing a target cell type
#' against all other cell types within a specific diagnosis condition. It calculates
#' Cohen's d effect sizes and performs Welch's t-tests for each gene.
#'
#' The analysis follows the methodology used in spatial transcriptomics studies
#' where target cell types are compared against reference populations to identify
#' cell type-specific gene expression patterns.
#'
#' @examples
#' \dontrun{
#' # Example with simulated data
#' set.seed(123)
#' # Create simulated expression data
#' expression_data <- matrix(rnorm(1000 * 50), nrow = 1000, ncol = 50)
#' rownames(expression_data) <- paste0("Gene_", 1:1000)
#' colnames(expression_data) <- paste0("Sample_", 1:50)
#'
#' # Create metadata
#' metadata <- data.frame(
#'   celltype = rep(c("Fibroblasts", "Epithelial", "Immune"), c(20, 20, 10)),
#'   Diagnosis = rep("IPF", 50)
#' )
#'
#' # Perform cell type comparison
#' result <- t_test_compare_cell(
#'   normalized_count = expression_data,
#'   meta_data = metadata,
#'   target_cell_type = "Fibroblasts",
#'   diagnosis = "IPF"
#' )
#' head(result)
#' }
#'
#' @importFrom stats t.test
#' @export
t_test_compare_cell <- function(normalized_count,
                              meta_data,
                              target_cell_type,
                              diagnosis,
                              celltype_col = "celltype",
                              diagnosis_col = "Diagnosis",
                              min_samples = 3,
                              correction_method = "Nominal") {

  # Input validation
  if (!is.matrix(normalized_count) && !is.data.frame(normalized_count)) {
    stop("normalized_count must be a matrix or data frame")
  }
  
  if (!correction_method %in% c("BH", "Nominal")) {
    stop("correction_method must be either 'BH' or 'Nominal'")
  }

  if (!is.data.frame(meta_data)) {
    stop("meta_data must be a data frame")
  }

  if (!celltype_col %in% colnames(meta_data)) {
    stop(paste("Column", celltype_col, "not found in meta_data"))
  }

  if (!diagnosis_col %in% colnames(meta_data)) {
    stop(paste("Column", diagnosis_col, "not found in meta_data"))
  }

  # Filter metadata for specific diagnosis
  meta_filtered <- meta_data[meta_data[[diagnosis_col]] == diagnosis, ]

  if (nrow(meta_filtered) == 0) {
    stop(paste("No samples found for diagnosis:", diagnosis))
  }

  # Get sample indices for target vs non-target cell types
  target_indices <- which(meta_filtered[[celltype_col]] == target_cell_type)
  non_target_indices <- which(meta_filtered[[celltype_col]] != target_cell_type)

  if (length(target_indices) < min_samples) {
    stop(paste("Insufficient target cell type samples. Found:", length(target_indices),
               "Required:", min_samples))
  }

  if (length(non_target_indices) < min_samples) {
    stop(paste("Insufficient non-target samples. Found:", length(non_target_indices),
               "Required:", min_samples))
  }

  # Get column names/indices for expression data
  target_samples <- rownames(meta_filtered)[target_indices]
  non_target_samples <- rownames(meta_filtered)[non_target_indices]

  # Match sample names between expression data and metadata
  expr_target_indices <- match(target_samples, colnames(normalized_count))
  expr_non_target_indices <- match(non_target_samples, colnames(normalized_count))

  # Remove NA matches
  expr_target_indices <- expr_target_indices[!is.na(expr_target_indices)]
  expr_non_target_indices <- expr_non_target_indices[!is.na(expr_non_target_indices)]

  if (length(expr_target_indices) == 0 || length(expr_non_target_indices) == 0) {
    stop("No matching samples found between expression data and metadata")
  }

  # Initialize results vectors
  n_genes <- nrow(normalized_count)
  stat_sig <- numeric(n_genes)
  d_values <- numeric(n_genes)
  mean_target <- numeric(n_genes)
  mean_non_target <- numeric(n_genes)

  # Perform analysis for each gene
  for (i in seq_len(n_genes)) {
    # Extract expression values for current gene
    target_expr <- as.numeric(normalized_count[i, expr_target_indices])
    non_target_expr <- as.numeric(normalized_count[i, expr_non_target_indices])

    # Remove NA values
    target_expr <- target_expr[!is.na(target_expr)]
    non_target_expr <- non_target_expr[!is.na(non_target_expr)]

    if (length(target_expr) >= 2 && length(non_target_expr) >= 2) {
      # Perform Welch's t-test (does not assume equal variances)
      tryCatch({
        t_test_result <- stats::t.test(target_expr, non_target_expr, var.equal = FALSE)
        stat_sig[i] <- t_test_result$p.value

        # Calculate Cohen's d
        pooled_sd <- sqrt(((length(target_expr) - 1) * var(target_expr) +
                          (length(non_target_expr) - 1) * var(non_target_expr)) /
                         (length(target_expr) + length(non_target_expr) - 2))

        if (pooled_sd > 0) {
          d_values[i] <- (mean(target_expr) - mean(non_target_expr)) / pooled_sd
        } else {
          d_values[i] <- 0
        }

        # Store means
        mean_target[i] <- mean(target_expr)
        mean_non_target[i] <- mean(non_target_expr)

      }, error = function(e) {
        stat_sig[i] <- NA
        d_values[i] <- NA
        mean_target[i] <- NA
        mean_non_target[i] <- NA
      })
    } else {
      stat_sig[i] <- NA
      d_values[i] <- NA
      mean_target[i] <- NA
      mean_non_target[i] <- NA
    }
  }

  # Create results data frame
  gene_ids <- if (!is.null(rownames(normalized_count))) {
    rownames(normalized_count)
  } else {
    paste0("Gene_", seq_len(n_genes))
  }

  results <- data.frame(
    index = seq_len(n_genes),
    gene_id = gene_ids,
    d = d_values,
    stat.sig = stat_sig,
    mean_target = mean_target,
    mean_non_target = mean_non_target,
    n_target = length(expr_target_indices),
    n_non_target = length(expr_non_target_indices),
    target_cell_type = target_cell_type,
    diagnosis = diagnosis,
    stringsAsFactors = FALSE
  )

  # Remove rows with NA values
  results <- results[!is.na(results$stat.sig), ]
  
  # Apply multiple testing correction if requested
  if (correction_method == "BH") {
    results$stat.sig <- stats::p.adjust(results$stat.sig, method = "BH")
    cat("Applied Benjamini-Hochberg FDR correction\n")
  }
  
  # Sort results by p-value (most significant first)
  results <- results[order(results$stat.sig), ]
  
  # Reset row names
  rownames(results) <- NULL
  
  cat("Successfully analyzed", nrow(results), "genes\n")
  if (correction_method == "BH") {
    cat("Significant genes (FDR < 0.05):", sum(results$stat.sig < 0.05, na.rm = TRUE), "\n")
  } else {
    cat("Significant genes (p < 0.05):", sum(results$stat.sig < 0.05, na.rm = TRUE), "\n")
  }

  return(results)
}