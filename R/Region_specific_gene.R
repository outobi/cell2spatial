#' Wilcoxon Comparison Test for Spatial Regions
#'
#' Performs Wilcoxon rank-sum tests for comparing gene expression between spatial
#' regions in transcriptomics data, specifically designed for region-specific analysis.
#'
#' @param expression_data A matrix or data frame of gene expression data with
#'   genes as rows and samples as columns. Row names should be gene symbols,
#'   column names should be sample IDs (S1, S2, S3, etc.)
#' @param sample_metadata A data frame with two columns:
#'   Column 1: Sample IDs (S1, S2, S3, etc.) matching column names in expression_data
#'   Column 2: Region types ("foci", "alveoli", "fibrosis", etc.)
#' @param target_region Character string specifying the target region type for comparison
#'   (e.g., "foci", "alveoli", "fibrosis")
#' @param sample_col Character string specifying the column name for sample IDs
#'   in sample_metadata (default: assumes first column)
#' @param region_col Character string specifying the column name for region types
#'   in sample_metadata (default: assumes second column)
#'
#' @return A data frame containing results for each gene:
#'   \item{gene_symbol}{Gene symbol from row names}
#'   \item{statistic}{The Wilcoxon test statistic}
#'   \item{p.value}{P-value from the test (or FDR-adjusted p-value if correction_method = "BH")}
#'   \item{method}{The method used for the test}
#'   \item{target_region}{Target region in comparison}
#'   \item{comparison_group}{Comparison group description}
#'   \item{n_target}{Sample size of target region}
#'   \item{n_other}{Sample size of other regions}
#'   \item{effect_size}{Rank-biserial correlation effect size}
#'   \item{median_target}{Median expression in target region}
#'   \item{median_other}{Median expression in other regions}
#'
#' @details
#' This function performs Wilcoxon rank-sum tests (Mann-Whitney U tests) to compare
#' gene expression levels between a target spatial region and all other regions combined.
#' It's designed for spatial transcriptomics analysis where you want to identify
#' region-specific gene expression patterns.
#'
#' The function processes each gene separately, comparing expression values from
#' samples in the target region versus samples from all other regions.
#'
#' @examples
#' \dontrun{
#' # Example with simulated data
#' set.seed(123)
#' # Expression data: 100 genes x 50 samples
#' expression_data <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
#' rownames(expression_data) <- paste0("Gene_", 1:100)
#' colnames(expression_data) <- paste0("S", 1:50)
#'
#' # Sample metadata
#' sample_metadata <- data.frame(
#'   sample_id = paste0("S", 1:50),
#'   region_type = rep(c("foci", "alveoli", "fibrosis"), c(20, 20, 10))
#' )
#'
#' # Compare foci vs other regions
#' result <- wilcox.comp(
#'   expression_data = expression_data,
#'   sample_metadata = sample_metadata,
#'   target_region = "foci"
#' )
#' head(result)
#' }
#'
#' @param correction_method Character string specifying multiple testing correction method.
#'   Options: "BH" for Benjamini-Hochberg FDR correction, "Nominal" for no correction (default: "Nominal")
#'
#' @importFrom stats wilcox.test p.adjust
#' @export
wilcox_compare_region <- function(expression_data,
                        sample_metadata,
                        target_region,
                        sample_col = NULL,
                        region_col = NULL,
                        correction_method = "Nominal") {

  # Input validation
  if (!is.matrix(expression_data) && !is.data.frame(expression_data)) {
    stop("expression_data must be a matrix or data frame")
  }

  if (!is.data.frame(sample_metadata)) {
    stop("sample_metadata must be a data frame")
  }
  
  if (!correction_method %in% c("BH", "Nominal")) {
    stop("correction_method must be either 'BH' or 'Nominal'")
  }

  if (ncol(sample_metadata) < 2) {
    stop("sample_metadata must have at least 2 columns")
  }

  # Set column names if not specified (assume first column is sample, second is region)
  if (is.null(sample_col)) {
    sample_col <- colnames(sample_metadata)[1]
  }
  if (is.null(region_col)) {
    region_col <- colnames(sample_metadata)[2]
  }

  if (!sample_col %in% colnames(sample_metadata)) {
    stop(paste("Column", sample_col, "not found in sample_metadata"))
  }

  if (!region_col %in% colnames(sample_metadata)) {
    stop(paste("Column", region_col, "not found in sample_metadata"))
  }

  if (!target_region %in% sample_metadata[[region_col]]) {
    stop(paste("Target region", target_region, "not found in", region_col))
  }

  # Get sample IDs for target and non-target regions
  target_samples <- sample_metadata[sample_metadata[[region_col]] == target_region, sample_col]
  other_samples <- sample_metadata[sample_metadata[[region_col]] != target_region, sample_col]

  # Match sample IDs with expression data column names
  target_indices <- match(target_samples, colnames(expression_data))
  other_indices <- match(other_samples, colnames(expression_data))

  # Remove NA matches
  target_indices <- target_indices[!is.na(target_indices)]
  other_indices <- other_indices[!is.na(other_indices)]

  if (length(target_indices) == 0) {
    stop("No target region samples found in expression data")
  }

  if (length(other_indices) == 0) {
    stop("No other region samples found in expression data")
  }

  cat("Target region samples:", length(target_indices), "\n")
  cat("Other region samples:", length(other_indices), "\n")

  # Initialize results vectors
  n_genes <- nrow(expression_data)
  gene_symbols <- if (!is.null(rownames(expression_data))) {
    rownames(expression_data)
  } else {
    paste0("Gene_", seq_len(n_genes))
  }

  results_list <- list()

  # Process each gene
  for (i in seq_len(n_genes)) {
    # Extract expression values for current gene
    target_expr <- as.numeric(expression_data[i, target_indices])
    other_expr <- as.numeric(expression_data[i, other_indices])

    # Remove NA values
    target_expr <- target_expr[!is.na(target_expr)]
    other_expr <- other_expr[!is.na(other_expr)]

    if (length(target_expr) >= 2 && length(other_expr) >= 2) {
      # Perform Wilcoxon test
      tryCatch({
        test_result <- stats::wilcox.test(
          x = target_expr,
          y = other_expr,
          alternative = "two.sided",
          exact = FALSE  # Use normal approximation
        )

        # Calculate effect size (rank-biserial correlation)
        effect_size <- calculate_rank_biserial_correlation(target_expr, other_expr)

        # Store results only if test was successful
        results_list[[i]] <- data.frame(
          gene_symbol = gene_symbols[i],
          statistic = as.numeric(test_result$statistic),
          p.value = test_result$p.value,
          method = test_result$method,
          target_region = target_region,
          comparison_group = "other_regions",
          n_target = length(target_expr),
          n_other = length(other_expr),
          effect_size = effect_size,
          median_target = median(target_expr),
          median_other = median(other_expr),
          stringsAsFactors = FALSE
        )

      }, error = function(e) {
        # Skip genes with test errors - do not add to results_list
        cat("Warning: Skipping", gene_symbols[i], "due to test error:", e$message, "\n")
        results_list[[i]] <<- NULL
      })
    } else {
      # Skip genes with insufficient data - do not add to results_list
      cat("Warning: Skipping", gene_symbols[i], "- insufficient data (target:", 
          length(target_expr), ", other:", length(other_expr), ")\n")
      results_list[[i]] <- NULL
    }
  }

  # Combine all results (only valid results, excluding NULLs)
  valid_results <- results_list[!sapply(results_list, is.null)]
  
  if (length(valid_results) == 0) {
    stop("No genes had sufficient data for analysis")
  }
  
  final_results <- do.call(rbind, valid_results)
  
  # Apply multiple testing correction if requested
  if (correction_method == "BH") {
    final_results$p.value <- stats::p.adjust(final_results$p.value, method = "BH")
    cat("Applied Benjamini-Hochberg FDR correction\n")
  }
  
  # Sort results by p-value (most significant first)
  final_results <- final_results[order(final_results$p.value), ]
  
  # Reset row names
  rownames(final_results) <- NULL
  
  cat("Successfully analyzed", nrow(final_results), "out of", n_genes, "genes\n")
  if (correction_method == "BH") {
    cat("Significant genes (FDR < 0.05):", sum(final_results$p.value < 0.05, na.rm = TRUE), "\n")
  } else {
    cat("Significant genes (p < 0.05):", sum(final_results$p.value < 0.05, na.rm = TRUE), "\n")
  }

  return(final_results)
}

#' Calculate rank-biserial correlation as effect size for Wilcoxon test
#'
#' @param x Numeric vector for group 1
#' @param y Numeric vector for group 2
#' @return Numeric value representing the rank-biserial correlation
#' @keywords internal
calculate_rank_biserial_correlation <- function(x, y) {
  n1 <- length(x)
  n2 <- length(y)

  # Calculate U statistic
  u1 <- sum(outer(x, y, ">")) + 0.5 * sum(outer(x, y, "=="))
  u2 <- n1 * n2 - u1

  # Rank-biserial correlation
  r <- (u1 - u2) / (n1 * n2)

  return(r)
}


#' T-Test Comparison for Spatial Regions
#'
#' Performs two-sample t-tests for comparing gene expression between spatial
#' regions in transcriptomics data, specifically designed for region-specific analysis.
#' Uses Cohen's d for effect size calculation.
#'
#' @param expression_data A matrix or data frame of gene expression data with
#'   genes as rows and samples as columns. Row names should be gene symbols,
#'   column names should be sample IDs (S1, S2, S3, etc.)
#' @param sample_metadata A data frame with two columns:
#'   Column 1: Sample IDs (S1, S2, S3, etc.) matching column names in expression_data
#'   Column 2: Region types ("foci", "alveoli", "fibrosis", etc.)
#' @param target_region Character string specifying the target region type for comparison
#'   (e.g., "foci", "alveoli", "fibrosis")
#' @param sample_col Character string specifying the column name for sample IDs
#'   in sample_metadata (default: assumes first column)
#' @param region_col Character string specifying the column name for region types
#'   in sample_metadata (default: assumes second column)
#' @param var_equal Logical. If TRUE, assumes equal variances (Student's t-test).
#'   If FALSE (default), uses Welch's t-test which does not assume equal variances.
#'
#' @return A data frame containing results for each gene:
#'   \item{gene_symbol}{Gene symbol from row names}
#'   \item{statistic}{The t-test statistic}
#'   \item{p.value}{P-value from the test (or FDR-adjusted p-value if correction_method = "BH")}
#'   \item{method}{The method used for the test}
#'   \item{target_region}{Target region in comparison}
#'   \item{comparison_group}{Comparison group description}
#'   \item{n_target}{Sample size of target region}
#'   \item{n_other}{Sample size of other regions}
#'   \item{cohens_d}{Cohen's d effect size}
#'   \item{mean_target}{Mean expression in target region}
#'   \item{mean_other}{Mean expression in other regions}
#'   \item{sd_target}{Standard deviation in target region}
#'   \item{sd_other}{Standard deviation in other regions}
#'
#' @details
#' This function performs two-sample t-tests to compare gene expression levels
#' between a target spatial region and all other regions combined. By default,
#' it uses Welch's t-test which does not assume equal variances.
#'
#' Cohen's d is calculated as the standardized mean difference:
#' d = (mean_target - mean_other) / pooled_sd
#'
#' For Welch's t-test, the pooled standard deviation is calculated as:
#' sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2) / (n1 + n2 - 2))
#'
#' Effect size interpretation (Cohen's guidelines):
#' - |d| < 0.2: negligible
#' - 0.2 <= |d| < 0.5: small
#' - 0.5 <= |d| < 0.8: medium
#' - |d| >= 0.8: large
#'
#' @examples
#' \dontrun{
#' # Example with simulated data
#' set.seed(123)
#' # Expression data: 100 genes x 50 samples
#' expression_data <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
#' rownames(expression_data) <- paste0("Gene_", 1:100)
#' colnames(expression_data) <- paste0("S", 1:50)
#'
#' # Sample metadata
#' sample_metadata <- data.frame(
#'   sample_id = paste0("S", 1:50),
#'   region_type = rep(c("foci", "alveoli", "fibrosis"), c(20, 20, 10))
#' )
#'
#' # Compare foci vs other regions
#' result <- t_test_compare_region(
#'   expression_data = expression_data,
#'   sample_metadata = sample_metadata,
#'   target_region = "foci"
#' )
#' head(result)
#' }
#'
#' @param correction_method Character string specifying multiple testing correction method.
#'   Options: "BH" for Benjamini-Hochberg FDR correction, "Nominal" for no correction (default: "Nominal")
#'
#' @importFrom stats t.test p.adjust
#' @export
t_test_compare_region <- function(expression_data,
                                   sample_metadata,
                                   target_region,
                                   sample_col = NULL,
                                   region_col = NULL,
                                   var_equal = FALSE,
                                   correction_method = "Nominal") {

  # Input validation
  if (!is.matrix(expression_data) && !is.data.frame(expression_data)) {
    stop("expression_data must be a matrix or data frame")
  }

  if (!is.data.frame(sample_metadata)) {
    stop("sample_metadata must be a data frame")
  }

  if (ncol(sample_metadata) < 2) {
    stop("sample_metadata must have at least 2 columns")
  }
  
  if (!correction_method %in% c("BH", "Nominal")) {
    stop("correction_method must be either 'BH' or 'Nominal'")
  }

  # Set column names if not specified (assume first column is sample, second is region)
  if (is.null(sample_col)) {
    sample_col <- colnames(sample_metadata)[1]
  }
  if (is.null(region_col)) {
    region_col <- colnames(sample_metadata)[2]
  }

  if (!sample_col %in% colnames(sample_metadata)) {
    stop(paste("Column", sample_col, "not found in sample_metadata"))
  }

  if (!region_col %in% colnames(sample_metadata)) {
    stop(paste("Column", region_col, "not found in sample_metadata"))
  }

  if (!target_region %in% sample_metadata[[region_col]]) {
    stop(paste("Target region", target_region, "not found in", region_col))
  }

  # Get sample IDs for target and non-target regions
  target_samples <- sample_metadata[sample_metadata[[region_col]] == target_region, sample_col]
  other_samples <- sample_metadata[sample_metadata[[region_col]] != target_region, sample_col]

  # Match sample IDs with expression data column names
  target_indices <- match(target_samples, colnames(expression_data))
  other_indices <- match(other_samples, colnames(expression_data))

  # Remove NA matches
  target_indices <- target_indices[!is.na(target_indices)]
  other_indices <- other_indices[!is.na(other_indices)]

  if (length(target_indices) == 0) {
    stop("No target region samples found in expression data")
  }

  if (length(other_indices) == 0) {
    stop("No other region samples found in expression data")
  }

  cat("Target region samples:", length(target_indices), "\n")
  cat("Other region samples:", length(other_indices), "\n")

  # Initialize results vectors
  n_genes <- nrow(expression_data)
  gene_symbols <- if (!is.null(rownames(expression_data))) {
    rownames(expression_data)
  } else {
    paste0("Gene_", seq_len(n_genes))
  }

  results_list <- list()

  # Process each gene
  for (i in seq_len(n_genes)) {
    # Extract expression values for current gene
    target_expr <- as.numeric(expression_data[i, target_indices])
    other_expr <- as.numeric(expression_data[i, other_indices])

    # Remove NA values
    target_expr <- target_expr[!is.na(target_expr)]
    other_expr <- other_expr[!is.na(other_expr)]

    if (length(target_expr) >= 2 && length(other_expr) >= 2) {
      # Perform t-test
      tryCatch({
        test_result <- stats::t.test(
          x = target_expr,
          y = other_expr,
          alternative = "two.sided",
          var.equal = var_equal
        )

        # Calculate Cohen's d effect size
        cohens_d <- calculate_cohens_d(target_expr, other_expr)

        # Store results only if test was successful
        results_list[[i]] <- data.frame(
          gene_symbol = gene_symbols[i],
          statistic = as.numeric(test_result$statistic),
          p.value = test_result$p.value,
          method = test_result$method,
          target_region = target_region,
          comparison_group = "other_regions",
          n_target = length(target_expr),
          n_other = length(other_expr),
          cohens_d = cohens_d,
          mean_target = mean(target_expr),
          mean_other = mean(other_expr),
          sd_target = sd(target_expr),
          sd_other = sd(other_expr),
          stringsAsFactors = FALSE
        )

      }, error = function(e) {
        # Skip genes with test errors - do not add to results_list
        cat("Warning: Skipping", gene_symbols[i], "due to test error:", e$message, "\n")
        results_list[[i]] <<- NULL
      })
    } else {
      # Skip genes with insufficient data - do not add to results_list
      cat("Warning: Skipping", gene_symbols[i], "- insufficient data (target:", 
          length(target_expr), ", other:", length(other_expr), ")\n")
      results_list[[i]] <- NULL
    }
  }

  # Combine all results (only valid results, excluding NULLs)
  valid_results <- results_list[!sapply(results_list, is.null)]
  
  if (length(valid_results) == 0) {
    stop("No genes had sufficient data for analysis")
  }
  
  final_results <- do.call(rbind, valid_results)
  
  # Apply multiple testing correction if requested
  if (correction_method == "BH") {
    final_results$p.value <- stats::p.adjust(final_results$p.value, method = "BH")
    cat("Applied Benjamini-Hochberg FDR correction\n")
  }
  
  # Sort results by p-value (most significant first)
  final_results <- final_results[order(final_results$p.value), ]
  
  # Reset row names
  rownames(final_results) <- NULL
  
  cat("Successfully analyzed", nrow(final_results), "out of", n_genes, "genes\n")
  if (correction_method == "BH") {
    cat("Significant genes (FDR < 0.05):", sum(final_results$p.value < 0.05, na.rm = TRUE), "\n")
  } else {
    cat("Significant genes (p < 0.05):", sum(final_results$p.value < 0.05, na.rm = TRUE), "\n")
  }

  return(final_results)
}


#' Compare Gene Expression Between Regions Using t-test
#'
#' @param x Numeric vector for group 1
#' @param y Numeric vector for group 2
#' @return Numeric value representing Cohen's d
#' @keywords internal
calculate_cohens_d <- function(x, y) {
  n1 <- length(x)
  n2 <- length(y)
  
  mean_x <- mean(x)
  mean_y <- mean(y)
  sd_x <- sd(x)
  sd_y <- sd(y)
  
  # Pooled standard deviation
  pooled_sd <- sqrt(((n1 - 1) * sd_x^2 + (n2 - 1) * sd_y^2) / (n1 + n2 - 2))
  
  # Cohen's d
  d <- (mean_x - mean_y) / pooled_sd
  
  return(d)
}