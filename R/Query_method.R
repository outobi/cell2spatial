#' Gene-wise Z-score Normalization
#'
#' Performs z-score normalization (standardization) for each gene across samples.
#' This transformation converts each gene's expression values to have mean=0 and sd=1
#' across all samples, following the formula: z = (x - mean)/sd
#'
#' @param expression_matrix A numeric matrix with genes as rows and samples as columns.
#'   Row names should be gene identifiers and column names should be sample identifiers.
#' @param handle_zero_variance Logical. If TRUE (default), genes with zero variance 
#'   (constant expression across all samples) will be set to 0. If FALSE, these genes
#'   will result in NaN values.
#' @param na_action Character string specifying how to handle missing values.
#'   Options: "fail" (default, stop if NA found), "omit" (exclude NAs from calculations),
#'   "zero" (replace NAs with 0 after normalization).
#'
#' @return A numeric matrix of the same dimensions as input, with each gene 
#'   normalized to have mean=0 and standard deviation=1 across samples.
#'
#' @details
#' The function applies z-score normalization gene-wise (row-wise), which is useful
#' for spatial transcriptomics data where you want to compare relative expression
#' patterns of each gene across different regions or conditions.
#' 
#' The transformation formula for each gene i: z_i = (x_i - μ_i) / σ_i
#' where μ_i is the mean expression and σ_i is the standard deviation for gene i.
#'
#' @examples
#' # Create example expression data
#' set.seed(123)
#' expression_data <- matrix(rnorm(100*20, mean=5, sd=2), 
#'                          nrow=100, ncol=20)
#' rownames(expression_data) <- paste0("Gene_", 1:100)
#' colnames(expression_data) <- paste0("Sample_", 1:20)
#' 
#' # Apply gene-wise normalization
#' normalized_data <- normalization_by_gene(expression_data)
#' 
#' # Check that each gene now has mean~0 and sd~1
#' head(apply(normalized_data, 1, mean))  # Should be close to 0
#' head(apply(normalized_data, 1, sd))    # Should be close to 1
#'
#' @export
normalization_by_gene <- function(expression_matrix, 
                                  handle_zero_variance = TRUE,
                                  na_action = "fail") {
  
  # Input validation
  if (!is.matrix(expression_matrix) && !is.data.frame(expression_matrix)) {
    stop("expression_matrix must be a matrix or data.frame")
  }
  
  if (!is.numeric(as.matrix(expression_matrix))) {
    stop("expression_matrix must contain numeric values")
  }
  
  # Convert to matrix if data.frame
  if (is.data.frame(expression_matrix)) {
    expression_matrix <- as.matrix(expression_matrix)
  }
  
  # Handle missing values based on na_action
  if (na_action == "fail" && any(is.na(expression_matrix))) {
    stop("Missing values (NA) found in expression_matrix. Set na_action='omit' or 'zero' to handle them.")
  }
  
  # Create output matrix
  normalized_matrix <- expression_matrix
  
  # Calculate mean and standard deviation for each gene (row)
  gene_means <- apply(expression_matrix, 1, function(x) {
    if (na_action == "omit") {
      mean(x, na.rm = TRUE)
    } else {
      mean(x)
    }
  })
  
  gene_sds <- apply(expression_matrix, 1, function(x) {
    if (na_action == "omit") {
      sd(x, na.rm = TRUE)
    } else {
      sd(x)
    }
  })
  
  # Apply z-score normalization for each gene
  for (i in seq_len(nrow(normalized_matrix))) {
    if (handle_zero_variance && (gene_sds[i] == 0 || is.na(gene_sds[i]))) {
      # Handle genes with zero variance (constant expression)
      normalized_matrix[i, ] <- 0
      if (!is.na(gene_sds[i]) && gene_sds[i] == 0) {
        warning(paste("Gene", rownames(expression_matrix)[i], 
                     "has zero variance. Set to 0 after normalization."))
      }
    } else {
      # Standard z-score normalization: (x - mean) / sd
      normalized_matrix[i, ] <- (expression_matrix[i, ] - gene_means[i]) / gene_sds[i]
    }
  }
  
  # Handle NAs after normalization if requested
  if (na_action == "zero") {
    normalized_matrix[is.na(normalized_matrix)] <- 0
  }
  
  # Preserve original row and column names
  rownames(normalized_matrix) <- rownames(expression_matrix)
  colnames(normalized_matrix) <- colnames(expression_matrix)
  
  return(normalized_matrix)
}


#' Apply Multiple Normalization Methods
#'
#' Convenience function to apply different normalization methods to expression data.
#'
#' @param expression_matrix A numeric matrix with genes as rows and samples as columns.
#' @param method Character string specifying normalization method. Options:
#'   "zscore" (default): Gene-wise z-score normalization
#'   "minmax": Min-max scaling (0-1 range) per gene
#'   "robust": Median and MAD-based robust scaling per gene
#'   "quantile": Quantile normalization across samples
#'
#' @return A normalized matrix of the same dimensions as input.
#'
#' @examples
#' # Create example data
#' set.seed(123)
#' expr_data <- matrix(rnorm(50*10), nrow=50, ncol=10)
#' rownames(expr_data) <- paste0("Gene_", 1:50)
#' colnames(expr_data) <- paste0("Sample_", 1:10)
#' 
#' # Different normalization methods
#' zscore_norm <- normalize_expression(expr_data, method = "zscore")
#' minmax_norm <- normalize_expression(expr_data, method = "minmax")
#' robust_norm <- normalize_expression(expr_data, method = "robust")
#'
#' @export
normalize_expression <- function(expression_matrix, method = "zscore") {
  
  valid_methods <- c("zscore", "minmax", "robust", "quantile")
  if (!method %in% valid_methods) {
    stop(paste("method must be one of:", paste(valid_methods, collapse = ", ")))
  }
  
  switch(method,
    "zscore" = normalization_by_gene(expression_matrix),
    
    "minmax" = {
      # Min-max scaling per gene: (x - min) / (max - min)
      normalized_matrix <- expression_matrix
      for (i in seq_len(nrow(expression_matrix))) {
        gene_min <- min(expression_matrix[i, ], na.rm = TRUE)
        gene_max <- max(expression_matrix[i, ], na.rm = TRUE)
        if (gene_max == gene_min) {
          normalized_matrix[i, ] <- 0
        } else {
          normalized_matrix[i, ] <- (expression_matrix[i, ] - gene_min) / (gene_max - gene_min)
        }
      }
      normalized_matrix
    },
    
    "robust" = {
      # Robust scaling using median and MAD
      normalized_matrix <- expression_matrix
      for (i in seq_len(nrow(expression_matrix))) {
        gene_median <- median(expression_matrix[i, ], na.rm = TRUE)
        gene_mad <- mad(expression_matrix[i, ], na.rm = TRUE)
        if (gene_mad == 0) {
          normalized_matrix[i, ] <- 0
        } else {
          normalized_matrix[i, ] <- (expression_matrix[i, ] - gene_median) / gene_mad
        }
      }
      normalized_matrix
    },
    
    "quantile" = {
      # Simple quantile normalization
      # Note: This is a basic implementation
      sorted_expr <- apply(expression_matrix, 2, sort, na.last = TRUE)
      mean_quantiles <- apply(sorted_expr, 1, mean, na.rm = TRUE)
      
      normalized_matrix <- expression_matrix
      for (j in seq_len(ncol(expression_matrix))) {
        ranks <- rank(expression_matrix[, j], na.last = "keep")
        normalized_matrix[, j] <- mean_quantiles[ranks]
      }
      normalized_matrix
    }
  )
}


#' Calculate Gene Signature Score for Single Cell Type
#'
#' Calculates the average expression of a gene signature in a specific cell type
#' and disease condition from single-cell RNA-seq data.
#'
#' @param expression_matrix A numeric matrix with genes as rows and cells as columns.
#'   Should be normalized/transformed expression data (e.g., output from normalization_by_gene).
#' @param cell_type Character string specifying the cell type to analyze.
#' @param disease_condition Character string specifying the disease condition 
#'   (e.g., "IPF", "Control", etc.).
#' @param metadata A data.frame containing single-cell metadata with at least two columns:
#'   - celltype: Cell type annotations for each cell
#'   - Diagnosis: Disease condition for each cell
#'   Rows should correspond to columns in expression_matrix.
#' @param celltype_column Character string specifying the column name in metadata 
#'   containing cell type information. Default: "celltype".
#' @param diagnosis_column Character string specifying the column name in metadata 
#'   containing diagnosis information. Default: "Diagnosis".
#'
#' @return A numeric value representing the average gene signature score for the 
#'   specified cell type and disease condition.
#'
#' @details
#' The function:
#' 1. Identifies cells matching the specified cell type and disease condition
#' 2. Extracts expression data for those cells
#' 3. Calculates the sum of all gene expression values
#' 4. Divides by the number of cells to get the average signature score
#' 
#' Formula: signature_score = sum(all_gene_expression) / number_of_cells
#'
#' @examples
#' # Create example data
#' set.seed(123)
#' expr_data <- matrix(rnorm(50*200), nrow=50, ncol=200)
#' rownames(expr_data) <- paste0("Gene_", 1:50)
#' colnames(expr_data) <- paste0("Cell_", 1:200)
#' 
#' # Create metadata (all samples from same condition)
#' metadata <- data.frame(
#'   celltype = rep(c("T_cell", "B_cell", "Macrophage"), length.out=200),
#'   Diagnosis = "IPF"  # All samples from IPF condition
#' )
#' 
#' # Calculate signature score
#' score <- calculate_signature_score(expr_data, "T_cell", "IPF", metadata)
#'
#' @export
calculate_signature_score <- function(expression_matrix, 
                                     cell_type, 
                                     disease_condition, 
                                     metadata,
                                     celltype_column = "celltype",
                                     diagnosis_column = "Diagnosis") {
  
  # Input validation
  if (!is.matrix(expression_matrix) && !is.data.frame(expression_matrix)) {
    stop("expression_matrix must be a matrix or data.frame")
  }
  
  if (!is.data.frame(metadata)) {
    stop("metadata must be a data.frame")
  }
  
  if (!celltype_column %in% colnames(metadata)) {
    stop(paste("Column", celltype_column, "not found in metadata"))
  }
  
  if (!diagnosis_column %in% colnames(metadata)) {
    stop(paste("Column", diagnosis_column, "not found in metadata"))
  }
  
  if (ncol(expression_matrix) != nrow(metadata)) {
    stop("Number of columns in expression_matrix must equal number of rows in metadata")
  }
  
  # Find cells matching the specified cell type and disease condition
  cell_indices <- which(metadata[[celltype_column]] == cell_type & 
                       metadata[[diagnosis_column]] == disease_condition)
  
  if (length(cell_indices) == 0) {
    warning(paste("No cells found for cell type:", cell_type, 
                 "and disease condition:", disease_condition))
    return(NA)
  }
  
  # Extract expression data for matching cells
  cell_expression_data <- expression_matrix[, cell_indices, drop = FALSE]
  
  # Calculate signature score: sum of all gene expression / number of cells
  signature_score <- sum(cell_expression_data) / length(cell_indices)
  
  return(signature_score)
}


#' Calculate Gene Signature Scores for Multiple Cell Types
#'
#' Calculates gene signature scores for multiple cell types in a specific disease condition.
#' This is equivalent to your original loop-based approach but with better error handling.
#'
#' @param expression_matrix A numeric matrix with genes as rows and cells as columns.
#' @param cell_types Character vector of cell type names to analyze.
#' @param disease_condition Character string specifying the disease condition.
#' @param metadata A data.frame containing single-cell metadata.
#' @param celltype_column Character string specifying the column name for cell types.
#' @param diagnosis_column Character string specifying the column name for diagnosis.
#'
#' @return A data.frame with columns:
#'   - Number: Sequential number (1 to length of cell_types)
#'   - Cell.type: Cell type names
#'   - signature.score: Calculated signature scores
#'
#' @examples
#' # Create example data
#' set.seed(123)
#' expr_data <- matrix(rnorm(50*300), nrow=50, ncol=300)
#' rownames(expr_data) <- paste0("Gene_", 1:50)
#' colnames(expr_data) <- paste0("Cell_", 1:300)
#' 
#' # Create metadata (all samples from same condition)
#' metadata <- data.frame(
#'   celltype = rep(c("T_cell", "B_cell", "Macrophage"), each=100),
#'   Diagnosis = "IPF"  # All samples from IPF condition
#' )
#' 
#' cell_types <- c("T_cell", "B_cell", "Macrophage")
#' 
#' # Calculate scores for all cell types
#' scores <- calculate_multiple_signature_scores(expr_data, cell_types, "IPF", metadata)
#'
#' @export
calculate_multiple_signature_scores <- function(expression_matrix, 
                                              cell_types, 
                                              disease_condition, 
                                              metadata,
                                              celltype_column = "celltype",
                                              diagnosis_column = "Diagnosis") {
  
  # Input validation
  if (!is.character(cell_types) || length(cell_types) == 0) {
    stop("cell_types must be a non-empty character vector")
  }
  
  # Initialize results vector
  signature_scores <- numeric(length(cell_types))
  
  # Calculate signature score for each cell type
  for (i in seq_along(cell_types)) {
    signature_scores[i] <- calculate_signature_score(
      expression_matrix = expression_matrix,
      cell_type = cell_types[i],
      disease_condition = disease_condition,
      metadata = metadata,
      celltype_column = celltype_column,
      diagnosis_column = diagnosis_column
    )
  }
  
  # Create results data.frame
  results <- data.frame(
    Number = seq_along(cell_types),
    Cell.type = cell_types,
    signature.score = signature_scores,
    stringsAsFactors = FALSE
  )
  
  return(results)
}


#' Query Gene Signature Across Conditions and Cell Types
#'
#' Advanced function to query gene signatures across multiple cell types and conditions.
#' Provides comprehensive analysis similar to your original approach.
#'
#' @param expression_matrix A numeric matrix with genes as rows and cells as columns.
#' @param cell_types Character vector of cell type names to analyze.
#' @param disease_conditions Character vector of disease conditions to compare.
#' @param metadata A data.frame containing single-cell metadata.
#' @param celltype_column Character string specifying the column name for cell types.
#' @param diagnosis_column Character string specifying the column name for diagnosis.
#' @param include_stats Logical. If TRUE, includes additional statistics. Default: FALSE.
#'
#' @return A data.frame with signature scores for each cell type and condition combination.
#'
#' @examples
#' # Create example data
#' set.seed(123)
#' expr_data <- matrix(rnorm(50*400), nrow=50, ncol=400)
#' rownames(expr_data) <- paste0("Gene_", 1:50)
#' colnames(expr_data) <- paste0("Cell_", 1:400)
#' 
#' # Create metadata
#' metadata <- data.frame(
#'   celltype = rep(c("T_cell", "B_cell"), each=200),
#'   Diagnosis = rep(c("IPF", "Control"), length.out=400)
#' )
#' 
#' cell_types <- c("T_cell", "B_cell")
#' conditions <- c("IPF", "Control")
#' 
#' # Comprehensive query
#' results <- query_signature_comprehensive(expr_data, cell_types, conditions, metadata)
#'
#' @export
query_signature_comprehensive <- function(expression_matrix, 
                                        cell_types, 
                                        disease_conditions, 
                                        metadata,
                                        celltype_column = "celltype",
                                        diagnosis_column = "Diagnosis",
                                        include_stats = FALSE) {
  
  # Initialize results list
  all_results <- list()
  
  # Calculate scores for each condition
  for (condition in disease_conditions) {
    condition_results <- calculate_multiple_signature_scores(
      expression_matrix = expression_matrix,
      cell_types = cell_types,
      disease_condition = condition,
      metadata = metadata,
      celltype_column = celltype_column,
      diagnosis_column = diagnosis_column
    )
    
    # Add condition information
    condition_results$Condition <- condition
    all_results[[condition]] <- condition_results
  }
  
  # Combine all results
  final_results <- do.call(rbind, all_results)
  rownames(final_results) <- NULL
  
  # Reorder columns
  final_results <- final_results[, c("Number", "Cell.type", "Condition", "signature.score")]
  
  # Add statistics if requested
  if (include_stats) {
    # Calculate additional metrics for each cell type-condition combination
    final_results$cell_count <- NA
    final_results$mean_expression <- NA
    
    for (i in seq_len(nrow(final_results))) {
      cell_indices <- which(metadata[[celltype_column]] == final_results$Cell.type[i] & 
                           metadata[[diagnosis_column]] == final_results$Condition[i])
      
      if (length(cell_indices) > 0) {
        final_results$cell_count[i] <- length(cell_indices)
        cell_data <- expression_matrix[, cell_indices, drop = FALSE]
        final_results$mean_expression[i] <- mean(cell_data)
      }
    }
  }
  
  return(final_results)
}