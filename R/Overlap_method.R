#' Calculate Region-Specific Gene Overlap with Cell-Type Gene Signatures
#'
#' Performs hypergeometric test to assess enrichment or depletion of region-specific
#' genes within cell-type-specific gene signatures from single-cell RNA-seq data.
#'
#' @param region_signature_genes Character vector of gene symbols from region-specific
#'   signature (e.g., genes from fibroblast foci, alveoli, etc.).
#' @param cell_signature_list A list where each element contains a character vector
#'   of gene symbols for a specific cell type. Names of the list elements should be
#'   cell type names.
#' @param universe_size Integer specifying the total number of genes in the universe
#'   (background). Default: 24470 (typical for human transcriptome).
#' @param output_file Character string specifying the output file path for Excel export.
#'   If NULL (default), no file is exported.
#'
#' @return A data frame containing overlap analysis results:
#'   \item{cell.type}{Cell type names from the input list}
#'   \item{overlap_count}{Number of overlapping genes between region signature and cell signature}
#'   \item{expected_overlap}{Expected number of overlapping genes by chance}
#'   \item{direction}{Either "enriched" or "depleted" based on comparison with expectation}
#'   \item{p.value}{P-value from hypergeometric test}
#'   \item{-log10(p.value)}{Negative log10 transformed p-value (always positive)}
#'   \item{p.value.with.direction}{Signed -log10(p.value): positive for enrichment, negative for depletion}
#'
#' @details
#' This function performs hypergeometric enrichment testing to determine whether
#' region-specific genes are significantly over-represented (enriched) or
#' under-represented (depleted) in cell-type-specific gene signatures.
#'
#' The hypergeometric test is performed for each cell type signature:
#' - If overlap > expected: tests for enrichment (upper tail test)
#' - If overlap < expected: tests for depletion (lower tail test)
#'
#' The hypergeometric probability is calculated as:
#' P(X = k) where X ~ Hypergeometric(N, K, n)
#' - N = universe_size (total genes)
#' - K = size of cell signature
#' - n = size of region signature
#' - k = observed overlap
#'
#' @examples
#' \dontrun{
#' # Example with fibroblast foci signature
#' ff_genes <- c("COL1A1", "COL1A2", "ACTA2", "FN1", "POSTN")
#' 
#' # Cell type signatures (list of gene vectors)
#' cell_signatures <- list(
#'   Fibroblasts = c("COL1A1", "COL3A1", "FN1", "DCN", "LUM"),
#'   Macrophages = c("CD68", "CD163", "MSR1", "MRC1", "MARCO"),
#'   T_cells = c("CD3D", "CD3E", "CD8A", "CD4", "IL7R")
#' )
#' 
#' # Calculate overlap
#' result <- calculate_signature_overlap(
#'   region_signature_genes = ff_genes,
#'   cell_signature_list = cell_signatures,
#'   universe_size,
#'   output_file = "ff_overlap_score.xlsx"
#' )
#' 
#' # View results
#' print(result)
#' }
#'
#' @importFrom stats phyper
#' @export
calculate_signature_overlap <- function(region_signature_genes,
                                       cell_signature_list,
                                       universe_size,
                                       output_file = NULL) {
  
  # Input validation
  if (!is.character(region_signature_genes) || length(region_signature_genes) == 0) {
    stop("region_signature_genes must be a non-empty character vector")
  }
  
  if (!is.list(cell_signature_list) || length(cell_signature_list) == 0) {
    stop("cell_signature_list must be a non-empty list")
  }
  
  if (!is.numeric(universe_size) || universe_size <= 0) {
    stop("universe_size must be a positive number")
  }
  
  # Get cell type names
  cell_types <- names(cell_signature_list)
  if (is.null(cell_types) || any(cell_types == "")) {
    cell_types <- paste0("CellType_", seq_along(cell_signature_list))
    warning("cell_signature_list elements are not named. Using default names: ", 
            paste(cell_types, collapse = ", "))
  }
  
  n_cell_types <- length(cell_signature_list)
  
  # Initialize result vectors
  overlap_count <- numeric(n_cell_types)
  expected_overlap <- numeric(n_cell_types)
  direction <- character(n_cell_types)
  p_value <- numeric(n_cell_types)
  
  # Size of region signature
  k <- length(region_signature_genes)
  
  cat("Region signature size:", k, "genes\n")
  cat("Universe size:", universe_size, "genes\n")
  cat("Testing overlap with", n_cell_types, "cell type signatures...\n\n")
  
  # Calculate overlap for each cell type
  for (i in seq_along(cell_signature_list)) {
    cell_sig <- cell_signature_list[[i]]
    
    # Calculate overlap
    overlap_genes <- intersect(region_signature_genes, cell_sig)
    overlap_count[i] <- length(overlap_genes)
    
    # Hypergeometric test parameters
    m <- length(cell_sig)  # Size of cell signature
    n <- universe_size - m  # Genes not in cell signature
    
    # Expected overlap by chance
    expected_overlap[i] <- (k * m) / universe_size
    
    # Determine direction and calculate p-value
    if (overlap_count[i] > expected_overlap[i]) {
      direction[i] <- "enriched"
      # Upper tail test for enrichment
      p_value[i] <- stats::phyper(overlap_count[i] - 1, m, n, k, lower.tail = FALSE)
      
    } else if (overlap_count[i] < expected_overlap[i]) {
      direction[i] <- "depleted"
      # Lower tail test for depletion
      p_value[i] <- stats::phyper(overlap_count[i], m, n, k, lower.tail = TRUE)
      
    } else {
      # Overlap equals expectation (rare case)
      direction[i] <- "neutral"
      p_value[i] <- 1.0
    }
    
    # Handle edge case where overlap is 0
    if (is.null(overlap_count[i]) || is.na(overlap_count[i])) {
      overlap_count[i] <- 0
    }
  }
  
  # Create results data frame
  results <- data.frame(
    cell.type = cell_types,
    overlap_count = overlap_count,
    expected_overlap = round(expected_overlap, 2),
    direction = direction,
    p.value = p_value,
    stringsAsFactors = FALSE
  )
  
  # Add -log10(p.value)
  results$`-log10(p.value)` <- -log10(results$p.value)
  
  # Add signed p-value (positive for enrichment, negative for depletion)
  signed_log_p <- results$`-log10(p.value)`
  signed_log_p[results$direction == "depleted"] <- -signed_log_p[results$direction == "depleted"]
  results$p.value.with.direction <- signed_log_p
  
  # Print summary
  cat("Overlap analysis complete!\n")
  cat("Enriched cell types (p < 0.05):", sum(results$direction == "enriched" & results$p.value < 0.05), "\n")
  cat("Depleted cell types (p < 0.05):", sum(results$direction == "depleted" & results$p.value < 0.05), "\n")
  cat("Total significant (p < 0.05):", sum(results$p.value < 0.05), "\n\n")
  
  # Export to Excel if requested
  if (!is.null(output_file)) {
    if (requireNamespace("rio", quietly = TRUE)) {
      rio::export(results, output_file)
      cat("Results exported to:", output_file, "\n")
    } else {
      warning("Package 'rio' not available. Results not exported. Install with: install.packages('rio')")
    }
  }
  
  return(results)
}


#' Calculate Multiple Region Signature Overlaps
#'
#' Convenience wrapper to calculate overlaps for multiple region signatures
#' against the same set of cell-type signatures.
#'
#' @param region_signature_list A named list where each element contains a character
#'   vector of gene symbols for a specific region (e.g., "fibroblast_foci", "alveoli").
#' @param cell_signature_list A list where each element contains a character vector
#'   of gene symbols for a specific cell type.
#' @param universe_size Integer specifying the total number of genes in the universe.
#'   Required parameter, no default value.
#' @param export_dir Character string specifying directory path for exporting results.
#'   If NULL (default), no files are exported.
#'
#' @return A named list of data frames, one for each region signature, containing
#'   overlap analysis results.
#'
#' @examples
#' \dontrun{
#' # Multiple region signatures
#' region_sigs <- list(
#'   fibroblast_foci = c("COL1A1", "COL1A2", "ACTA2"),
#'   alveoli = c("SFTPC", "SFTPA1", "LAMP3"),
#'   blood_vessel = c("VWF", "PECAM1", "CDH5")
#' )
#' 
#' # Cell signatures
#' cell_sigs <- list(
#'   Fibroblasts = c("COL1A1", "FN1", "DCN"),
#'   AT2_cells = c("SFTPC", "SFTPA1", "SFTPB"),
#'   Endothelial = c("VWF", "PECAM1", "CDH5")
#' )
#' 
#' # Calculate all overlaps
#' all_results <- calculate_multiple_overlaps(
#'   region_signature_list = region_sigs,
#'   cell_signature_list = cell_sigs,
#'   universe_size,      
#'   export_dir = "overlap_results"
#' )
#' }
#'
#' @export
calculate_multiple_overlaps <- function(region_signature_list,
                                       cell_signature_list,
                                       universe_size,
                                       export_dir = NULL) {
  
  if (!is.list(region_signature_list) || length(region_signature_list) == 0) {
    stop("region_signature_list must be a non-empty list")
  }
  
  region_names <- names(region_signature_list)
  if (is.null(region_names) || any(region_names == "")) {
    stop("region_signature_list elements must be named")
  }
  
  # Create export directory if needed
  if (!is.null(export_dir) && !dir.exists(export_dir)) {
    dir.create(export_dir, recursive = TRUE)
    cat("Created export directory:", export_dir, "\n")
  }
  
  # Process each region signature
  all_results <- list()
  
  for (region_name in region_names) {
    cat("\n========================================\n")
    cat("Processing region:", region_name, "\n")
    cat("========================================\n")
    
    # Set output file path if export directory is specified
    output_file <- if (!is.null(export_dir)) {
      file.path(export_dir, paste0(region_name, "_overlap_score.xlsx"))
    } else {
      NULL
    }
    
    # Calculate overlap
    results <- calculate_signature_overlap(
      region_signature_genes = region_signature_list[[region_name]],
      cell_signature_list = cell_signature_list,
      universe_size = universe_size,
      output_file = output_file
    )
    
    # Add region name to results
    results$region <- region_name
    
    # Store in list
    all_results[[region_name]] <- results
  }
  
  cat("\n========================================\n")
  cat("All overlap analyses complete!\n")
  cat("========================================\n")
  
  return(all_results)
}