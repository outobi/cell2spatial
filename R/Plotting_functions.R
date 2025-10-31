#' Create Signature Score Bubble Plot
#'
#' Creates a bubble plot to visualize gene signature scores across cell types and conditions.
#' This function converts wide format data to long format and creates publication-ready plots.
#'
#' @param signature_data A data.frame in wide format with cell types as rows and 
#'   conditions as columns. Values represent signature scores.
#' @param plot_type Character string specifying plot type. Options:
#'   "positive" (default): For positive signature scores (Figure 1D style)
#'   "negative": For negative signature scores (Figure S3A style)
#' @param size_limits Numeric vector of length 2 specifying size scale limits.
#'   Default: c(0, 160) for positive, c(-30, 0) for negative.
#' @param size_range Numeric vector of length 2 specifying bubble size range.
#'   Default: c(1, 10) for positive, c(10, 1) for negative.
#' @param size_breaks Numeric vector specifying legend breaks for size scale.
#' @param sample_column Character string specifying column name containing cell types.
#'   Default: "Sample" (assumes cell types are in a column named "Sample").
#' @param title Character string for plot title. Default: NULL (no title).
#' @param x_label Character string for x-axis label. Default: "Condition".
#' @param y_label Character string for y-axis label. Default: "Cell Type".
#'
#' @return A ggplot object that can be displayed or saved.
#'
#' @details
#' This function:
#' 1. Converts wide format data to long format using reshape2::melt()
#' 2. Creates a bubble plot where bubble size represents signature score magnitude
#' 3. Provides different scaling for positive vs negative signature scores
#' 4. Uses publication-ready formatting with bold text and clean theme
#'
#' @examples
#' # Create example signature data (wide format)
#' cell_types <- c("T Cells", "B Cells", "Macrophages", "Fibroblasts")
#' conditions <- c("Control alveoli", "IPF blood vessel", "IPF fibroblast foci")
#' 
#' # Positive signature scores
#' positive_data <- data.frame(
#'   Sample = cell_types,
#'   `Control alveoli` = c(50, 30, 80, 120),
#'   `IPF blood vessel` = c(60, 40, 90, 140),
#'   `IPF fibroblast foci` = c(70, 35, 95, 160),
#'   check.names = FALSE
#' )
#' 
#' # Create positive signature plot
#' p1 <- plot_signature_scores(positive_data, plot_type = "positive")
#' 
#' # Negative signature scores
#' negative_data <- data.frame(
#'   Sample = cell_types,
#'   `Control alveoli` = c(-10, -5, -20, -25),
#'   `IPF blood vessel` = c(-15, -8, -22, -28),
#'   `IPF fibroblast foci` = c(-12, -6, -18, -30),
#'   check.names = FALSE
#' )
#' 
#' # Create negative signature plot
#' p2 <- plot_signature_scores(negative_data, plot_type = "negative")
#'
#' @export
plot_signature_scores <- function(signature_data, 
                                 plot_type = "positive",
                                 size_limits = NULL,
                                 size_range = NULL, 
                                 size_breaks = NULL,
                                 sample_column = "Sample",
                                 title = NULL,
                                 x_label = "Condition",
                                 y_label = "Cell Type") {
  
  # Load required packages (add these to NAMESPACE imports)
  if (!requireNamespace("reshape2", quietly = TRUE)) {
    stop("Package 'reshape2' is required but not installed. Please install it with install.packages('reshape2')")
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but not installed. Please install it with install.packages('ggplot2')")
  }
  
  # Input validation
  if (!is.data.frame(signature_data)) {
    stop("signature_data must be a data.frame")
  }
  
  if (!sample_column %in% colnames(signature_data)) {
    stop(paste("Column", sample_column, "not found in signature_data"))
  }
  
  if (!plot_type %in% c("positive", "negative")) {
    stop("plot_type must be either 'positive' or 'negative'")
  }
  
  # Calculate data range from numeric columns (excluding sample column)
  numeric_columns <- sapply(signature_data, is.numeric)
  if (sample_column %in% names(signature_data)) {
    numeric_columns[sample_column] <- FALSE  # Exclude sample column from range calculation
  }
  
  if (sum(numeric_columns) == 0) {
    stop("No numeric columns found in signature_data")
  }
  
  numeric_data <- signature_data[, numeric_columns, drop = FALSE]
  data_range <- range(numeric_data, na.rm = TRUE)
  data_min <- data_range[1]
  data_max <- data_range[2]
  
  # Set default parameters based on plot type and actual data range
  if (is.null(size_limits)) {
    if (plot_type == "positive") {
      size_limits <- c(max(0, data_min), data_max)
    } else {
      size_limits <- c(data_min, min(0, data_max))
    }
  }
  
  if (is.null(size_range)) {
    size_range <- if (plot_type == "positive") c(1, 10) else c(10, 1)
  }
  
  if (is.null(size_breaks)) {
    # Create 4-5 breaks evenly distributed across the size_limits range
    break_min <- size_limits[1]
    break_max <- size_limits[2]
    
    if (plot_type == "positive") {
      # For positive: start from 0 (or data min if > 0) to data max
      n_breaks <- 5
      size_breaks <- seq(break_min, break_max, length.out = n_breaks)
      size_breaks <- round(size_breaks, 1)  # Round to 1 decimal place
    } else {
      # For negative: start from data min to 0 (or data max if < 0)
      n_breaks <- 5
      size_breaks <- seq(break_min, break_max, length.out = n_breaks)
      size_breaks <- round(size_breaks, 1)  # Round to 1 decimal place
    }
  }
  
  # Convert data frame from "wide" format to "long" format
  pcm <- reshape2::melt(signature_data, id = sample_column)
  pcm[[sample_column]] <- factor(pcm[[sample_column]], levels = unique(pcm[[sample_column]]))
  
  # Create the plot
  p <- ggplot2::ggplot(pcm, ggplot2::aes_string(y = sample_column, x = "variable")) + 
    ggplot2::geom_point(ggplot2::aes_string(size = "value", fill = "variable"), shape = 21) +
    ggplot2::scale_size_continuous(
      limits = size_limits, 
      range = size_range, 
      breaks = size_breaks,
      name = "Signature\nScore"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.key = ggplot2::element_blank(), 
      axis.text.x = ggplot2::element_text(
        colour = "black", 
        size = 12, 
        face = "bold", 
        angle = 45, 
        vjust = 1, 
        hjust = 1
      ), 
      axis.text.y = ggplot2::element_text(
        colour = "black", 
        face = "bold", 
        size = 11
      ),
      axis.title.x = ggplot2::element_text(
        colour = "black", 
        face = "bold", 
        size = 12
      ),
      axis.title.y = ggplot2::element_text(
        colour = "black", 
        face = "bold", 
        size = 12
      ),
      panel.grid.major = ggplot2::element_line(colour = "grey90"),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = title,
      x = x_label,
      y = y_label
    )
  
  # Add color scale to distinguish conditions
  n_conditions <- length(unique(pcm$variable))
  p <- p + ggplot2::scale_fill_manual(
    values = rainbow(n_conditions, alpha = 0.7),
    name = "Condition"
  )
  
  return(p)
}


#' Create Positive Signature Score Plot (Figure 1D Style)
#'
#' Convenience function to create a bubble plot for positive signature scores.
#' This replicates the style of Figure 1D from your publication.
#'
#' @param signature_data A data.frame in wide format with cell types and conditions.
#' @param title Character string for plot title. Default: "Positive Signature Scores".
#' @param custom_breaks Numeric vector for custom size scale breaks. 
#'   If NULL (default), breaks will be automatically calculated from data range.
#' @param max_value Numeric value for maximum size limit. 
#'   If NULL (default), will be automatically determined from data maximum.
#'
#' @return A ggplot object.
#'
#' @examples
#' # Create example data
#' data <- data.frame(
#'   Sample = c("T Cells", "B Cells", "Macrophages"),
#'   `Control alveoli` = c(50, 80, 120),
#'   `IPF fibroblast foci` = c(90, 110, 150),
#'   check.names = FALSE
#' )
#' 
#' p1 <- plot_positive_signature(data)
#'
#' @export
plot_positive_signature <- function(signature_data, 
                                   title = "Positive Signature Scores",
                                   custom_breaks = NULL,
                                   max_value = NULL) {
  
  # Calculate size_limits if max_value is provided
  size_limits <- if (!is.null(max_value)) c(0, max_value) else NULL
  
  return(plot_signature_scores(
    signature_data = signature_data,
    plot_type = "positive",
    size_limits = size_limits,
    size_breaks = custom_breaks,
    title = title
  ))
}


#' Create Negative Signature Score Plot (Figure S3A Style)
#'
#' Convenience function to create a bubble plot for negative signature scores.
#' This replicates the style of Figure S3A from your publication.
#'
#' @param signature_data A data.frame in wide format with cell types and conditions.
#' @param title Character string for plot title. Default: "Negative Signature Scores".
#' @param custom_breaks Numeric vector for custom size scale breaks.
#'   If NULL (default), breaks will be automatically calculated from data range.
#' @param min_value Numeric value for minimum size limit. 
#'   If NULL (default), will be automatically determined from data minimum.
#'
#' @return A ggplot object.
#'
#' @examples
#' # Create example data
#' data <- data.frame(
#'   Sample = c("T Cells", "B Cells", "Macrophages"),
#'   `Control alveoli` = c(-10, -15, -25),
#'   `IPF fibroblast foci` = c(-12, -18, -30),
#'   check.names = FALSE
#' )
#' 
#' p2 <- plot_negative_signature(data)
#'
#' @export
plot_negative_signature <- function(signature_data, 
                                   title = "Negative Signature Scores",
                                   custom_breaks = NULL,
                                   min_value = NULL) {
  
  # Calculate size_limits if min_value is provided
  size_limits <- if (!is.null(min_value)) c(min_value, 0) else NULL
  
  return(plot_signature_scores(
    signature_data = signature_data,
    plot_type = "negative", 
    size_limits = size_limits,
    size_breaks = custom_breaks,
    title = title
  ))
}