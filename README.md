# cell2spatial

[![R build status](https://github.com/outobi/cell2spatial/workflows/R-CMD-check/badge.svg)](https://github.com/your-username/cell2spatial/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/cell2spatial)](https://CRAN.R-project.org/package=cell2spatial)

## Overview

cell2spatial is an R package designed for comprehensive multiomics integrative analysis of spatial transcriptomics and single cell transcriptomics data to infer relative cell type enrichment in a histopathological region, with a focus on Idiopathic Pulmonary Fibrosis (IPF) research. It is also a powerful modality-agonist integration method that can combine LC/MS spatial proteomics with single cell transcriptomics. The package provides tools for region/cell type specific gene/protein extraction, query expression scoring analysis, overlap enrichment analysis and statistical comparisons. It can be smoothly expanded to other disease indications with spatial heterogeity. 

## Key Features

- **Spatial Region-Specific Differential Expression**: Compare gene/protein expression between different spatial regions (e.g., fibrotic foci vs. alveolar regions) using both parametric (t-test) and non-parametric (Wilcoxon) approaches. Calculates appropriate effect sizes (Cohen's d or rank-biserial correlation). It can be applied to both spatial transcriptomics like GeoMx or Visium (manual region combition) and spatial proteomics (LC/MS) data like laser capture microdissection.

- **Cell Type-Specific Differential Expression**: Identify genes that are differentially expressed in specific cell types compared to all other cell types within a diagnosis group. Uses Welch's t-test with Cohen's d effect sizes to quantify the magnitude of differences.

- **Gene Signature Scoring**: Query of region-specific differential gene expression in distinct cell types from scRNAseq to compute z-scores based on relative expression levels. Supports both single signature scoring and batch processing of multiple signatures.

- **Overlap Enrichment Analysis**: Determine whether two gene signatures have statistically significant overlap using hypergeometric tests. Useful for comparing region-specific signatures and cell type-specific signatures.

- **Flexible Multiple Testing Correction**: Choose between nominal p-values or automatic Benjamini-Hochberg FDR correction, giving you control over when and how to apply multiple testing correction in your analysis pipeline to extract region-specific and cell-type specific genes.

- **Data Normalization**: Built-in functions for z-score normalization and gene-wise scaling to prepare expression data for downstream analysis.

- **Publication-Ready Visualizations**: Generate customizable bubble plots for signature scores with automatic data-adaptive scaling.

## Installation

### Development Version

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install cell2spatial from GitHub
devtools::install_github("outobi/cell2spatial")
```

### Dependencies

The package requires R >= 4.0.0 and depends on several CRAN packages:

```r
install.packages(c("dplyr", "rstatix", "readxl", "rio", "ggplot2", 
                   "reshape2", "pheatmap", "RColorBrewer"))
```

## Core Functions

### Statistical Comparison Functions

#### `t_test_compare_cell()`
**Purpose**: Identify genes differentially expressed in a specific cell type compared to all other cell types in a diagnosis group.

**Use Case**: When you have single-cell and you want to find gene features that are specifically enriched in one cell type (e.g., identifying fibroblast-specific genes in IPF).

**Method**: Welch's t-test (does not assume equal variances) comparing target cell type vs. all others within a diagnosis group.

**Output**: Data frame with Cohen's d effect sizes, p-values (or FDR), mean expression values, and sample sizes. You can use the gene with max d effect sizes and p value/FDR < 0.05 to identify cell type specific genes. Try to avoid the case that one gene belongs to two or more cell types, which may comfound the following analysis.

**Key Parameters**:
- `normalized_count`: Expression matrix
- `meta_data`: Metadata with cell type and diagnosis columns
- `target_cell_type`: Cell type to analyze (e.g., "Fibroblasts")
- `diagnosis`: Condition to subset (e.g., "IPF")
- `correction_method`: "Nominal" or "BH" for FDR correction

---

#### `wilcox_compare_region()`
**Purpose**: Compare gene expression between a target spatial region and all other regions using non-parametric statistics.

**Use Case**: When comparing spatial regions with potentially non-normal distributions or when you want a robust, assumption-free test (e.g., comparing fibrotic foci vs. healthy alveolar regions).

**Method**: Wilcoxon rank-sum test (Mann-Whitney U) with rank-biserial correlation as effect size measure.

**Output**: Data frame with test statistics, p-values, effect sizes, and summary statistics for each gene. In the following analysis, try to modify cutoff power (effect sizes and p-values)to: 1 try to avoid the case that one gene belongs to two or more regions. 2 boost the number of region-specific genes in each region.

**When to Use**: 
- Small sample sizes per region
- Non-normal data distributions
- Presence of outliers


---

#### `t_test_compare_region()`
**Purpose**: Compare gene expression between spatial regions using parametric statistics.

**Use Case**: Similar to `wilcox_compare_region()` but assumes data follows normal distribution. Provides more statistical power when assumptions are met.

**Method**: Student's t-test (equal variance) or Welch's t-test (unequal variance) with Cohen's d effect size.

**Output**: Data frame with t-statistics, p-values, Cohen's d, means, and standard deviations. You can define the gene with max d effect sizes and p value/FDR < 0.05 as region specific genes. This can avoid the case that one gene belongs to two or more regions.

**When to Use**:
- Larger sample sizes (e.g. n ≥ 30 per group)
- Approximately normal distributions


---

### Overlap Enrichment Functions

#### `calculate_signature_overlap()`
**Purpose**: Test whether two gene signature sets have statistically significant overlap.

**Use Case**: 
- Comparing a region-specific gene signature set and multiple cell-type specific gene signatures (e.g., do fibroblast cells, immule cells, and epithelial cells enriched in foci areas?)

**Method**: Hypergeometric test (Fisher's exact test for enrichment) with fold enrichment calculation.

**Output**: Single row data frame with:
- `overlap`: Number of genes in common
- `p.value`: Hypergeometric p-value
- `enrichment_type`: "enriched" or "depleted"
- Contingency table values

**Key Parameters**:
- `region_signature_genes`: First gene list (character vector)
- `query_signature_genes`: Second gene list (character vector)
- `universe_size`: Total number of genes in background (integer)

**Statistical Note**: Uses `phyper()` to calculate the probability of observing the overlap by chance given the universe size.

---

#### `calculate_multiple_overlaps()`
**Purpose**: Batch process overlap enrichment for multiple cell type specific gene signatures against multiple region specific gene signature.

**Use Case**: Compare multiple region-specific signatures against multiple cell-type-specific signatures. This is useful when you want to compare multiple regions.

**Output**: Data frame with one row per comparison, including all statistics from `calculate_signature_overlap()`.

---

### Gene Signature Scoring Functions

#### `calculate_signature_score()`
**Purpose**: Query of region-specific differential gene expression in distinct cell types from scRNAseq to compute z-scores based on relative expression levels. 
**Use Case**: 
- Scoring region-derived signature summed expression in all cell type 
- Quantifying cell type relative enrichment in spatial regions

**Method**: 
1. Subset expression matrix to region signature genes
2. Calculate z-scores if not already normalized
3. Compute mean z-score across all cells in the cell type for the given region
4. Return enrichment score per cell type

**Output**: Named numeric vector of enrichment scores (one per cell type for the given region).

**Interpretation**: 
- Positive scores: Signature is up-regulated
- Negative scores: Signature is down-regulated
- Magnitude: Strength of enrichment

---

#### `calculate_multiple_signature_scores()`
**Purpose**: Score multiple region gene signatures simultaneously across all cell types.

**Use Case**: Efficiently calculate scores for many cell type signatures at once (e.g., scoring all cell types from single-cell atlas).

**Method**: Vectorized version of `calculate_signature_score()` applied to a list of signatures.

**Output**: Matrix with cell type as rows and region as columns, filled with enrichment scores.

**Advantages**: 
- Much faster than looping
- Returns organized matrix for downstream analysis
- Easy to visualize with heatmaps

---

### Data Normalization Functions

#### `normalization_by_gene()`
**Purpose**: Perform gene-wise z-score normalization across samples.

**Method**: For each gene, subtract mean and divide by standard deviation across all samples: `z = (x - mean) / sd`

**Use Case**: Prepare expression data for signature scoring or when genes have different scales.

---

#### `normalize_expression()`
**Purpose**: General normalization wrapper supporting multiple methods.

**Methods Available**: Z-score (default), quantile normalization, log transformation.

---

### Visualization Functions

#### `plot_signature_scores()`
**Purpose**: Create publication-ready bubble plots showing signature enrichment across samples/regions.

**Features**:
- Automatic data-adaptive scaling
- Size indicates magnitude
- Customizable themes and labels

**Use Case**: Visualize which cell type signatures are enriched in which spatial regions.

## Data Structure

The package is designed to work with spatial transcriptomics/proteomics data (e.g., GeoMx DSP, Visium, CosMx) and single-cell RNA-seq data. It expects:

### Expression Data Format
- **Matrix/Data Frame**: Genes as rows, samples as columns
- **Row Names**: Gene symbols (e.g., "COL1A1", "ACTA2")
- **Column Names**: Sample identifiers (e.g., "S1", "S2", "Sample_001")
- **Values**: Normalized expression counts (not raw counts)
- **Example**:
```r
           S1    S2    S3    S4
COL1A1   12.5  15.2  10.8  14.3
ACTA2     8.3   9.1   7.5   8.8
MMP7      5.2   6.1   4.8   5.9
```

### Metadata Format
Required columns depend on the analysis:

**For Cell Type Analysis** (`t_test_compare_cell`):
- Sample identifiers matching expression data columns
- `celltype` column: Cell type annotations (e.g., "Fibroblasts", "AT1", "Macrophages")
- `Diagnosis` column: Disease status (e.g., "IPF", "Control", "COPD")

**For Region Comparison** (`wilcox_compare_region`, `t_test_compare_region`):
- Column 1: Sample identifiers (e.g., "S1", "S2")
- Column 2: Region types (e.g., "foci", "alveoli", "fibrosis")

**Example Metadata**:
```r
# For cell type analysis
           celltype  Diagnosis
S1      Fibroblasts       IPF
S2      Fibroblasts   Control
S3             AT1        IPF
S4       Macrophages      IPF

# For region comparison
  sample_id region_type
1        S1        foci
2        S2     alveoli
3        S3    fibrosis
```

### Gene Signature Lists
- **Format**: Character vectors of gene symbols
- **Source**: Can be derived from single-cell RNA-seq marker genes, literature, or pathway databases
- **Example**:
```r
fibroblast_signature <- c("COL1A1", "COL1A2", "ACTA2", "FN1", "VIM")
epithelial_signature <- c("EPCAM", "CDH1", "KRT8", "KRT18")
```

## Complete Analysis Workflow

This section demonstrates a typical cell2spatial analysis pipeline from data loading to result visualization.

### Workflow Overview
1. **Load and prepare data** - Import expression matrices and metadata
2. **Identify region-specific genes** - Compare spatial regions to find differentially expressed genes
3. **Score cell type signatures** - Project single-cell signatures onto spatial data
4. **Perform overlap enrichment** - Test whether spatial signatures overlap with cell type markers
5. **Visualize results** - Create publication-ready figures

### Step-by-Step Example

```r
library(cell2spatial)
library(rio)  # For Excel export

# Load your data
load("your_spatial_data.RData")
load("your_scrna_data.RData")
# Assumes you have:
# Spatial data:
# - expression_data: genes × samples matrix (spatial proteomics/transcriptomics)
# - sample_metadata: sample annotations with region information
# Single-cell data:
# - scRNA_expression_data: genes × cells matrix
# - scRNA_metadata: cell annotations with celltype and diagnosis columns

# ============================================
# 1. Identify region-specific genes from spatial data
# ============================================
region_signatures <- list()
for (region in unique(sample_metadata[, 2])) {
  # Option 1: Use Wilcoxon test with nominal p-values (default)
  results_wilcox <- wilcox_compare_region(
    expression_data = expression_data,
    sample_metadata = sample_metadata,
    target_region = region,
    correction_method = "Nominal"
  )
  
  # Option 2: Use t-test with Benjamini-Hochberg FDR correction
  results_ttest <- t_test_compare_region(
    expression_data = expression_data,
    sample_metadata = sample_metadata,
    target_region = region,
    correction_method = "BH"
  )
  
  # Export results
  rio::export(results_wilcox, paste0("Region_", region, "_wilcox_nominal.xlsx"))
  rio::export(results_ttest, paste0("Region_", region, "_ttest_FDR.xlsx"))
  
  # Extract region-specific genes (FDR < 0.05, |d| > 0.5)
  region_signatures[[region]] <- results_ttest[results_ttest$p.value < 0.05 & 
                                                abs(results_ttest$d) > 0.5, "gene_id"]
}

# Store individual region signatures for convenience
sig_genes_foci <- region_signatures$foci
sig_genes_alveoli <- region_signatures$alveoli
sig_genes_fibrosis <- region_signatures$fibrosis

# ============================================
# 2. Query method - Score region signatures in single-cell data
# ============================================
# Calculate enrichment scores for region-specific genes across cell types
signature_scores <- calculate_multiple_signature_scores(
  expression_matrix = scRNA_expression_data,  # Single-cell expression matrix
  signature_list = region_signatures,  # List of region-specific gene signatures
  normalize = TRUE
)

# Export signature scores
rio::export(signature_scores, "region_signature_scores_in_cell_types.csv")

# ============================================
# 3. Identify cell type-specific genes from single-cell data
# ============================================
cell_type_signatures <- list()
for (cell_type in unique(scRNA_metadata$celltype)) {
  # Use t-test with FDR correction to identify cell type markers
  cell_results <- t_test_compare_cell(
    normalized_count = scRNA_expression_data,  # Single-cell expression data
    meta_data = scRNA_metadata,  # Single-cell metadata
    target_cell_type = cell_type,
    diagnosis = "IPF",  # or "All" if you want pan-disease markers
    correction_method = "BH"
  )
  
  # Extract cell type-specific genes (FDR < 0.05, 0.5 < |d| < 10)
  cell_type_signatures[[cell_type]] <- cell_results[cell_results$stat.sig < 0.05 & 
                                                      abs(cell_results$d) > 0.5 & 
                                                      abs(cell_results$d) < 10, "gene_id"]
  
  # Export cell type results
  rio::export(cell_results, paste0("Cell_type_", cell_type, "_FDR.xlsx"))
}

# ============================================
# 4. Overlap enrichment analysis - Compare regions with cell types
# ============================================
# Test multiple regions against multiple cell types
all_overlaps <- list()
for (region_name in names(region_signatures)) {
  region_genes <- region_signatures[[region_name]]
  
  # Skip if no genes found for this region
  if (length(region_genes) == 0) next
  
  # Test overlap with all cell type signatures
  overlap_results <- calculate_multiple_overlaps(
    region_signature_list = list(region = region_genes),
    query_signature_genes = cell_type_signatures,
    universe_size = nrow(expression_data)
  )
  
  overlap_results$region <- region_name
  all_overlaps[[region_name]] <- overlap_results
}

# Combine all overlap results
all_overlaps_df <- do.call(rbind, all_overlaps)

# Filter for significant enrichments
sig_overlaps <- all_overlaps_df[all_overlaps_df$p.value < 0.05 & 
                                 all_overlaps_df$enrichment_type == "enriched", ]

# Display and export results
print(sig_overlaps)
rio::export(sig_overlaps, "significant_region_celltype_overlaps.xlsx")
rio::export(all_overlaps_df, "all_region_celltype_overlaps.xlsx")

# ============================================
# 5. Visualization - Plot signature scores and overlap results
# ============================================

# Plot 1: Query method results - Region signature scores across cell types to dplot positive enrichment cell type by default
plot_signature_scores(
  signature_data = signature_scores,
  title = "Region-Specific Gene Signatures Enrichment Across Cell Types"
)
ggsave("region_signature_scores_heatmap.pdf", width = 10, height = 8)

# Plot 2: Overlap enrichment results - Using plot_signature_scores
# Prepare overlap matrix for plotting
overlap_matrix <- reshape2::dcast(sig_overlaps, region ~ query_signature, value.var = "overlap")
rownames(overlap_matrix) <- overlap_matrix$region
overlap_matrix$region <- NULL

# Plot overlap counts using plot_signature_scores function to plot positive enrichment cell type by default
plot_signature_scores(
  signature_data = overlap_matrix,
  title = "Gene Overlap Counts: Regions vs Cell Types"
)
ggsave("region_celltype_overlap_bubble.pdf", width = 12, height = 8)




---
### Understanding Effect Sizes

#### Cohen's d (parametric tests)
- **Small effect**: d ≈ 0.2
- **Medium effect**: d ≈ 0.5
- **Large effect**: d ≈ 0.8 or higher
- **Interpretation**: Number of standard deviations between group means
- **Sign**: Positive = target group higher, Negative = target group lower

#### Rank-Biserial Correlation (non-parametric tests)
- **Range**: -1 to +1
- **Interpretation**: Similar to Cohen's d but uses ranks instead of means
- **Thresholds**: |r| > 0.3 (small), > 0.5 (medium), > 0.7 (large)
- **Advantage**: Robust to outliers and non-normal distributions


### Reading Overlap Enrichment Results

From `calculate_signature_overlap()`:

- **overlap**: Number of genes shared between two signatures
- **p.value**: Probability of observing this overlap by chance
  - p < 0.05: Significant enrichment or depletion
- **enrichment_type**: 
  - "enriched": More overlap than expected by chance
  - "depleted": Less overlap than expected (rare)
- **fold_enrichment**: Observed overlap / expected overlap
  - > 1: Enriched
  - < 1: Depleted
  - = 1: No enrichment

**Example Interpretation**:
```
overlap = 45, p.value = 1.2e-10, enrichment_type = "enriched", fold_enrichment = 3.2
```
"45 genes are shared between signatures, which is 3.2 times more than expected by chance (p = 1.2e-10). This indicates strong biological overlap."

---

## Analysis Best Practices

### Data Preparation

1. **Normalization**: Always use normalized data (not raw counts)
   - For spatial data: Use built-in platform normalization (GeoMx Q3, Visium SCTransform)
   - For signature scoring: Use z-score normalization (`normalization_by_gene()`)

2. **Quality Control**: Remove low-quality samples before analysis
   - Check for technical artifacts
   - Remove samples with very low gene detection
   - Ensure adequate biological replicates (≥3 per group recommended)

3. **Gene Filtering**: Consider pre-filtering genes
   - Remove genes with zero variance
   - Filter for expressed genes (present in ≥50% of samples)
   - For overlap tests, ensure universe_size reflects filtered gene count

### Statistical Test Selection

**Use Wilcoxon (`wilcox_compare_region`) when:**
- Sample size is small (n < 30 per group)
- Data distribution is skewed
- Outliers are present
- You want robust, assumption-free results

**Use t-test (`t_test_compare_region`) when:**
- Sample size is adequate (n ≥ 30 per group)
- Data is approximately normally distributed
- You need Cohen's d for meta-analysis
- You want maximum statistical power

**Use Welch's t-test (default in `t_test_compare_cell`) when:**
- Group variances are unequal
- This is the safer parametric choice

### Multiple Correction Options

All statistical comparison functions support flexible multiple testing correction via the `correction_method` parameter:

#### **`correction_method = "Nominal"`** (default)
- Returns unadjusted p-values
- **When to Use**: 
  - Exploratory analysis
  - When you want to apply custom correction later
  - When combining results across multiple analyses
  - For effect size-focused analyses (p-values are secondary)
- **Advantage**: Maximum flexibility, no information loss

#### **`correction_method = "BH"`**
- Applies Benjamini-Hochberg FDR correction automatically
- Controls false discovery rate at specified level (typically 0.05)
- **When to Use**:
  - Final analysis for publication
  - When you need immediate corrected values
  - Standard differential expression analysis
- **Interpretation**: FDR < 0.05 means 5% of significant genes are expected to be false positives

### Functions Supporting Correction
- `t_test_compare_cell()` - Cell type-specific differential expression
- `wilcox_compare_region()` - Wilcoxon tests between regions
- `t_test_compare_region()` - t-tests between regions

### Best Practices
1. **Start with Nominal**: Begin exploratory analysis with nominal p-values
2. **Apply FDR for Final Results**: Use BH correction for results you'll publish
3. **Report Effect Sizes**: Always report effect sizes (Cohen's d, rank-biserial) alongside p-values
4. **Set Thresholds**: Consider both statistical significance (p < 0.05) and biological significance (effect size > 0.5) to identify region- and cell type-specific genes.




### Signature Scoring Recommendations

1. **Signature Size**: Optimal 10-100 genes
   - Too small (<10): Unstable scores, high variance
   - Too large (>200): Diluted signal, includes noise

2. **Gene Quality**: Use high-confidence markers
   - Cell type-specific genes (not broadly expressed)
   - Genes with strong effect sizes in source data
   - Avoid housekeeping genes

3. **Interpretation**: Compare scores across samples
   - Scores are relative (z-scores)
   - Compare within same experiment/batch
   - Use same normalization method throughout

### Troubleshooting Common Issues

#### Low Statistical Power
- **Symptom**: Many genes show large effect sizes but non-significant p-values
- **Solutions**: 
  - Increase sample size
  - Use less stringent FDR threshold (0.10)
  - Focus on effect size rankings
  - Validate top candidates with independent data

#### Too Many Significant Genes
- **Symptom**: Thousands of genes pass FDR < 0.05
- **Solutions**:
  - Check data normalization (may have batch effects)
  - Apply stricter effect size thresholds (|d| > 0.8)
  - Use more stringent FDR (< 0.01)
  - Verify groups are biologically comparable

#### No Overlap Enrichment Detected
- **Symptom**: Overlap p-values are all > 0.05
- **Solutions**:
  - Check universe size (should match filtered gene count)
  - Verify gene identifier consistency (symbols vs. Ensembl IDs)
  - Use larger gene lists (top 200-500 genes)
  - Consider biological differences between platforms
  - some region just do not have representative up-regulated genes, we do not recommend using this workflow in this region then.

#### Inconsistent Results Between Tests
- **Symptom**: Wilcoxon and t-test give different results
- **Explanations**:
  - Different assumptions: Wilcoxon tests ranks, t-test tests means
  - Outliers: Affect t-tests more than Wilcoxon
  - Non-normal data: Violates t-test assumptions
- **Action**: Examine data distributions, report both results, prefer Wilcoxon for robustness

---

## Documentation

Comprehensive documentation is available for all functions:

```r
# View function documentation
?t_test_compare_cell
?wilcox_compare_region
?t_test_compare_region
?calculate_signature_overlap
?calculate_signature_score

# View package vignettes
browseVignettes("cell2spatial")
```

## Citation

If you use cell2spatial in your research, please cite:

```
Wang, F., Jin, L., Wang, X., Cui, B., Yang, Y., Duggan, L., Schwartz Sterman, A., Lloyd, S. M., Hazelwood, L. A., Chaudhary, N., Bawa, B., Phillips, L. A., He, Y., & Tian, Y. (2025). Novel Integration of Spatial and Single-Cell Omics Data Sets Enables Deeper Insights into IPF Pathogenesis. Proteomes, 13(1), 3. https://doi.org/10.3390/proteomes13010003
```

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Feel free to fork the repository and submit pull requests.

## Issues

Please report bugs and feature requests at: https://github.com/outobi/cell2spatial/issues

## Development

This package is actively developed for IPF research applications. Another implement in Rheumatoid Athritis is also available in Proteomes 2025, 13(2), 17; https://doi.org/10.3390/proteomes13020017. For questions or collaboration opportunities, please contact Fei Wang at wf199209@gmail.com.