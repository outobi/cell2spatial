# cell2spatial

[![R build status](https://github.com/ut/cell2spatial/workflows/R-CMD-check/badge.svg)](https://github.com/your-username/cell2spatial/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/cell2spatial)](https://CRAN.R-project.org/package=cell2spatial)

## Overview

cell2spatial is an R package designed for comprehensive multiomics integrative analysis of spatial transcriptomics/proteomics and single cell transcriptomics data to infer relative cell type enrichment in a histopathological region, with a focus on Idiopathic Pulmonary Fibrosis (IPF) research. It is a powerful modality-agonist integration method that can combine LC/MS spatial proteomics with single cell transcriptomics. The package provides tools for region/cell type specific gene extraction, query expression scoring analysis, overlap enrichment analysis and statistical comparisons. It can be smoothly expanded to other disease indications with spatial heterogeity. 

## Key Features

- **Cell Type-Specific Differential Expression**: Identify genes that are differentially expressed in specific cell types compared to all other cell types within a diagnosis group. Uses Welch's t-test with Cohen's d effect sizes to quantify the magnitude of differences.

- **Spatial Region Comparison**: Compare gene expression between different spatial regions (e.g., fibrotic foci vs. alveolar regions) using both parametric (t-test) and non-parametric (Wilcoxon) approaches. Calculates appropriate effect sizes (Cohen's d or rank-biserial correlation).

- **Gene Signature Scoring**: Calculate enrichment scores for gene signatures (e.g., from single-cell RNA-seq) across spatial regions. Supports both single signature scoring and batch processing of multiple signatures.

- **Overlap Enrichment Analysis**: Determine whether two gene signatures have statistically significant overlap using hypergeometric tests. Useful for comparing region-specific signatures or validating findings across different spatial modalities.

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

## Quick Start

```r
library(cell2spatial)

# Example: Cell type-specific analysis with nominal p-values
cell_analysis <- t_test_compare_cell(normalized_count = count_data,
                                     meta_data = metadata,
                                     target_cell_type = "AT1",
                                     diagnosis = "IPF",
                                     correction_method = "Nominal")

# Example: Cell type analysis with FDR correction
cell_analysis_fdr <- t_test_compare_cell(normalized_count = count_data,
                                         meta_data = metadata,
                                         target_cell_type = "AT1",
                                         diagnosis = "IPF",
                                         correction_method = "BH")

# Example: Region comparison with Wilcoxon test
region_wilcox <- wilcox_compare_region(expression_data = expr_data,
                                       sample_metadata = sample_meta,
                                       target_region = "foci",
                                       correction_method = "BH")

# Example: Region comparison with t-test
region_ttest <- t_test_compare_region(expression_data = expr_data,
                                      sample_metadata = sample_meta,
                                      target_region = "foci",
                                      correction_method = "Nominal")
```

## Core Functions

### Statistical Comparison Functions

#### `t_test_compare_cell()`
**Purpose**: Identify genes differentially expressed in a specific cell type compared to all other cell types.

**Use Case**: When you have single-cell or deconvolved spatial data annotated with cell types, and you want to find genes that are specifically enriched or depleted in one cell type (e.g., identifying fibroblast-specific genes in IPF).

**Method**: Welch's t-test (does not assume equal variances) comparing target cell type vs. all others within a diagnosis group.

**Output**: Data frame with Cohen's d effect sizes, p-values (or FDR), mean expression values, and sample sizes. You can use the gene with max d effect sizes and p value/FDR < 0.05 to identify cell type specific genes. Avoid the case that one gene belongs to two or more cell types.

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

**Output**: Data frame with test statistics, p-values, effect sizes, and summary statistics for each gene. Modify cutoff power to 1 avoid the case that one gene belongs to two or more regions. 2 boost the number of region-specific genes.

**When to Use**: 
- Small sample sizes per region
- Non-normal data distributions
- Presence of outliers


---

#### `t_test_compare_region()`
**Purpose**: Compare gene expression between spatial regions using parametric statistics.

**Use Case**: Similar to `wilcox_compare_region()` but assumes data follows normal distribution. Provides more statistical power when assumptions are met.

**Method**: Student's t-test (equal variance) or Welch's t-test (unequal variance) with Cohen's d effect size.

**Output**: Data frame with t-statistics, p-values, Cohen's d, means, and standard deviations. You can use the gene with max d effect sizes and p value/FDR < 0.05 to identify region specific genes. Avoid the case that one gene belongs to two or more regions.

**When to Use**:
- Larger sample sizes (n ≥ 30 per group)
- Approximately normal distributions


---

### Overlap Enrichment Functions

#### `calculate_signature_overlap()`
**Purpose**: Test whether two gene signatures have statistically significant overlap.

**Use Case**: 
- Comparing region-specific gene signatures and cell-type specific gene signatures (e.g., do fibroblast cells enriched in foci areas)

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
**Purpose**: Batch process overlap enrichment for multiple cell type specific gene signatures against a region specific gene signature.

**Use Case**: Compare multiple region-specific signatures against a single reference (e.g., testing multiple cell markers against a region gene set).

**Output**: Data frame with one row per comparison, including all statistics from `calculate_signature_overlap()`.

---

### Gene Signature Scoring Functions

#### `calculate_signature_score()`
**Purpose**: Calculate enrichment scores for a gene signature across samples based on expression levels.

**Use Case**: 
- Scoring single-cell-derived signatures in region specific gene set
- Quantifying cell type abundance in spatial regions

**Method**: 
1. Subset expression matrix to region signature genes
2. Calculate z-scores if not already normalized
3. Compute mean z-score across all cells in the cell type for each region
4. Return enrichment score per cell type

**Output**: Named numeric vector of enrichment scores (one per cell type for each region).

**Interpretation**: 
- Positive scores: Signature is up-regulated
- Negative scores: Signature is down-regulated
- Magnitude: Strength of enrichment

---

#### `calculate_multiple_signature_scores()`
**Purpose**: Score multiple region gene signatures simultaneously across all cell types.

**Use Case**: Efficiently calculate scores for many cell type signatures at once (e.g., scoring all 31 cell types from single-cell atlas).

**Method**: Vectorized version of `calculate_signature_score()` applied to a list of signatures.

**Output**: Matrix with samples as rows and signatures as columns, filled with enrichment scores.

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
- Color indicates enrichment direction (red = up, blue = down)
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
load("your_geomx_data.RData")
# Assumes you have:
# - expression_data: genes × samples matrix
# - metadata: sample annotations
# - scRNA_signatures: list of cell type marker genes

# ============================================
# 1. Identify region-specific genes
# ============================================
for (cell_type in unique(metadata$celltype)) {
  # Option 1: Use nominal p-values (default)
  results <- t_test_compare_cell(
    normalized_count = expression_data,
    meta_data = metadata,
    target_cell_type = cell_type,
    diagnosis = "IPF",
    correction_method = "Nominal"
  )
  
  # Option 2: Apply Benjamini-Hochberg FDR correction
  results_fdr <- t_test_compare_cell(
    normalized_count = expression_data,
    meta_data = metadata,
    target_cell_type = cell_type,
    diagnosis = "IPF",
    correction_method = "BH"
  )
  
  # Export results
  rio::export(results, paste0("Cell_type_", cell_type, "_nominal.xlsx"))
  rio::export(results_fdr, paste0("Cell_type_", cell_type, "_FDR.xlsx"))
}

# 2. Perform overlap enrichment between regions
# Extract significant genes (p < 0.05 with effect size threshold)
sig_genes_foci <- results_fdr[results_fdr$stat.sig < 0.05 & abs(results_fdr$d) > 0.5, "gene_id"]

# Compare with another region's signature
overlap_result <- calculate_signature_overlap(
  region_signature_genes = sig_genes_foci,
  query_signature_genes = other_region_genes,
  universe_size = nrow(expression_data)
)

# 3. Score cell type signatures in spatial data
# Calculate enrichment for all cell type signatures
signature_scores <- calculate_multiple_signature_scores(
  expression_matrix = expression_data,
  signature_list = scRNA_signatures,  # List of gene vectors
  normalize = TRUE
)

# 4. Test overlap between spatial and cell type signatures
# Which cell types are enriched in fibrotic foci?
overlap_results <- calculate_multiple_overlaps(
  region_signature_list = list(
    foci = sig_genes_foci,
    alveoli = sig_genes_alveoli,
    fibrosis = sig_genes_fibrosis
  ),
  query_signature_genes = scRNA_signatures$Fibroblasts,
  universe_size = nrow(expression_data)
)

# Filter for significant overlaps
sig_overlaps <- overlap_results[overlap_results$p.value < 0.05 & 
                                overlap_results$enrichment_type == "enriched", ]
print(sig_overlaps)

# 5. Visualize signature scores
plot_signature_scores(
  signature_data = signature_scores,
  title = "Cell Type Enrichment Across Spatial Regions"
)

# Export results
rio::export(sig_overlaps, "significant_overlaps.xlsx")
rio::export(signature_scores, "signature_scores.csv")
```

### Common Analysis Scenarios

#### Scenario 1: Single Region vs. All Others
```r
# Compare fibrotic foci to all other regions
foci_genes <- wilcox_compare_region(
  expression_data = expr_data,
  sample_metadata = sample_meta,
  target_region = "foci",
  correction_method = "BH"
)

# Extract significant genes with strong effects
sig_foci <- foci_genes[foci_genes$p.value < 0.05 & 
                       abs(foci_genes$effect_size) > 0.5, ]
```

#### Scenario 2: Cell Type Enrichment in Disease
```r
# Find fibroblast-specific genes in IPF
ipf_fibroblasts <- t_test_compare_cell(
  normalized_count = expr_data,
  meta_data = metadata,
  target_cell_type = "Fibroblasts",
  diagnosis = "IPF",
  correction_method = "BH"
)

# Compare to control fibroblasts
ctrl_fibroblasts <- t_test_compare_cell(
  normalized_count = expr_data,
  meta_data = metadata,
  target_cell_type = "Fibroblasts",
  diagnosis = "Control",
  correction_method = "BH"
)

# Find IPF-specific changes
ipf_specific <- ipf_fibroblasts[!(ipf_fibroblasts$gene_id %in% 
                                   ctrl_fibroblasts$gene_id[ctrl_fibroblasts$stat.sig < 0.05]), ]
```

#### Scenario 3: Validate Spatial Findings with Single-Cell
```r
# Get region-specific genes from spatial data
spatial_foci_genes <- wilcox_compare_region(...)$gene_id[1:200]

# Get cell type markers from single-cell
sc_fibroblast_markers <- scRNA_signatures$Fibroblasts

# Test overlap
validation <- calculate_signature_overlap(
  region_signature_genes = spatial_foci_genes,
  query_signature_genes = sc_fibroblast_markers,
  universe_size = nrow(expression_data)
)

if(validation$p.value < 0.05 & validation$enrichment_type == "enriched") {
  cat("Spatial foci genes are enriched for fibroblast markers!\n")
  cat("Overlap:", validation$overlap, "genes\n")
}
```

## Multiple Testing Correction

### Why This Matters
When testing thousands of genes simultaneously, some will appear significant by chance (false positives). Multiple testing correction adjusts p-values to control the false discovery rate.

### Correction Options

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
4. **Set Thresholds**: Consider both statistical significance (p < 0.05) and biological significance (effect size > 0.5)

---

## Interpreting Results

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

### Statistical Significance Guidelines

1. **P-value Thresholds**:
   - p < 0.05: Statistically significant
   - p < 0.01: Highly significant
   - p < 0.001: Very highly significant

2. **FDR Thresholds** (when using `correction_method = "BH"`):
   - FDR < 0.05: 5% expected false positives (standard)
   - FDR < 0.10: 10% expected false positives (more lenient for discovery)
   - FDR < 0.01: 1% expected false positives (very stringent)

3. **Combined Criteria** (recommended):
   - **High confidence**: p < 0.05 AND |effect size| > 0.5
   - **Strong candidates**: p < 0.01 AND |effect size| > 0.8
   - **Discovery mode**: p < 0.10 AND |effect size| > 0.3

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

### Multiple Testing Strategy

1. **Exploratory Phase**: Use `correction_method = "Nominal"`
   - Examine effect size distributions
   - Identify patterns and trends
   - Generate hypotheses

2. **Validation Phase**: Use `correction_method = "BH"`
   - Apply to final gene lists for publication
   - Report both nominal and FDR-corrected values
   - Consider pre-filtering to reduce multiple testing burden

3. **Reporting**: Always report:
   - Effect sizes with confidence intervals
   - Both nominal and corrected p-values
   - Number of tests performed
   - Filtering criteria applied

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

This package is actively developed for IPF research applications. For questions or collaboration opportunities, please contact Fei Wang at wf199209@gmail.com.