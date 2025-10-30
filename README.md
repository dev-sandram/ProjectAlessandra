## ProjectAlessandra

ProjectAlessandra is an R package designed to provide a unified, reproducible framework for **multi-omics data analysis** --- including bulk RNA-seq, clinical, ATAC-seq, variant, and single-cell integration datasets - and benchmarking. Each analytical workflow has been implemented in **two interchangeable paradigms**: - A **data.table** version for optimized performance on large datasets. - A **data.frame** version for compatibility and simplicity.

The key features are: - Dual implementation for clarity (data.frame) and performance (data.table) - Complete coverage of realistic multi-omics analysis steps, including bulk RNA-seq summarization and QC, clinical lab classification and temporal matching of vitals, ATAC-seq peak filtering and gene overlap detection, variant mapping to genes, multi-cohort integration and high-variance gene selection, single-cell cluster integration and visualization - Built-in benchmarking suite comparing computational efficiency between paradigms - Fully documented, reproducible functions for education and research

The package was developed as part of a complete data integration and exploration pipeline, focusing on **generalizable, reusable functions** with consistent interfaces and documentation.

## Installation

You can install the development version of ProjectAlessandra from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("dev-sandram/ProjectAlessandra")
```

Once installed:

``` r
# library(ProjectAlessandra)
```

## Example

This is a basic example which shows you how to solve a common problem: A minimal RNA-seq summarization workflow:

``` r
library(ProjectAlessandra)

# Summarize bulk RNA-seq counts
result <- bulk_counts_summary_dt(
  counts_path = "project_oct25/bulk_counts_long.csv",
  meta_path = "project_oct25/sample_metadata.csv"
)

#Outputs
filtered_data <- result$filtered_data
gene_mean_median <- result$gene_mean_median
gene_condition_means <- result$gene_condition_means

head(gene_mean_median)
```

# Documentation

The full workflow, including all 12 core tasks and the final single-cell integration step, is detailed in the vignette:

``` r
vignette("workflow_analysis", package = "ProjectAlessandra")
```

# Citation and Acknowledgements

Developed by Alessandra Monterosso, as part of the October 2025 multi-omics integration project --- a didactic and research exercise in modular R workflow design.

Built using: data.table, ggplot2, microbenchmark, usethis, devtools, roxygen2
