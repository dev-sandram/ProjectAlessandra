# 🧬 ProjectAlessandra

**ProjectAlessandra** is an R package designed to build, benchmark, and demonstrate modular data analysis workflows for **multi-omics datasets** — including RNA-seq, ATAC-seq, genomic variants, and clinical metadata — using both base **data.frame** and high-performance **data.table** paradigms.

The package was created as part of a didactic and research exercise (October 2025) to:

- Develop **modular R scripts** for realistic omics workflows  
- Handle **heterogeneous biomedical datasets** reproducibly  
- Compare **data.frame vs data.table** performance in analytical tasks  
- Provide **step-by-step examples** for teaching, benchmarking, and exploration  

---

## 🚀 Features

Each function in the package corresponds to a realistic step of an omics data analysis pipeline, covering:

| Category | Example Functions | Description |
|-----------|------------------|--------------|
| **RNA-seq QC and summarization** | `bulk_counts_summary_dt()`, `bulk_counts_qc_dt()` | Summarize and transform raw RNA-seq counts |
| **Metadata integration** | `annotate_counts_dt()`, `combine_cohorts_dt()` | Merge counts and metadata for cohort-level analyses |
| **Benchmarking and performance** | `subset_counts_dt()` | Compare indexed vs non-indexed data.table operations |
| **Clinical data processing** | `classify_labs_dt()`, `match_vitals_dt()` | Integrate and correlate lab and vital sign data |
| **ATAC-seq and variants** | `top_peaks_dt()`, `atac_to_gene_dt()`, `variants_to_genes_dt()` | Map peaks and variants to genes efficiently |
| **Exploratory filtering** | `gene_stats_filter_dt()` | Identify upregulated or high-variance genes |
| **Single-cell integration** | `final_revision_dt()` | Quantify cell type distribution across clusters/tissues |

Each task is implemented in two versions:
- a **data.frame** version for clarity and learning  
- a **data.table** version for scalability and benchmarking  

---

## 🧩 Installation

You can install the package directly from GitHub:

```r
# install.packages("devtools")
devtools::install_github("dev-sandram/ProjectAlessandra")
Then load it as usual:

r
library(ProjectAlessandra)
📘 Usage Overview
A typical workflow involves:

r
# Summarize RNA-seq counts
result <- bulk_counts_summary_dt(
  "project_oct25/data/bulk_counts_long.csv",
  "project_oct25/data/sample_metadata.csv"
)

# Quality control transformation
qc_data <- bulk_counts_qc_dt("project_oct25/data/bulk_counts_long.csv")

# Match clinical labs and vitals
matched <- match_vitals_dt(
  "project_oct25/data/clinical_labs.csv",
  "project_oct25/data/vitals_time_series.csv"
)

# Map ATAC-seq peaks to genes
peaks_to_genes <- atac_to_gene_dt(
  "project_oct25/data/atac_peaks.bed.csv",
  "project_oct25/data/gene_annotation.bed.csv"
)
Detailed examples and end-to-end workflows are available in the vignette:
vignettes/workflow_analysis.Rmd

Benchmarking
Each function includes two parallel implementations (*_df() and *_dt()) allowing reproducible benchmarking via microbenchmark:

r
library(microbenchmark)
bench <- microbenchmark(
  df_version = bulk_counts_summary_df(...),
  dt_version = bulk_counts_summary_dt(...),
  times = 10
)
autoplot(bench)
Across all workflows, data.table implementations demonstrate 1–2 orders of magnitude faster execution, particularly for joins and grouped aggregations.

Workflow Coverage
Task	Topic	Function
1	Bulk RNA-seq summary	bulk_counts_summary_dt()
2	QC columns	bulk_counts_qc_dt()
3	Subset + indexing	subset_counts_dt()
4	Annotate & summarize	annotate_counts_dt()
5	Clinical lab classification	classify_labs_dt()
6	Match labs ↔ vitals	match_vitals_dt()
7	Top ATAC peaks	top_peaks_dt()
8	Gene stats + filtering	gene_stats_filter_dt()
9	Reshape wide ↔ long	wide_long_wide_dt()
10	ATAC peaks ↔ genes	atac_to_gene_dt()
11	Variants ↔ genes	variants_to_genes_dt()
12	Multi-cohort merge	combine_cohorts_dt()
Final	Single-cell integration	final_revision_dt()

Folder Structure
ProjectAlessandra/
├── DESCRIPTION               # Package metadata
├── NAMESPACE                 # Function exports
├── R/                        # Core R source functions
├── man/                      # Auto-generated documentation (.Rd)
├── scripts/                  # Supplementary analysis or development scripts
├── vignettes/                # Detailed workflows and tutorials (.Rmd)
├── project_oct25/            # Project data and resources
├── Outputs/                  # Generated CSVs and plots (e.g. single-cell summaries)
├── README.md                 # You are here
└── rstudio.Rproj             # RStudio project file
Final Revision: Single-cell integration
The function final_revision_dt() merges Seurat integration clusters with cell type and tissue annotations, producing:
- Cluster–cell type–tissue summary tables
- Normalized frequency matrices
- Bar plots of cell type composition (normal vs tumor)

This closes the workflow from bulk and clinical omics to single-cell interpretation.

Example outputs are available in Outputs/:
- counts_per_cluster_celltype.csv
- normalized_percent_cluster_celltype_tissue.csv
- plot4_1.jpeg

Dependencies
Package	Purpose
data.table	Fast data manipulation
ggplot2	Visualization
microbenchmark	Performance comparison
knitr, rmarkdown	Reporting
dplyr (optional)	Readability in comparisons

📖 Citation
If you use this package or adapt its structure, please cite:

Alessandra, A. (2025). ProjectAlessandra: modular R workflows for multi-omics analysis using data.table.
GitHub Repository: https://github.com/your-username/ProjectAlessandra

💡 Author
Alessandra [Monterosso]
[alessandramonte050@gmail.com]

📜 License
This project is released under the MIT License.
