# ==========================================================
# TASK 8 – data.frame version
# ==========================================================

counts <- read.csv("project_oct25/bulk_counts_long.csv", stringsAsFactors = FALSE)
meta   <- read.csv("project_oct25/sample_metadata.csv", stringsAsFactors = FALSE)

merged <- merge(counts, meta, by = "sample_id")

stats_by_gene_condition <- aggregate(count ~ gene + condition, data = merged,
                                     FUN = function(x) c(mean = mean(x), median = median(x),
                                                         Q1 = quantile(x, 0.25, type = 2),
                                                         Q3 = quantile(x, 0.75, type = 2)))

stats_df <- data.frame(
  gene = stats_by_gene_condition$gene,
  condition = stats_by_gene_condition$condition,
  mean_count = stats_by_gene_condition$count[, "mean"],
  median_count = stats_by_gene_condition$count[, "median"],
  Q1 = stats_by_gene_condition$count[, "Q1.25%"],
  Q3 = stats_by_gene_condition$count[, "Q3.75%"]
)

treated <- subset(stats_df, condition == "treated")[, c("gene", "mean_count")]
control <- subset(stats_df, condition == "control")[, c("gene", "mean_count")]
names(treated)[2] <- "treated_mean"
names(control)[2] <- "control_mean"

means_wide <- merge(treated, control, by = "gene")
kept_genes <- subset(means_wide, treated_mean >= 2 * control_mean)

