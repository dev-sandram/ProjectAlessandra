# ==========================================================
# TASK 9 – data.frame version
# ==========================================================

counts_wide <- read.csv("project_oct25/bulk_counts_wide.csv", stringsAsFactors = FALSE)
counts_long <- reshape(counts_wide, varying = names(counts_wide)[-1],
                       v.names = "count", timevar = "sample_id",
                       times = names(counts_wide)[-1], idvar = "gene",
                       direction = "long")

meta <- read.csv("project_oct25/sample_metadata.csv", stringsAsFactors = FALSE)
merged <- merge(counts_long, meta, by = "sample_id")

totals_per_sample <- aggregate(count ~ sample_id, data = merged, sum)
names(totals_per_sample)[2] <- "total_count"

merged <- merge(merged, totals_per_sample, by = "sample_id")

gene_condition_means <- aggregate(count ~ gene + condition, data = merged, mean)
counts_condition_wide <- reshape(gene_condition_means,
                                 timevar = "condition", idvar = "gene",
                                 direction = "wide")

cat("\n--- Tabella finale: media dei conteggi per gene e condizione ---\n")
print(head(counts_condition_wide, 10))
