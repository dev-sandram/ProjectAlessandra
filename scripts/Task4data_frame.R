# ==========================================================
# TASK 4 – data.frame version
# ==========================================================

# 1. Import
counts <- read.csv("project_oct25/bulk_counts_long.csv", stringsAsFactors = FALSE)
meta   <- read.csv("project_oct25/sample_metadata.csv", stringsAsFactors = FALSE)

# 2. Join
merged_data <- merge(counts, meta, by = "sample_id")

# 3. Total counts per patient
patient_totals <- aggregate(count ~ patient_id, data = merged_data, FUN = sum)
names(patient_totals)[2] <- "total_count"

# 4. Mean count per gene and condition
gene_means <- aggregate(count ~ gene + condition, data = merged_data, FUN = mean)
names(gene_means)[3] <- "mean_count"

# 5. Top 10 genes by condition
conditions <- unique(gene_means$condition)
top10_by_condition <- data.frame()

for (cond in conditions) {
  subset_cond <- subset(gene_means, condition == cond)
  subset_cond <- subset_cond[order(-subset_cond$mean_count), ]
  top10 <- head(subset_cond, 10)
  top10_by_condition <- rbind(top10_by_condition, top10)
}


