cohortA <- read.csv("project_oct25/cohortA_samples.csv", stringsAsFactors = FALSE)
cohortB <- read.csv("project_oct25/cohortB_samples.csv", stringsAsFactors = FALSE)
counts  <- read.csv("project_oct25/bulk_counts_long.csv", stringsAsFactors = FALSE)

cohortA$cohort <- "A"
cohortB$cohort <- "B"

combined_cohorts <- rbind(cohortA, cohortB)
combined_cohorts <- combined_cohorts[order(combined_cohorts$cohort,
                                           combined_cohorts$condition,
                                           combined_cohorts$sample_id), ]

merged <- merge(counts, combined_cohorts, by = "sample_id", all.x = TRUE)

gene_variance <- aggregate(count ~ gene, data = merged, FUN = var, na.rm = TRUE)
names(gene_variance)[2] <- "variance"
top100_genes <- head(gene_variance[order(-gene_variance$variance), "gene"], 100)

top100_data <- subset(merged, gene %in% top100_genes)
mean_counts <- aggregate(count ~ gene + cohort + condition,
                         data = top100_data, FUN = mean, na.rm = TRUE)
names(mean_counts)[4] <- "mean_count"

