# Carichiamo solo i pacchetti necessari
library(ggplot2)

# Funzione di normalizzazione ID cellule
normalize_cell_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("_X_", "_", x)
  x <- gsub("_Y_", "_", x)
  x <- gsub("_\\.", ".", x)
  return(x)
}

# ==========================================
# TASK 1.1
# ==========================================

# Leggi i file in data.frame
integration_df <- read.csv("project_oct25/annotated_GSM3516673_normal_annotated_GSM3516672_tumor_SeuratIntegration.csv")
clustering_df  <- read.csv("project_oct25/nt_combined_clustering.output.csv")

# Normalizza gli ID
integration_df$cell_clean <- normalize_cell_id(integration_df$cell)
clustering_df$cell_clean  <- normalize_cell_id(clustering_df$cell)

# Join sulle cellule
combined_df <- merge(integration_df, clustering_df, by = "cell_clean", all = FALSE)

# Salva file Task 1
write.csv(combined_df, "Outputs/combined_cell_cluster_celltype.csv", row.names = FALSE)

# ==========================================
# TASK 2.1
# ==========================================

counts_cluster_celltype <- aggregate(cell_clean ~ integration_cluster + cell_type,
                                     data = combined_df, FUN = length)
colnames(counts_cluster_celltype)[3] <- "cell_count"

write.csv(counts_cluster_celltype, "Outputs/counts_per_cluster_celltype.csv", row.names = FALSE)

# ==========================================
# TASK 3.1
# ==========================================

# Conta per cluster, cell_type, sample_type
tab_ct_st <- aggregate(cell_clean ~ integration_cluster + cell_type + sample_type,
                       data = combined_df, FUN = length)
colnames(tab_ct_st)[4] <- "cell_count"

# Totale per cluster
totals_cluster <- aggregate(cell_clean ~ integration_cluster,
                            data = combined_df, FUN = length)
colnames(totals_cluster)[2] <- "cluster_total_cells"

# Merge per aggiungere totale cluster
tab_ct_st <- merge(tab_ct_st, totals_cluster, by = "integration_cluster", all.x = TRUE)

# Percentuale dentro cluster
tab_ct_st$pct_within_cluster <- (tab_ct_st$cell_count / tab_ct_st$cluster_total_cells) * 100

# Totale per (cluster, cell_type)
totals_cluster_celltype <- aggregate(cell_clean ~ integration_cluster + cell_type,
                                     data = combined_df, FUN = length)
colnames(totals_cluster_celltype)[3] <- "cluster_celltype_total"

# Merge + percentuale sample_type dentro cell_type del cluster
tab_ct_st <- merge(tab_ct_st, totals_cluster_celltype,
                   by = c("integration_cluster", "cell_type"), all.x = TRUE)

tab_ct_st$pct_within_cluster_celltype <- (tab_ct_st$cell_count / tab_ct_st$cluster_celltype_total) * 100

# Salva Task 3
write.csv(tab_ct_st, "Outputs/summary_cluster_celltype_sampletype.csv", row.names = FALSE)

# ==========================================
# TASK 4.1
# ==========================================

plot_data <- as.data.frame(table(combined_df$integration_cluster,
                                 combined_df$cell_type,
                                 combined_df$sample_type))
colnames(plot_data) <- c("cluster", "cell_type", "tissue", "count")

# Totale per cluster e tessuto
totals <- aggregate(count ~ cluster + tissue, data = plot_data, sum)
colnames(totals)[3] <- "total"

# Merge + percentuali
plot_data <- merge(plot_data, totals, by = c("cluster", "tissue"))
plot_data$percent <- (plot_data$count / plot_data$total) * 100

# Plot
ggplot(plot_data, aes(x = cluster, y = percent, fill = cell_type)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ tissue) +
  labs(title = "Distribution of cell types in clusters (Normal vs Tumor)",
       x = "Integration Cluster",
       y = "Cell Percentage (%)") +
  theme_minimal()


# ==========================================
# TASK 5.1
# ==========================================
# Conta celle
counts <- as.data.frame(table(combined_df$integration_cluster,
                              combined_df$cell_type,
                              combined_df$sample_type))
colnames(counts) = c("cluster", "cell_type", "tissue", "count")

# Totale per cluster e tessuto
totals <- aggregate(count ~ cluster + tissue, data = counts, sum)
colnames(totals)[3] <- "total"

# Merge e calcolo % normalizzata
counts_norm <- merge(counts, totals, by = c("cluster", "tissue"))
counts_norm$percent <- round((counts_norm$count / counts_norm$total) * 100, 2)

# Salva Task 5
write.csv(counts_norm, "Outputs/normalized_percent_cluster_celltype_tissue.csv", row.names = FALSE)
