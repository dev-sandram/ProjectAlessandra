#Final revision
#Goal: Associate cell type to integration clusters taking tract of their belongings to normal or tumor tissue.
#Data:
#annotated_GSM3516673_normal_annotated_GSM3516672_tumor_SeuratIntegration.csv contains the following columns:
#• cell = cell id
#• integration_cluster = integration cluster

#nt_combined_clustering.output.csv contains the following columns:
#• cell = cell id
#• cell_type = predicted cell type
#• sample_type = N for normal tissue neat to tumor, and T for tumor tissue

#SETUP PACCHETTI E PERCORSI
# Setup iniziale
library(data.table)
library(ggplot2)

# ==========================================
# LETTURA FILE DI INPUT
# ==========================================
integration_dt <- fread("project_oct25/annotated_GSM3516673_normal_annotated_GSM3516672_tumor_SeuratIntegration.csv")
clustering_dt  <- fread("project_oct25/nt_combined_clustering.output.csv")

# ==========================================
# FUNZIONE NORMALIZZAZIONE ID CELLULARI (necessaria)
#perchè i nomi sono diversi: rimuove anche _X_ o _Y_ in mezzo agli ID
# ==========================================
normalize_cell_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)            # toglie spazi iniziali/finali
  x <- gsub("_X_", "_", x)  # rimuove sequenze _X_ nel mezzo
  x <- gsub("_Y_", "_", x)  # rimuove sequenze _Y_ nel mezzo
  x <- gsub("_\\.", ".", x) # rimuove underscore prima del punto
  return(x)
}

# Applichiamo la funzione di normalizzazione
integration_dt[, cell_clean := normalize_cell_id(cell)]
clustering_dt[, cell_clean := normalize_cell_id(cell)]

# ==========================================
# TASK 1.1
# ==========================================
#provide a new file where cell type, cells and integration clusters are combined
#Unione (join) sui cell id
combined <- merge(integration_dt, clustering_dt, by = "cell_clean", all = FALSE)  #all FALSE evita duplicazioni indesiderate

fwrite(combined,
       file = "Outputs/combined_celltype_integrationcluster.csv")

# ==========================================
# TASK 2.1
# ==========================================
#provide a file where the number of each cell type is indicated for each cluster
#Per ogni integration_cluster contare quante cells di ciascun cell_type ci sono
#(tabella cluster × cell_type con conti).
counts_cluster_celltype <- combined[, .N, by = .(integration_cluster, cell_type)]
setnames(counts_cluster_celltype, "N", "cell_count")

fwrite(counts_cluster_celltype,
       file = "Outputs/celltype_counts_per_cluster.csv")

# ==========================================
# TASK 3.1
# ==========================================
#provide summary table showing cell types associated to integration clusters
#and also their association to normal and tumor tissue
#Costruire una tabella che per ogni (integration_cluster, cell_type, sample_type) riporta:
#cell_count (numero di celle),
#pct_within_cluster = percentuale di quel cell_type dentro il cluster
#(cioè cell_count / total_cells_in_cluster * 100),
#pct_within_celltype = percentuale di quelle celle di quel tipo che sono N o T
#(utile per vedere normal vs tumor per quello specifico cell_type nel cluster).
# Conta celle per cluster × cell_type × sample_type

# 2. Conta righe per (cluster, cell_type, sample_type)
tab_ct_st <- combined[, .N, by = .(integration_cluster, cell_type, sample_type)]
setnames(tab_ct_st, "N", "cell_count")

# Totale cells per cluster, calcola totale per cluster (per normalizzare dentro il cluster)
totals_cluster <- combined[, .N, by = integration_cluster]
setnames(totals_cluster, "N", "cluster_total_cells")

# Unisci il totale per cluster
tab_ct_st <- merge(tab_ct_st, totals_cluster, by = "integration_cluster")

# Percentuale all'interno del cluster
tab_ct_st[, pct_within_cluster := (cell_count / cluster_total_cells) * 100]

# Totale celle per cluster × cell_type, calcolo percentuale dentro la coppia (cluster, cell_type) split by sample_type
#    Prima calcolo il totale per (cluster, cell_type) (sommando tutte le sample_type)
totals_cluster_celltype <- combined[, .N, by = .(integration_cluster, cell_type)]
setnames(totals_cluster_celltype, "N", "cluster_celltype_total")

tab_ct_st <- merge(tab_ct_st, totals_cluster_celltype,
                   by = c("integration_cluster", "cell_type"))

# Percentuale di sample_type all'interno di quel tipo cellulare nel cluster
tab_ct_st[, pct_within_cluster_celltype := (cell_count / cluster_celltype_total) * 100]

fwrite(tab_ct_st,
       file = "Outputs/summary_cluster_celltype_tissue.csv")


# ==========================================
# TASK 4.1
# ==========================================
#Provide a plot describing the distribution of the cell type in normal/tumor tissue
#given the integration clusters.

# Usiamo il file che hai già creato nel Task 1
# (quello con le colonne: cell_type, integration_cluster, sample_type)
plot_data <- as.data.frame(table(combined$integration_cluster,
                                 combined$cell_type,
                                 combined$sample_type))
names(plot_data) <- c("cluster", "cell_type", "tissue", "count")

totals <- aggregate(count ~ cluster + tissue, data = plot_data, sum)
names(totals)[3] <- "total"

plot_data <- merge(plot_data, totals, by = c("cluster", "tissue"))
plot_data$percent <- (plot_data$count / plot_data$total) * 100

# Facciamo il grafico
ggplot(plot_data, aes(x = cluster, y = percent, fill = cell_type)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ tissue) +
  labs(title = "Distribution of cell types in clusters",
       x = "Integration cluster",
       y = "Cells percentage (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ==========================================
# TASK 5.1
# ==========================================
#Provide a normalized % for the cell types in each of the integration clusters
#for normal and tumor specimen.
counts_norm <- combined[, .N, by = .(cluster = integration_cluster,
                                     cell_type,
                                     tissue = sample_type)]

totals_norm <- counts_norm[, .(total = sum(N)), by = .(cluster, tissue)]

counts_norm <- merge(counts_norm, totals_norm, by = c("cluster", "tissue"))
counts_norm[, percent := round((N / total) * 100, 2)]

fwrite(counts_norm,
       file = "Outputs/normalized_celltype_percentages_by_tissue.csv")

