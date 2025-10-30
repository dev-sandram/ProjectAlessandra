# ==========================================================
# TASK 1 – data.frame version
# Stesso obiettivo ma usando solo funzioni base R
# ==========================================================

# 1. Import
counts <- read.csv("project_oct25/bulk_counts_long.csv")
meta   <- read.csv("project_oct25/sample_metadata.csv")

# 2. Join counts + metadata by sample_id
merged_data <- merge(counts, meta, by = "sample_id")

# 3. Filter for treated samples and GENE_00* genes
treated_data <- subset(merged_data,
                       condition == "treated" & grepl("^GENE_00", gene))

# 4. Compute mean and median per gene
gene_summary <- aggregate(
  count ~ gene,
  data = treated_data,
  FUN  = function(x) c(mean = mean(x), median = median(x))
)

# 5. Simplify format
gene_summary_df <- data.frame(
  gene = gene_summary$gene,
  mean_count   = gene_summary$count[, "mean"],
  median_count = gene_summary$count[, "median"]
)

# 6 Calcolo della media dei conteggi per ciascun gene e condizione
# aggregate() calcola una funzione (mean) per gruppi
gene_condition_means_f <- aggregate(
  count ~ gene + condition,   # formula: raggruppa per gene e condition
  data = merged_data,         # tabella di partenza
  FUN = mean,                 # funzione da applicare
  na.rm = TRUE                # ignora eventuali valori mancanti
)

# 7 Ordina i risultati per gene e condition (per leggibilità)
gene_condition_means_f <- gene_condition_means_f[order(gene_condition_means_f$gene,
                                                   gene_condition_means_f$condition), ]

