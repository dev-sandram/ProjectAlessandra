# ==========================================================
# TASK 2 – data.frame version
# ==========================================================

# 1. Import
counts <- read.csv("project_oct25/bulk_counts_long.csv", stringsAsFactors = FALSE)

# 2. Add log2(count + 1) column
counts$log2_count <- log2(counts$count + 1)

# 3. Add binary flag 'high' (count > 100)
counts$high <- counts$count > 100

# 4. Overwrite 'high' using gene-wise median threshold
#    Qui usiamo tapply() per calcolare la mediana per gene,
#    poi combiniamo i risultati in un vettore logico della stessa lunghezza
medians_by_gene <- tapply(counts$count, counts$gene, median)
counts$high <- counts$count > medians_by_gene[counts$gene]
