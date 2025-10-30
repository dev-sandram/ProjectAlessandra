#Goal: Combine cohorts safely.
#Data: cohortA_samples.csv, cohortB_samples.csv
#Tasks:
#• rbindlist(list(A, B), use.names=TRUE, fill=TRUE) and verify column alignment.
#• Order by cohort, condition, sample_id.
#• Join back to bulk_counts_long and compute per-cohort, per-condition mean counts of the top 100 most variable genes.

# Carichiamo la libreria
library(data.table)

# Leggiamo i file, specificando i percorsi
cohortA <- fread("project_oct25/cohortA_samples.csv")
cohortB <- fread("project_oct25/cohortB_samples.csv")
counts  <- fread("project_oct25/bulk_counts_long.csv")

# Aggiungiamo una colonna 'cohort' per tenere traccia della provenienza
cohortA[, cohort := "A"]
cohortB[, cohort := "B"]

# Uniamo le due coorti
# use.names = TRUE --> allinea le colonne con lo stesso nome
# fill = TRUE      --> se mancano colonne in uno dei due file, le crea e le riempie con NA
combined_cohorts <- rbindlist(list(cohortA, cohortB), use.names = TRUE, fill = TRUE)

# Ordiniamo per coorte, condizione e sample_id
setorder(combined_cohorts, cohort, condition, sample_id)

# Uniamo per sample_id
merged_per_sampleid <- merge(counts, combined_cohorts, by = "sample_id", all.x = TRUE)

# PARTE 3: Troviamo i 100 geni più variabili
# Calcoliamo la varianza dei conteggi per ciascun gene
gene_variance <- merged_per_sampleid[, .(variance = var(count, na.rm = TRUE)), by = gene]

# Ordiniamo per varianza decrescente e prendiamo i primi 100
top100_genes <- gene_variance[order(-variance)][1:100, gene]
print(head(top100_genes, 5))

# Calcoliamo la media dei conteggi per coorte e condizione
#           solo per i top 100 geni
top100_data <- merged_per_sampleid[gene %in% top100_genes]

mean_counts <- top100_data[, .(
  mean_count = mean(count, na.rm = TRUE)
), by = .(gene, cohort, condition)]
