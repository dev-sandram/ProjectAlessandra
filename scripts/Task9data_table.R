#Goal: Go wide‚ to long‚ to wide for downstream plotting. Data: bulk_counts_wide.csv
#Tasks:
# • Convert the matrix to long, add per-sample totals, then back to a gene per condition table with mean counts.

# =====================================================
# TASK 9: Go wide → long → wide for downstream plotting
# =====================================================

# Carichiamo le librerie
library(data.table)

# Leggiamo i dati (matrice larga: una riga per gene, una colonna per campione)
counts_wide <- fread("project_oct25/bulk_counts_wide.csv")


# PARTE 1: Convertiamo da formato wide → long
# -----------------------------------------------------
# Supponiamo che la prima colonna si chiami "gene"
counts_long <- melt(counts_wide, id.vars = "gene",
                    variable.name = "sample_id", value.name = "count")
#melt()	scioglie la tabella da larga a lunga.
#counts_wide	è la tabella di partenza (quella larga).
#id.vars = "gene"	dice a R: mantieni la colonna “gene” così com’è (non la trasformare).
#Tutte le altre colonne (sample_1, sample_2, ecc.) verranno “sciolte”.
#variable.name = "sample_id"	dà un nome alla nuova colonna che conterrà i nomi delle vecchie colonne (cioè i nomi dei campioni)
#value.name = "count"	dà un nome alla colonna che conterrà i valori numerici (cioè i conteggi per gene e campione).

# Aggiungiamo le informazioni di condizione per ogni sample
# Carichiamo i metadati
meta <- fread("project_oct25/sample_metadata.csv")

# Uniamo metadati ai conteggi
merged <- merge(counts_long, meta, by = "sample_id")

# Calcoliamo i totali per campione
totals_per_sample <- merged[, .(total_count = sum(count)), by = sample_id]

# Aggiungiamo questi totali alla tabella principale
merged <- merge(merged, totals_per_sample, by = "sample_id")

# Torniamo da long → wide,
#           ma ora con colonne per condizione (treated, control)
#           e valori = media dei conteggi per gene
gene_condition_means <- merged[, .(mean_count = mean(count)),
                               by = .(gene, condition)]

# Da long a wide: una riga per gene, colonne = condizioni
counts_condition_wide <- dcast(gene_condition_means,
                               gene ~ condition,
                               value.var = "mean_count")

