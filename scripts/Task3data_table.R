#Goal: Speed up joins/lookups.
#Data: sample_metadata.csv, bulk_counts_long.csv
#Tasks:
#• setkey() on sample_metadata by sample_id; equi-join into the long counts.
#• Add a secondary index on (gene, sample_id) in the counts table, then benchmark a subset query before/after.

# 1 Carico la libreria
library(data.table)

# 2 Leggo i file CSV come data.table
bulk_counts <- fread("project_oct25/bulk_counts_long.csv")
sample_meta   <- fread("project_oct25/sample_metadata.csv")

# PARTE 1: Imposto la chiave su sample_metadata
# setkey() serve per ordinare la tabella e prepararla a join più veloci
setkey(sample_meta, sample_id)

# PARTE 2: Faccio una join tra metadata e counts
# Dopo aver impostato la chiave, la join è automatica e più veloce
join_data <- sample_meta[bulk_counts, on = "sample_id"] #prende meta e attacca counts in base al sample ID

# PARTE 4: Benchmark — confronto prima e dopo l’indice
# Esempio: voglio estrarre tutte le righe per un certo gene
gene_call <- "GENE_0051"
sample_chosen <- "S20"

# Prima (senza usare l’indice)
system.time({
  subset_no_index <- bulk_counts[gene == gene_call & sample_id == sample_chosen]
})

# PARTE 3: Creo un indice secondario su (gene, sample_id)
# L’indice serve per accelerare query che filtrano per gene e sample
setindex(bulk_counts, gene, sample_id)

# Dopo (l’indice ora è attivo)
system.time({
  subset_with_index <- bulk_counts[gene == gene_call & sample_id == sample_chosen]
})


result <- subset_counts_dt("project_oct25/bulk_counts_long.csv","project_oct25/sample_metadata.csv","GENE_0051","S20")
