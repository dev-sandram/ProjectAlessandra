#Goal: Filter, summarize, and group bulk counts. Data: bulk_counts_long.csv, sample_metadata.csv Tasks:
#• Keep counts for samples in condition == "treated" and genes starting with "GENE_00".
#• Compute mean and median count by gene (valore centrale).
#• Join sample metadata and compute per-condition mean counts by gene in one pipeline.

library(data.table)

# 1. Import
bulk_counts <- fread("project_oct25/bulk_counts_long.csv")
sample_meta   <- fread("project_oct25/sample_metadata.csv")

# 2. Join counts + metadata by sample_id
#imposta sample_id come chiave sull'oggetto counts e meta quindi counts viene ordinato per sample_id
#effettua una modifica in place (non crea una copia)
#questo permette join molto più rapidi
setkey(bulk_counts, sample_id) #setkey ordina la tabella e rende più veloce l'unione di tabelle
setkey(sample_meta, sample_id)
join_data <- bulk_counts[sample_meta, nomatch = 0]
#In data.table la forma X[i] può essere usata per i join:
#X = tabella da cui prendi i dati (qui counts)
#i = tabella che serve come query/indice (qui meta)
#Quando ci sono chiavi compatibili (sample_id), counts[meta] interpreta meta come la tabella di lookup
#e per ciascuna riga di meta estrae le righe corrispondenti da counts
#nomatch = 0 significa: escludi le righe di meta che non trovano corrispondenza in counts.
#Questo equivale a un inner join (solo righe with match).

# 3. Filter for treated samples and GENE_00* genes
filtered_data <- join_data[condition == "treated" & grepl("^GENE_00", gene)]
#filtrare stringhe in base a un pattern (cioè una “regola” di testo).
#Quindi grepl("^GENE_00", gene) restituisce:
#TRUE per tutti i geni che iniziano con GENE_00
#FALSE per tutti gli altri.

# 4. Compute mean and median per gene
gene_mean_median <- filtered_data[, .(
  mean_count   = mean(count),
  median_count = median(count)
), by = gene]

# Punto 5 corretto: join + per-condition mean by gene (una pipeline)

# Calcola la media dei conteggi per ogni combinazione (gene, condition)
# tutto in una singola pipeline: unisco metadata e counts e poi raggruppo
gene_condition_means <- bulk_counts[
  sample_meta,                # uso sample_meta come "i" (lookup)
  on = "sample_id"            # join su sample_id (non serve setkey se usi on=)
][ , .(
 mean_count = mean(count)   # media dei conteggi (ignora NA se ci sono)
), by = .(gene, condition) ]               # raggruppa per gene e per condition


#CON FUNCTION:
result <- bulk_counts_summary_dt("project_oct25/bulk_counts_long.csv","project_oct25/sample_metadata.csv")
filtered_data <- result$filtered_data
gene_mean_median <- result$gene_mean_median
gene_condition_means <- result$gene_condition_means




