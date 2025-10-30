#Goal: Annotate counts with sample and patient info. Data: bulk_counts_long.csv, sample_metadata.csv Tasks:
#• Join and compute per-patient total counts.
#• Find the top 10 genes by average count within each condition.

library(data.table)

# 1. Import
bulk_counts <- fread("project_oct25/bulk_counts_long.csv")
sample_meta   <- fread("project_oct25/sample_metadata.csv")

# 2. Join counts + metadata
setkey(bulk_counts, sample_id)
setkey(sample_meta, sample_id)
join_data <- bulk_counts[sample_meta, nomatch = 0]

# 3. Total counts per patient
patient_tot <- join_data[, .(total_count = sum(count)), by = patient_id]

# 4. Mean count per gene and condition
gene_means <- join_data[, .(mean_count = mean(count)), by = .(gene, condition)]

# 5. Top 10 genes (highest mean) within each condition
top10 <- gene_means[
  order(condition, -mean_count)
][, head(.SD, 10), by = condition]

#by = condition → significa: “ripeti il comando successivo separatamente per ogni valore di condition”.
#.SD → significa “Subset of Data”: è la porzione della tabella relativa al gruppo attuale
#head(.SD, 10) → prende le prime 10 righe del gruppo (cioè i 10 geni con i valori più alti di mean_count).

result <- annotate_counts_dt("project_oct25/bulk_counts_long.csv","project_oct25/sample_metadata.csv")
patient_tot <- result$patient_tot
top10 <- result$top10
