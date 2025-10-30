#Goal: Multi-column operations per group.
#Data: bulk_counts_long.csv, sample_metadata.csv 
#Tasks:
#• Compute per-condition robust summary stats for each gene: mean, median, Q1/Q3.
#• Return only genes where treated mean ‚â• 2√ó control mean.

library(data.table)
counts <- fread("project_oct25/bulk_counts_long.csv")
meta   <- fread("project_oct25/sample_metadata.csv")

#unione di counts e metadata per avere la colonna 'condition' insieme ai conteggi
merged <- counts[meta, on = "sample_id"]

#Calcoliamo le statistiche per gene e condition: mean, median, Q1 (25%) e Q3 (75%)
stats_by_gene_condition <- merged[, .(
  mean_count   = mean(count),
  median_count = median(count),
  Q1           = quantile(count, 0.25, type = 2), 
  Q3           = quantile(count, 0.75, type = 2)
), by = .(gene, condition)]

#Per applicare la regola "treated mean >= 2 * control mean" serve una tabella wide con le medie:
#separiamo le medie per condition e poi facciamo un merge per gene.
treated_means <- stats_by_gene_condition[condition == "treated", .(gene, treated_mean = mean_count)]
control_means <- stats_by_gene_condition[condition == "control", .(gene, control_mean = mean_count)]

#Uniamo le due tabelle delle medie per avere una riga per gene con entrambe le medie
means_wide <- merge(treated_means, control_means, by = "gene", all = FALSE)
# all = FALSE -> keep only genes present in both (necessario per confronto)

#Applichiamo il filtro: kept genes dove treated_mean >= 2 * control_mean
kept_genes <- means_wide[treated_mean >= 2 * control_mean]
