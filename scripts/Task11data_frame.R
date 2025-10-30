# =========================================================
# TASK 11 (versione base R con data.frame)
# Goal:
# - Mappare varianti (SNP) ai geni
# - Contare varianti ad alto impatto (HIGH) per gene e per sample
# =========================================================

# 1️⃣ Carico i dati -----------------------------------------------------------
variants <- read.csv("project_oct25/variants.csv")
genes    <- read.csv("project_oct25/gene_annotation.bed.csv")

# 2️⃣ Creo intervalli 1-bp per le varianti -----------------------------------
# Ogni variante è un punto, quindi creo due colonne identiche (start, end)
variants$start <- variants$pos
variants$end   <- variants$pos

# 3 Trovo le varianti che si sovrappongono a ciascun gene ------------------
# (equivalente a foverlaps)
# Facciamo un ciclo semplice: per ogni variante controlliamo se rientra nell’intervallo del gene
# Questo è molto meno efficiente, ma chiaro e comprensibile
overlaps_list <- list()

for (i in seq_len(nrow(variants))) {
  var_chr   <- variants$chr[i]
  var_start <- variants$start[i]
  var_end   <- variants$end[i]
  
  # Trova i geni sullo stesso cromosoma che si sovrappongono a quella variante?????
  #???????OVERLAPPING_GENES 0 OBS
  overlapping_genes <- subset(genes, 
                              chr == var_chr & 
                                start <= var_end & 
                                end >= var_start)
  
  # Se ci sono geni che si sovrappongono, li aggiungo alla lista
  if (nrow(overlapping_genes) > 0) {
    tmp <- data.frame(
      sample_id = variants$sample_id[i],
      gene = overlapping_genes$gene,
      chr = var_chr,
      pos = var_start,
      impact = variants$impact[i],
      stringsAsFactors = FALSE
    )
    overlaps_list[[length(overlaps_list) + 1]] <- tmp
  }
}

# Combino tutti i risultati in un unico data.frame
if (length(overlaps_list) > 0) {
  overlaps <- do.call(rbind, overlaps_list)
} else {
  overlaps <- data.frame()
}

# 4️⃣ Normalizzo la colonna impact a maiuscolo --------------------------------
if (nrow(overlaps) > 0) {
  overlaps$impact_upper <- toupper(overlaps$impact)
}

# 5️⃣ Filtro solo le varianti di alto impatto ---------------------------------
high_overlaps <- subset(overlaps, impact_upper == "HIGH")

# 6️⃣ Conteggio delle varianti HIGH per gene e per sample ---------------------
if (nrow(high_overlaps) > 0) {
  high_counts_by_gene_sample <- aggregate(
    x = list(high_variant_count = high_overlaps$impact_upper),
    by = list(gene = high_overlaps$gene, sample_id = high_overlaps$sample_id),
    FUN = length
  )
} else {
  high_counts_by_gene_sample <- data.frame()
}

# 7️⃣ Conteggio totale delle varianti HIGH per gene --------------------------
if (nrow(high_overlaps) > 0) {
  high_counts_by_gene <- aggregate(
    x = list(total_high_variants = high_overlaps$impact_upper),
    by = list(gene = high_overlaps$gene),
    FUN = length
  )
  
  # Ordina in modo decrescente
  high_counts_by_gene <- high_counts_by_gene[order(-high_counts_by_gene$total_high_variants), ]
  
  # Geni con almeno una variante HIGH
  genes_with_high <- unique(high_counts_by_gene$gene)
} else {
  high_counts_by_gene <- data.frame()
  genes_with_high <- character(0)
}

# 8️⃣ Salvo i risultati su file CSV ------------------------------------------
write.csv(high_counts_by_gene_sample, "project_oct25/high_variants_by_gene_sample.csv", row.names = FALSE)
write.csv(high_counts_by_gene, "project_oct25/high_variants_by_gene_total.csv", row.names = FALSE)
write.csv(data.frame(gene = genes_with_high), "project_oct25/genes_with_high_variants.csv", row.names = FALSE)

cat("\nAnalisi completata. File salvati in project_oct25/.\n")

