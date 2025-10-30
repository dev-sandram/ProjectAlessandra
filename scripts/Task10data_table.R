#Goal: ATAC-to-gene mapping.
#Data: atac_peaks.bed.csv, gene_annotation.bed.csv 
#Tasks:
#• setkey() both on (chr, start, end) and intersect peaks with gene bodies.
#• Count peaks per gene.
#• Compute overlap length (bp) per peak-gene pair and then sum per gene.
#• Return the top 20 genes by total overlapped bp.

# Carichiamo la libreria
library(data.table)

# Leggiamo i file
peaks <- fread("project_oct25/atac_peaks.bed.csv")
genes <- fread("project_oct25/gene_annotation.bed.csv")

# Imposta le chiavi per (chr, start, end)
setkey(peaks, chr, start, end)
setkey(genes, chr, start, end)

# Trova intersezioni tra picchi e geni
# (non-equi join: regioni che si sovrappongono)
# La condizione di overlap: 
# start_peak <= end_gene e end_peak >= start_gene
overlaps <- foverlaps(peaks, genes, 
                      by.x = c("chr", "start", "end"), 
                      by.y = c("chr", "start", "end"), 
                      type = "any", nomatch = 0L)

# Calcola la lunghezza di overlap (in basi)
# L'overlap fra due intervalli è: min(end1, end2) - max(start1, start2)
overlaps[, overlap_bp := pmin(end, i.end) - pmax(start, i.start)] 

# Mantieni solo le righe con overlap positivo
overlaps <- overlaps[overlap_bp > 0]

# Conta quanti picchi per gene
peaks_per_gene <- overlaps[, .N, by = gene]
setnames(peaks_per_gene, "N", "num_peaks")

# Somma la lunghezza totale di overlap per gene
overlap_sum_per_gene <- overlaps[, .(total_overlap_bp = sum(overlap_bp)), by = gene]

# Ordina per overlap totale (decrescente) e prendi top 20
# -----------------------------------------------------
top20_genes <- overlap_sum_per_gene[order(-total_overlap_bp)][1:20]
