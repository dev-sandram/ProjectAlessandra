peaks <- read.csv("project_oct25/atac_peaks.bed.csv", stringsAsFactors = FALSE)
genes <- read.csv("project_oct25/gene_annotation.bed.csv", stringsAsFactors = FALSE)

# overlap calcolato “a mano”
overlaps_list <- lapply(1:nrow(peaks), function(i) {
  p <- peaks[i, ]
  gsub <- subset(genes, chr == p$chr & start <= p$end & end >= p$start)
  if (nrow(gsub) == 0) return(NULL)
  gsub$overlap_bp <- pmin(p$end, gsub$end) - pmax(p$start, gsub$start)
  gsub <- gsub[gsub$overlap_bp > 0, ]
  gsub$peak_id <- i
  gsub
})
overlaps <- do.call(rbind, overlaps_list)

peaks_per_gene <- aggregate(overlap_bp ~ gene, data = overlaps, FUN = length)
names(peaks_per_gene)[2] <- "num_peaks"

overlap_sum_per_gene <- aggregate(overlap_bp ~ gene, data = overlaps, sum)
names(overlap_sum_per_gene)[2] <- "total_overlap_bp"
top20_genes <- head(overlap_sum_per_gene[order(-overlap_sum_per_gene$total_overlap_bp), ], 20)
