# ==========================================================
# 🧩 TASK 1 — Bulk RNA counts summary
# ==========================================================

#' Summarize bulk RNA counts with metadata
#'
#' This function reads in raw bulk RNA-seq counts and metadata,
#' it merges bulk RNA-seq count data with sample metadata by "sample_id",
#' filters treated samples whose gene names start with "GENE_00",
#' and computes multiple summary statistics:
#' mean/median counts per gene as well as per-condition mean counts.
#' Designed for quick quality control and exploratory summaries using data.table's fast aggregation capabilities.
#'
#' @param counts_path Path to CSV file with columns: gene, sample_id, count.
#' @param meta_path Path to CSV file with columns: sample_id, condition, batch, patient ID, timepoint.
#' @return List with: filtered_data, gene_mean_median, gene_condition_means
#' \itemize{
#'   \item gene_mean_median: Mean and median counts by gene (treated only)
#'   \item gene_condition_means: Mean counts by gene and condition
#' }
#' @import data.table
#' @export
#'
bulk_counts_summary_dt <- function(counts_path, meta_path) {
  bulk_counts <- fread(counts_path)
  sample_meta <- fread(meta_path)

  setkey(bulk_counts, sample_id)
  setkey(sample_meta, sample_id)

  join_data <- bulk_counts[sample_meta, on = "sample_id", nomatch = 0]

  filtered_data <- join_data[condition == "treated" & grepl("^GENE_00", gene)]

  gene_mean_median <- filtered_data[, .(
    mean_count   = mean(count),
    median_count = median(count)
  ), by = gene]

  gene_condition_means <- bulk_counts[
    sample_meta,
    on = "sample_id"
  ][ , .(
    mean_count = mean(count)
  ), by = .(gene, condition) ]

  return(list(filtered_data = filtered_data, gene_mean_median = gene_mean_median, gene_condition_means = gene_condition_means))
}

# ==========================================================
# 🧩 TASK 2 — QC-style derived columns
# ==========================================================

#' Add QC-style derived columns
#'
#' Computes log2-transformed counts and adds a binary \code{high} flag
#' indicating whether a count is above the gene-wise median.
#' This mimics typical preprocessing steps in QC pipelines to visualize
#' count distributions and detect outliers.
#'
#' @param counts_path Path to CSV file with columns: gene, sample_id, count.
#' @return Data.table with added columns: log2_count, high.
#' @import data.table
#' @export
bulk_counts_qc_dt <- function(counts_path) {
  counts <- fread(counts_path)
  counts[, log2_count := log2(count)]
  counts[, high := count > 100]
  counts[, high := count > median(count), by = gene]

  return(counts)
}

# ==========================================================
# 🧩 TASK 3 — Subset by gene and sample
# ==========================================================

#' Subset counts data using secondary index
#'
#' Joins metadata and benchmarks subsetting by gene and sample with
#' and without secondary indexing — useful for optimizing workflows
#' with large RNA-seq tables.
#'
#' @param counts_path Path to CSV file with columns: gene, sample_id, count.
#' @param meta_path Path to CSV file with columns: sample_id, condition, batch, patient ID, timepoint.
#' @param gene_name Gene name to subset
#' @param sample_chosen Sample ID to subset
#' @return Subset of counts for the given gene and sample.
#' @import data.table
#' @export
subset_counts_dt <- function(counts_path, meta_path, gene_call, sample_chosen) {
  bulk_counts <- fread(counts_path)
  sample_meta <- fread(meta_path)

  setkey(sample_meta, sample_id)

  join_data <- sample_meta[bulk_counts, on = "sample_id"]

  time_no_index <- system.time({
    subset_no_index <- bulk_counts[gene == gene_call & sample_id == sample_chosen]
  })

  setindex(bulk_counts, gene, sample_id)

  time_with_index <- system.time({
    subset_with_index <- bulk_counts[gene == gene_call & sample_id == sample_chosen]
})

  dt_benchmark <- data.table(
    test = c("time_no_index", "time_with_index"),
    user = c(time_no_index["user.self"], time_with_index["user.self"]),
    system = c(time_no_index["sys.self"], time_with_index["sys.self"]),
    elapsed = c(time_no_index["elapsed"], time_with_index["elapsed"])
  )

 return(dt_benchmark)
  }

# ==========================================================
# 🧩 TASK 4 — Annotate counts and summarize
# ==========================================================

#' Annotate counts with metadata and compute summaries
#'
#' Adds patient metadata to bulk counts, computes per-patient total counts,
#' per-gene means by condition, and identifies top 10 genes per condition
#' ranked by average expression.
#'
#' @param counts_path Path to CSV file with columns: gene, sample_id, count.
#' @param meta_path Path to CSV file with columns: sample_id, condition, batch, patient ID, timepoint.
#' @return List with patient totals, and top10 genes.
#' @import data.table
#' @export
annotate_counts_dt <- function(counts_path, meta_path) {

  bulk_counts <- fread(counts_path)
  sample_meta <- fread(meta_path)


  setkey(bulk_counts, sample_id)
  setkey(sample_meta, sample_id)
  join_data <- bulk_counts[sample_meta, nomatch = 0]


  patient_tot <- join_data[, .(total_count = sum(count)), by = patient_id]


  gene_means <- join_data[, .(mean_count = mean(count)), by = .(gene, condition)]


  top10 <- gene_means[
    order(condition, -mean_count)
  ][, head(.SD, 10), by = condition]

  return(list(patient_tot = patient_tot, top10 = top10))
}

# ==========================================================
# 🧩 TASK 5 — Classify lab values
# ==========================================================

#' Classify lab values against reference intervals
#'
#' Merges patient lab results with reference intervals and assigns each test
#' as "normal" or "out_of_range", and summarizes abnormalities per patient and per lab.
#'
#' @param labs_path Path to CSV file with columns: patient_id, time_iso, lab, value.
#' @param ref_path Path to CSV file with columns: lab, sex,lower, upper.
#' @return List with merged_labs, abnormal_by_patient, abnormal_by_lab.
#' @import data.table
#' @export
classify_labs_dt <- function(labs_path, ref_path) {

  labs <- fread(labs_path)
  ref  <- fread(ref_path)

  ref_unique <- unique(ref[, .(lab, lower, upper)])

  merged_labs <- merge(labs, ref_unique, by = "lab")

  merged_labs[, status := ifelse(value >= lower & value <= upper, "normal", "out_of_range")]

  abnormal_by_patient <- merged_labs[, .(
    total_tests = .N,
    out_of_range = sum(status == "out_of_range")
  ), by = patient_id]

  abnormal_by_lab <- merged_labs[, .(
    total_tests = .N,
    out_of_range = sum(status == "out_of_range")
  ), by = lab]

  return(list(merged_labs = merged_labs, abnormal_by_patient = abnormal_by_patient, abnormal_by_lab = abnormal_by_lab))
}

# ==========================================================
# 🧩 TASK 6 — Match vitals to labs
# ==========================================================

#' Nearest-time matching of vitals and labs
#'
#' Matches nearest HR/SBP readings, vital signs, to each lab time and computes correlations.
#'
#' Performs "rolling join" to match the closest available heart rate (HR)
#' and systolic blood pressure (SBP) to each lab test times for the same patient.
#' Computes per-patient correlations between CRP (C-reactive protein) and vitals.
#'
#' @param labs_path Path to CSV file with columns: patient_id, time_iso, lab, value.
#' @param vitals_path Path to CSV file with columns: time_iso
#' @return List with labs_with_vitals and correlations (CRP vs HR/SBP)
#' @import data.table
#' @export
match_vitals_dt <- function(labs_path, vitals_path) {

  labs <- fread(labs_path)
  vitals <- fread(vitals_path)

  labs[, time_iso := as.POSIXct(time_iso)]
  vitals[, time_iso := as.POSIXct(time_iso)]

  setorder(labs, patient_id, time_iso)
  setorder(vitals, patient_id, time_iso)


  labs[, lab_time := time_iso]
  setkey(labs, patient_id, time_iso)
  setkey(vitals, patient_id, time_iso)


  vitals_hr <- vitals[vital == "HR", .(patient_id, time_iso, value)]
  setnames(vitals_hr, "value", "nearest_HR")

  vitals_hr[, hr_time := time_iso]
  setkey(vitals_hr, patient_id, time_iso)

  labs_with_hr <- vitals_hr[labs, roll = "nearest"]

  labs_with_hr[, hr_lag_minutes := as.numeric(difftime(lab_time, hr_time, units = "mins"))]

  vitals_sbp <- vitals[vital == "SBP", .(patient_id, time_iso, value)]
  setnames(vitals_sbp, "value", "nearest_SBP")

  vitals_sbp[, sbp_time := time_iso]
  setkey(vitals_sbp, patient_id, time_iso)

  labs_with_vitals <- vitals_sbp[labs_with_hr, roll = "nearest"]
  labs_with_vitals[, sbp_lag_minutes := as.numeric(difftime(lab_time, sbp_time, units = "mins"))]


  crp_data <- labs_with_vitals[lab == "CRP"]


  cor_crp_hr <- crp_data[!is.na(nearest_HR), .(
    correlation_CRP_HR = cor(value, nearest_HR, use = "complete.obs")
  ), by = patient_id]

  cor_crp_sbp <- crp_data[!is.na(nearest_SBP), .(
    correlation_CRP_SBP = cor(value, nearest_SBP, use = "complete.obs")
  ), by = patient_id]

  return(list(labs_with_hr = labs_with_hr, labs_with_vitals=labs_with_vitals, cor_crp_hr = cor_crp_hr, cor_crp_sbp = cor_crp_sbp))
}

# ==========================================================
# 🧩 TASK 7 — Top peaks by score
# ==========================================================

#' Extract peaks on chromosome and return top N by score
#'
#' Filters ATAC peaks by chromosome and genomic interval and returns
#' the top 50 by score — useful for quick locus-specific exploration.
#'
#' @param peaks_path Path to CSV file with columns:chr, start, end, peak_id, score
#' @param chr_sel selected chromosome string e.g.:"chr2"
#' @param start_min lower bound of range e.g.:2000000
#' @param start_max upper bound of range e.g.:4000000
#' @return Data.table of top peaks
#' @import data.table
#' @export
top_peaks_dt <- function(peaks_path, chr_sel, start_min , start_max) {

  peaks <- fread("project_oct25/atac_peaks.bed.csv")

  subset_peaks <- peaks[chr == chr_sel & start >= start_min & start <= start_max]

  subset_peaks <- setorder(subset_peaks, -score)

  top50_peaks <- head(subset_peaks, 50)

  return(top50_peaks)
}

# ==========================================================
# 🧩 TASK 8 — Gene stats and filtering
# ==========================================================

#' Compute per-condition gene stas and filter based on fold change.
#'
#' Calculates mean, median, and quartile statistics for each gene under
#' treated and control conditions, then filters genes where the treated
#' mean is at least twice the control mean (treated_mean ≥ 2 × control_mean).
#'
#' @param counts_path Path to CSV file with columns: gene, sample_id, count.
#' @param meta_path Path to CSV file with columns: sample_id, condition, batch, patient ID, timepoint.
#' @return List with stats_by_gene_condition and kept_genes
#' @import data.table
#' @export

gene_stats_filter_dt <- function(counts_path, meta_path) {
  counts <- fread(counts_path)
  meta   <- fread(meta_path)

  merged <- counts[meta, on = "sample_id"]

  stats_by_gene_condition <- merged[, .(
    mean_count   = mean(count),
    median_count = median(count),
    Q1           = quantile(count, 0.25, type = 2),
    Q3           = quantile(count, 0.75, type = 2)
  ), by = .(gene, condition)]

  treated_means <- stats_by_gene_condition[condition == "treated", .(gene, treated_mean = mean_count)]
  control_means <- stats_by_gene_condition[condition == "control", .(gene, control_mean = mean_count)]

  means_wide <- merge(treated_means, control_means, by = "gene", all = FALSE)

  kept_genes <- means_wide[treated_mean >= 2 * control_mean]

  return(list(stats_by_gene_condition = stats_by_gene_condition, kept_genes = kept_genes))
}

# ==========================================================
# 🧩 TASK 9 — Wide → Long → Wide
# ==========================================================

#' Convert wide counts to long and back, computing mean per condition
#'
#'Demonstrates data reshaping with \code{melt()} and \code{dcast()}.
#' Converts wide-format counts (one column per sample) to long format,
#' merges with metadata, computes mean per gene and condition,
#' and returns a condition-wide table again.
#'
#' @param counts_wide_path Path to CSV file with columns genes and samples
#' @return Data.table wide by condition with mean counts
#' @import data.table
#' @export

  wide_long_wide_dt <- function(counts_wide_path, meta_path) {

  counts_wide <- fread(counts_wide_path)

  counts_long <- melt(counts_wide, id.vars = "gene",
                      variable.name = "sample_id", value.name = "count")

  meta <- fread(meta_path)

  merged <- merge(counts_long, meta, by = "sample_id")

  totals_per_sample <- merged[, .(total_count = sum(count)), by = sample_id]

  merged <- merge(merged, totals_per_sample, by = "sample_id")

  gene_condition_means <- merged[, .(mean_count = mean(count)),
                                 by = .(gene, condition)]

  counts_condition_wide <- dcast(gene_condition_means,
                                 gene ~ condition,
                                 value.var = "mean_count")
  return(counts_condition_wide)
}

# ==========================================================
# 🧩 TASK 10 — ATAC peaks to genes
# ==========================================================

#' Map ATAC-seq peaks to genes using genomic overlaps
#'
#' Uses \code{foverlaps()} to detect overlaps between peak intervals
#' and gene coordinates, computes overlap lengths, and summarizes
#' peaks per gene and total overlap length. Returns top 20 genes
#' with the most total overlapping base pairs.
#'
#' @param peaks_path Path to CSV file with columns: chr, start, end, peak_id, score
#' @param genes_path Path to CSV file with columns: chr, start, end, gene
#' @return List with overlap tables, peaks_per_gene and top 20 genes by overlap
#' @import data.table
#' @export

atac_to_gene_dt <- function(peaks_path, genes_path) {

  peaks <- fread(peaks_path)
  genes <- fread(genes_path)

  setkey(peaks, chr, start, end)
  setkey(genes, chr, start, end)

  overlaps <- foverlaps(peaks, genes,
                        by.x = c("chr", "start", "end"),
                        by.y = c("chr", "start", "end"),
                        type = "any", nomatch = 0L)

  overlaps[, overlap_bp := pmin(end, i.end) - pmax(start, i.start)]

  overlaps <- overlaps[overlap_bp > 0]

  peaks_per_gene <- overlaps[, .N, by = gene]
  setnames(peaks_per_gene, "N", "num_peaks")

  overlap_sum_per_gene <- overlaps[, .(total_overlap_bp = sum(overlap_bp)), by = gene]

  top20_genes <- overlap_sum_per_gene[order(-total_overlap_bp)][1:20]
  return(list(overlaps=overlaps, peaks_per_gene=peaks_per_gene, top20_genes=top20_genes))
}

# ==========================================================
# 🧩 TASK 11 — Variants to genes
# ==========================================================

#' Map genetic variants to genes and count HIGH-impact variants
#'
#' Performs overlap join between variant coordinates and gene annotations,
#' counts variants per gene, and identifies genes with HIGH-impact variants.
#'
#' @param variants_path  Path to CSV file with columns: sample_id, chr, pos, ref, alt, impact
#' @param genes_path  Path to CSV file with columns: chr, start, end, gene
#' @return data.table of high-impact genes
#' @import data.table
#' @export
variants_to_genes_dt <- function(variants_path, genes_path) {

  variants <- fread(variants_path)
  genes    <- fread(genes_path)

  variants[, start := pos]  #creo per ogni variante un intervallo 1bp
  variants[, end   := pos]  #foverlaps () lavora con intervalli start-end quindi gli servono

  setkey(variants, chr, start, end)
  setkey(genes,    chr, start, end)

  overlaps <- foverlaps(variants, genes, type = "any", nomatch = 0L)

  overlaps[, impact_upper := toupper(impact)]

  high_overlaps <- overlaps[impact_upper == "HIGH"]


  high_counts_by_gene_sample <- high_overlaps[, .(
    high_variant_count = .N
  ), by = .(gene, sample_id)]


  high_counts_by_gene <- high_overlaps[, .(
    total_high_variants = .N
  ), by = gene][order(-total_high_variants)]


  genes_with_high <- unique(high_counts_by_gene$gene)
  return(genes_with_high = genes_with_high)
}

# ==========================================================
# 🧩 TASK 12 — Combine cohorts
# ==========================================================

#' Combine cohorts safely and compute variability summaries
#'
#' This function merges two cohort-level sample metadata tables
#' and joins them to a long-format RNA-seq counts table to perform a per-cohort, per-condition
#' expression summary of the most 100 variable genes.
#'
#' @param cohortA_path Path to CSV file with columns: sample_id, condition, batch, patient_id, timepoint, cohort
#' @param cohortB_path Path to CSV file with columns: sample_id, condition, batch, patient_id, timepoint, cohort
#' @param counts_path Path to CSV file with columns: gene, sample_id, count.
#' @return List with combined cohort info and top genes
#' @import data.table
#' @export

combine_cohorts_dt <- function(cohortA_path, cohortB_path, counts_path) {

  cohortA <- fread(cohortA_path)
  cohortB <- fread(cohortB_path)
  counts  <- fread(counts_path)

  cohortA[, cohort := "A"]
  cohortB[, cohort := "B"]

  combined_cohorts <- rbindlist(list(cohortA, cohortB), use.names = TRUE, fill = TRUE)

  setorder(combined_cohorts, cohort, condition, sample_id)

  merged_per_sampleid <- merge(counts, combined_cohorts, by = "sample_id", all.x = TRUE)

  gene_variance <- merged_per_sampleid[, .(variance = var(count, na.rm = TRUE)), by = gene]

  top100_genes <- gene_variance[order(-variance)][1:100, gene]

  top100_data <- merged_per_sampleid[gene %in% top100_genes]

  mean_counts <- top100_data[, .(
    mean_count = mean(count, na.rm = TRUE)
  ), by = .(gene, cohort, condition)]

  return(top100_data)
}


# ==========================================================
# 🧩 FINAL REVISION TASKS (1.1 → 5.1)
# ==========================================================

#' Combine integration and clustering data
#'
#' provide a new file where cell type, cells and integration clusters are combined,
#' provide a file where the number of each cell type is indicated for each cluster,
#' provide summary table showing cell types associated to integration clusters and also their
#' association to normal and tumor tissue.
#' Provide a plot describing the distribution of the cell type in normal/tumor tissue given the
#' integration clusters.
#' Provide a normalized % for the cell types in each of the integration clusters for normal and
#' tumor specimen.
#' @param integration_path Path to CSV file with columns: cell, integration_cluster
#' @param clustering_path Path to CSV file with columns: cell, cell_type, sample_type
#' @return files created in Output/
#' @import data.table
#' @export
 final_revision_dt <- function(integration_path, clustering_path) {
  integration_dt <- fread(integration_path)
  clustering_dt  <- fread(clustering_path)

  normalize_cell_id <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    x <- gsub("_X_", "_", x)
    x <- gsub("_Y_", "_", x)
    x <- gsub("_\\.", ".", x)
    return(x)
  }

  integration_dt[, cell_clean := normalize_cell_id(cell)]
  clustering_dt[, cell_clean := normalize_cell_id(cell)]

  # ==========================================
  # TASK 1.1
  # ==========================================

  combined <- merge(integration_dt, clustering_dt, by = "cell_clean", all = FALSE)

  fwrite(combined,
         file = "Outputs/combined_celltype_integrationcluster.csv")

  # ==========================================
  # TASK 2.1
  # ==========================================

  counts_cluster_celltype <- combined[, .N, by = .(integration_cluster, cell_type)]
  setnames(counts_cluster_celltype, "N", "cell_count")

  fwrite(counts_cluster_celltype,
         file = "Outputs/celltype_counts_per_cluster.csv")

  # ==========================================
  # TASK 3.1
  # ==========================================

  # Conta celle per cluster × cell_type × sample_type
  tab_ct_st <- combined[, .N, by = .(integration_cluster, cell_type, sample_type)]
  setnames(tab_ct_st, "N", "cell_count")

  # Totale celle per cluster
  totals_cluster <- combined[, .N, by = integration_cluster]
  setnames(totals_cluster, "N", "cluster_total_cells")

  tab_ct_st <- merge(tab_ct_st, totals_cluster, by = "integration_cluster")

  # Percentuale all'interno del cluster
  tab_ct_st[, pct_within_cluster := (cell_count / cluster_total_cells) * 100]

  # Totale celle per cluster × cell_type
  totals_cluster_celltype <- combined[, .N, by = .(integration_cluster, cell_type)]
  setnames(totals_cluster_celltype, "N", "cluster_celltype_total")

  tab_ct_st <- merge(tab_ct_st, totals_cluster_celltype,
                     by = c("integration_cluster", "cell_type"))

  # Percentuale di sample_type all'interno di quel tipo cellulare nel cluster
  tab_ct_st[, pct_within_cluster_celltype := (cell_count / cluster_celltype_total) * 100]

  fwrite(tab_ct_st,
         file = "Outputs/summary_cluster_celltype_tissue.csv")


  # ==========================================
  # TASK 4.1
  # ==========================================

  plot_data <- as.data.frame(table(combined$integration_cluster,
                                   combined$cell_type,
                                   combined$sample_type))
  names(plot_data) <- c("cluster", "cell_type", "tissue", "count")

  totals <- aggregate(count ~ cluster + tissue, data = plot_data, sum)
  names(totals)[3] <- "total"

  plot_data <- merge(plot_data, totals, by = c("cluster", "tissue"))
  plot_data$percent <- (plot_data$count / plot_data$total) * 100

  plot4_1 <- ggplot(plot_data, aes(x = cluster, y = percent, fill = cell_type)) +
    geom_bar(stat = "identity") +
    facet_wrap(~ tissue) +
    labs(title = "Distribution of cell types in clusters",
         x = "Integration cluster",
         y = "Cells percentage (%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave( "Outputs/plot4_1.jpeg", plot = plot4_1, width = 3000, height = 1500, units = "px")


  # ==========================================
  # TASK 5.1
  # ==========================================
  counts_norm <- combined[, .N, by = .(cluster = integration_cluster,
                                       cell_type,
                                       tissue = sample_type)]

  totals_norm <- counts_norm[, .(total = sum(N)), by = .(cluster, tissue)]

  counts_norm <- merge(counts_norm, totals_norm, by = c("cluster", "tissue"))
  counts_norm[, percent := round((N / total) * 100, 2)]

  fwrite(counts_norm,
         file = "Outputs/normalized_celltype_percentages_by_tissue.csv")

}


