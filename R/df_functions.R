# ==========================================================
# df_functions.R
# Generalized and reusable base R (data.frame) functions
# for bulk RNA-seq processing and summarization
# ==========================================================

# ----------------------------------------------------------
# TASK 1 – Merge, filter, summarize, and compute per-condition means
# ----------------------------------------------------------
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
#' @return none
#' @export
bulk_counts_summary_df <- function(counts_path, meta_path) {
  # 1. Import
  counts <- read.csv(counts_path)
  meta   <- read.csv(meta_path)

  # 2. Join counts + metadata by sample_id
  merged_data <- merge(counts, meta, by = "sample_id")

  # 3. Filter for treated samples and GENE_00* genes
  treated_data <- subset(merged_data,
                         condition == "treated" & grepl("^GENE_00", gene))

  # 4. Compute mean and median per gene
  gene_summary <- aggregate(
    count ~ gene,
    data = treated_data,
    FUN  = function(x) c(mean = mean(x), median = median(x))
  )

  # 5. Simplify format
  gene_summary_df <- data.frame(
    gene = gene_summary$gene,
    mean_count   = gene_summary$count[, "mean"],
    median_count = gene_summary$count[, "median"]
  )

  # 6 Calcolo della media dei conteggi per ciascun gene e condizione
  # aggregate() calcola una funzione (mean) per gruppi
  gene_condition_means_f <- aggregate(
    count ~ gene + condition,   # formula: raggruppa per gene e condition
    data = merged_data,         # tabella di partenza
    FUN = mean,                 # funzione da applicare
    na.rm = TRUE                # ignora eventuali valori mancanti
  )

  # 7 Ordina i risultati per gene e condition (per leggibilità)
  gene_condition_means_f <- gene_condition_means_f[order(gene_condition_means_f$gene,
                                                         gene_condition_means_f$condition), ]

}

# ----------------------------------------------------------
# TASK 2 – Add log2(count + 1) and binary flags
# ----------------------------------------------------------
#' Add QC-style derived columns
#'
#' Computes log2-transformed counts and adds a binary \code{high} flag
#' indicating whether a count is above the gene-wise median.
#' This mimics typical preprocessing steps in QC pipelines to visualize
#' count distributions and detect outliers.
#'
#' @param counts_path Path to CSV file with columns: gene, sample_id, count.
#' @return none
#' @export
bulk_counts_qc_df <- function(counts_path) {
  # 1. Import
  counts <- read.csv(counts_path)

  # 2. Add log2(count + 1) column
  counts$log2_count <- log2(counts$count + 1)

  # 3. Add binary flag 'high' (count > 100)
  counts$high <- counts$count > 100

  # 4. Overwrite 'high' using gene-wise median threshold
  #    Qui usiamo tapply() per calcolare la mediana per gene,
  #    poi combiniamo i risultati in un vettore logico della stessa lunghezza
  medians_by_gene <- tapply(counts$count, counts$gene, median)
  counts$high <- counts$count > medians_by_gene[counts$gene]
}

#-----------------------------------------------------------
#Task 3 non ha la funzione perchè ha senso solo in data.table
#-----------------------------------------------------------

# ----------------------------------------------------------
# TASK 4 – Filter genes by thresholded normalized expression
# ----------------------------------------------------------
#' Annotate counts with metadata and compute summaries
#'
#' Adds patient metadata to bulk counts, computes per-patient total counts,
#' per-gene means by condition, and identifies top 10 genes per condition
#' ranked by average expression.
#'
#' @param counts_path Path to CSV file with columns: gene, sample_id, count.
#' @param meta_path Path to CSV file with columns: sample_id, condition, batch, patient ID, timepoint.
#' @return none
#' @export
annotate_counts_df <- function(counts_path, meta_path) {
  counts <- read.csv(counts_path)
  meta   <- read.csv(meta_path, stringsAsFactors = FALSE)

  # 2. Join
  merged_data <- merge(counts, meta, by = "sample_id")

  # 3. Total counts per patient
  patient_totals <- aggregate(count ~ patient_id, data = merged_data, FUN = sum)
  names(patient_totals)[2] <- "total_count"

  # 4. Mean count per gene and condition
  gene_means <- aggregate(count ~ gene + condition, data = merged_data, FUN = mean)
  names(gene_means)[3] <- "mean_count"

  # 5. Top 10 genes by condition
  conditions <- unique(gene_means$condition)
  top10_by_condition <- data.frame()

  for (cond in conditions) {
    subset_cond <- subset(gene_means, condition == cond)
    subset_cond <- subset_cond[order(-subset_cond$mean_count), ]
    top10 <- head(subset_cond, 10)
    top10_by_condition <- rbind(top10_by_condition, top10)}
}

# ----------------------------------------------------------
# TASK 5 – Summary table by condition and gene class
# ----------------------------------------------------------
#' Classify lab values against reference intervals
#'
#' Merges patient lab results with reference intervals and assigns each test
#' as "normal" or "out_of_range", and summarizes abnormalities per patient and per lab.
#'
#' @param labs_path Path to CSV file with columns: patient_id, time_iso, lab, value.
#' @param ref_path Path to CSV file with columns: lab, sex,lower, upper.
#' @return none
#' @export
classify_labs_df <- function(labs_path, ref_path) {
  # 1. Import
  labs <- read.csv(labs_path)
  ref  <- read.csv(ref_path)

   # 2. Keep unique reference intervals
  ref_unique <- unique(ref[, c("lab", "lower", "upper")])

  # 3. Join
  merged_labs <- merge(labs, ref_unique, by = "lab")

  # 4. Classify as normal/out_of_range
  merged_labs$status <- ifelse(merged_labs$value >= merged_labs$lower &
                                 merged_labs$value <= merged_labs$upper,
                               "normal", "out_of_range")

  # 5. Summaries
  abnormal_by_patient <- aggregate(status ~ patient_id, data = merged_labs,
                                   FUN = function(x) {
                                     total <- length(x)
                                     out   <- sum(x == "out_of_range")
                                     c(total_tests = total, out_of_range = out)
                                   })

  # split matrix-like column into two numeric columns
  abnormal_by_patient$total_tests  <- abnormal_by_patient$status[, "total_tests"]
  abnormal_by_patient$out_of_range <- abnormal_by_patient$status[, "out_of_range"]
  abnormal_by_patient$status <- NULL

  abnormal_by_lab <- aggregate(status ~ lab, data = merged_labs,
                               FUN = function(x) {
                                 total <- length(x)
                                 out   <- sum(x == "out_of_range")
                                 c(total_tests = total, out_of_range = out)
                               })

  abnormal_by_lab$total_tests  <- abnormal_by_lab$status[, "total_tests"]
  abnormal_by_lab$out_of_range <- abnormal_by_lab$status[, "out_of_range"]
  abnormal_by_lab$status <- NULL

}


# ----------------------------------------------------------
# TASK 6 – Compute fold change between treated and control
# ----------------------------------------------------------
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
#' @return none
#' @export
match_vitals_df <- function(labs_path, vitals_path) {
  # 1. Import
  labs   <- read.csv(labs_path)
  vitals <- read.csv(vitals_path)

  # 2. Convert to POSIXct
  labs$time_iso   <- as.POSIXct(labs$time_iso)
  vitals$time_iso <- as.POSIXct(vitals$time_iso)

  # 3. Join nearest HR/SBP manually (loop-based)
  nearest_match <- function(lab_df, vitals_df, vital_type) {
    vitals_sub <- vitals_df[vitals_df$vital == vital_type, ]
    results <- list()
    for (i in seq_len(nrow(lab_df))) {
      pid <- lab_df$patient_id[i]
      time <- lab_df$time_iso[i]
      v_sub <- vitals_sub[vitals_sub$patient_id == pid, ]
      if (nrow(v_sub) > 0) {
        idx <- which.min(abs(difftime(v_sub$time_iso, time, units = "mins")))
        results[[i]] <- v_sub[idx, c("value", "time_iso")]
      } else {
        results[[i]] <- data.frame(value = NA, time_iso = NA)
      }
    }
    out <- do.call(rbind, results)
    names(out) <- c(paste0("nearest_", vital_type), paste0(vital_type, "_time"))
    return(out)
  }

  # 4. Attach nearest HR and SBP
  hr_data  <- nearest_match(labs, vitals, "HR")
  sbp_data <- nearest_match(labs, vitals, "SBP")

  labs_with_vitals <- cbind(labs, hr_data, sbp_data)
  labs_with_vitals$hr_lag_minutes  <- as.numeric(difftime(labs_with_vitals$time_iso, labs_with_vitals$HR_time, units = "mins"))
  labs_with_vitals$sbp_lag_minutes <- as.numeric(difftime(labs_with_vitals$time_iso, labs_with_vitals$SBP_time, units = "mins"))

  # 5. Correlation per patient (CRP only) — robust version in base R
  crp_data <- subset(labs_with_vitals, lab == "CRP")

  # Split by patient_id
  patients_list <- split(crp_data, crp_data$patient_id)

  # Function to compute correlation safely
  safe_cor <- function(df, col_x = "value", col_y = "nearest_HR") {
    if (nrow(df) < 2) return(NA_real_)           # meno di 2 osservazioni -> NA
    x <- df[[col_x]]
    y <- df[[col_y]]
    # keep only pairs non-NA
    ok <- !is.na(x) & !is.na(y)
    if (sum(ok) < 2) return(NA_real_)            # meno di 2 coppie valide -> NA
    return(cor(x[ok], y[ok], use = "complete.obs"))
  }

  # Calcola correlazione CRP vs HR per ogni paziente
  cor_crp_hr <- data.frame(
    patient_id = names(patients_list),
    correlation_CRP_HR = sapply(patients_list, safe_cor, col_x = "value", col_y = "nearest_HR"),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  # Calcola correlazione CRP vs SBP per ogni paziente
  cor_crp_sbp <- data.frame(
    patient_id = names(patients_list),
    correlation_CRP_SBP = sapply(patients_list, safe_cor, col_x = "value", col_y = "nearest_SBP"),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

# ----------------------------------------------------------
# TASK 7 – Compute coefficient of variation (CV)
# ----------------------------------------------------------
#' Extract peaks on chromosome and return top N by score
#'
#' Filters ATAC peaks by chromosome and genomic interval and returns
#' the top 50 by score — useful for quick locus-specific exploration.
#'
#' @param peaks_path Path to CSV file with columns:chr, start, end, peak_id, score
#' @param chr_sel selected chromosome string e.g.:"chr2"
#' @param start_min lower bound of range e.g.:2000000
#' @param start_max upper bound of range e.g.:4000000
#' @return none
#' @export
top_peaks_df <- function(peaks_path, chr_sel, start_min , start_max) {
  peaks <- read.csv(peaks_path, stringsAsFactors = FALSE)
  subset_peaks <- subset(peaks, chr == chr_sel & start >= start_min & start <= start_max)
  subset_peaks <- subset_peaks[order(-subset_peaks$score), ]
  top50_peaks <- head(subset_peaks, 50)
}

# ----------------------------------------------------------
# TASK 8 – Identify most variable genes
# ----------------------------------------------------------
#' Compute per-condition gene stas and filter based on fold change.
#'
#' Calculates mean, median, and quartile statistics for each gene under
#' treated and control conditions, then filters genes where the treated
#' mean is at least twice the control mean (treated_mean ≥ 2 × control_mean).
#'
#' @param counts_path Path to CSV file with columns: gene, sample_id, count.
#' @param meta_path Path to CSV file with columns: sample_id, condition, batch, patient ID, timepoint.
#' @return none
#' @export
gene_stats_filter_df <- function(counts_path, meta_path) {
  counts <- read.csv(counts_path, stringsAsFactors = FALSE)
  meta   <- read.csv(meta_path, stringsAsFactors = FALSE)

  merged <- merge(counts, meta, by = "sample_id")

  stats_by_gene_condition <- aggregate(count ~ gene + condition, data = merged,
                                       FUN = function(x) c(mean = mean(x), median = median(x),
                                                           Q1 = quantile(x, 0.25, type = 2),
                                                           Q3 = quantile(x, 0.75, type = 2)))

  stats_df <- data.frame(
    gene = stats_by_gene_condition$gene,
    condition = stats_by_gene_condition$condition,
    mean_count = stats_by_gene_condition$count[, "mean"],
    median_count = stats_by_gene_condition$count[, "median"],
    Q1 = stats_by_gene_condition$count[, "Q1.25%"],
    Q3 = stats_by_gene_condition$count[, "Q3.75%"]
  )

  treated <- subset(stats_df, condition == "treated")[, c("gene", "mean_count")]
  control <- subset(stats_df, condition == "control")[, c("gene", "mean_count")]
  names(treated)[2] <- "treated_mean"
  names(control)[2] <- "control_mean"

  means_wide <- merge(treated, control, by = "gene")
  kept_genes <- subset(means_wide, treated_mean >= 2 * control_mean)
}

# ----------------------------------------------------------
# TASK 9 – Compute correlation matrix across samples
# ----------------------------------------------------------
#' Convert wide counts to long and back, computing mean per condition
#'
#'Demonstrates data reshaping with \code{melt()} and \code{dcast()}.
#' Converts wide-format counts (one column per sample) to long format,
#' merges with metadata, computes mean per gene and condition,
#' and returns a condition-wide table again.
#'
#' @param counts_wide_path Path to CSV file with columns genes and samples
#' @return none
#' @export
wide_long_wide_df <- function(counts_path, meta_path) {
  counts_wide <- read.csv(counts_path, stringsAsFactors = FALSE)
  counts_long <- reshape(counts_wide, varying = names(counts_wide)[-1],
                         v.names = "count", timevar = "sample_id",
                         times = names(counts_wide)[-1], idvar = "gene",
                         direction = "long")

  meta <- read.csv(meta_path, stringsAsFactors = FALSE)
  merged <- merge(counts_long, meta, by = "sample_id")

  totals_per_sample <- aggregate(count ~ sample_id, data = merged, sum)
  names(totals_per_sample)[2] <- "total_count"

  merged <- merge(merged, totals_per_sample, by = "sample_id")

  gene_condition_means <- aggregate(count ~ gene + condition, data = merged, mean)
  counts_condition_wide <- reshape(gene_condition_means,
                                   timevar = "condition", idvar = "gene",
                                   direction = "wide")
}

# ----------------------------------------------------------
# TASK 10 – PCA on normalized data
# ----------------------------------------------------------
#' Map ATAC-seq peaks to genes using genomic overlaps
#'
#' Uses \code{foverlaps()} to detect overlaps between peak intervals
#' and gene coordinates, computes overlap lengths, and summarizes
#' peaks per gene and total overlap length. Returns top 20 genes
#' with the most total overlapping base pairs.
#'
#' @param peaks_path Path to CSV file with columns: chr, start, end, peak_id, score
#' @param genes_path Path to CSV file with columns: chr, start, end, gene
#' @return none
#' @export
atac_to_gene_df <- function(peaks_path, genes_path) {
  peaks <- read.csv(peaks_path, stringsAsFactors = FALSE)
  genes <- read.csv(genes_path, stringsAsFactors = FALSE)

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
}

# ----------------------------------------------------------
# TASK 11 – Plot gene expression distribution
# ----------------------------------------------------------
#' Map genetic variants to genes and count HIGH-impact variants
#'
#' Performs overlap join between variant coordinates and gene annotations,
#' counts variants per gene, and identifies genes with HIGH-impact variants.
#'
#' @param variants_path  Path to CSV file with columns: sample_id, chr, pos, ref, alt, impact
#' @param genes_path  Path to CSV file with columns: chr, start, end, gene
#' @return none
#' @export
variants_to_genes_df <- function(variants_path, genes_path) {

  variants <- read.csv(variants_path)
  genes    <- read.csv(genes_path)

  # Creo intervalli 1-bp per le varianti -----------------------------------
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

  # Normalizzo la colonna impact a maiuscolo --------------------------------
  if (nrow(overlaps) > 0) {
    overlaps$impact_upper <- toupper(overlaps$impact)
  }

  # Filtro solo le varianti di alto impatto ---------------------------------
  high_overlaps <- subset(overlaps, impact_upper == "HIGH")

  # Conteggio delle varianti HIGH per gene e per sample ---------------------
  if (nrow(high_overlaps) > 0) {
    high_counts_by_gene_sample <- aggregate(
      x = list(high_variant_count = high_overlaps$impact_upper),
      by = list(gene = high_overlaps$gene, sample_id = high_overlaps$sample_id),
      FUN = length
    )
  } else {
    high_counts_by_gene_sample <- data.frame()
  }

  # Conteggio totale delle varianti HIGH per gene --------------------------
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

  # Salvo i risultati su file CSV ------------------------------------------
  write.csv(high_counts_by_gene_sample, "project_oct25/high_variants_by_gene_sample.csv", row.names = FALSE)
  write.csv(high_counts_by_gene, "project_oct25/high_variants_by_gene_total.csv", row.names = FALSE)
  write.csv(data.frame(gene = genes_with_high), "project_oct25/genes_with_high_variants.csv", row.names = FALSE)
}

# ----------------------------------------------------------
# TASK 12 – Compute summary statistics
# ----------------------------------------------------------
#' Combine cohorts safely and compute variability summaries
#'
#' This function merges two cohort-level sample metadata tables
#' and joins them to a long-format RNA-seq counts table to perform a per-cohort, per-condition
#' expression summary of the most 100 variable genes.
#'
#' @param cohortA_path Path to CSV file with columns: sample_id, condition, batch, patient_id, timepoint, cohort
#' @param cohortB_path Path to CSV file with columns: sample_id, condition, batch, patient_id, timepoint, cohort
#' @param counts_path Path to CSV file with columns: gene, sample_id, count.
#' @return none
#' @export
combine_cohorts_df <- function(cohortA_path, cohortB_path, counts_path) {
  cohortA <- read.csv(cohortA_path, stringsAsFactors = FALSE)
  cohortB <- read.csv(cohortB_path, stringsAsFactors = FALSE)
  counts  <- read.csv(counts_path, stringsAsFactors = FALSE)

  cohortA$cohort <- "A"
  cohortB$cohort <- "B"

  combined_cohorts <- rbind(cohortA, cohortB)
  combined_cohorts <- combined_cohorts[order(combined_cohorts$cohort,
                                             combined_cohorts$condition,
                                             combined_cohorts$sample_id), ]

  merged <- merge(counts, combined_cohorts, by = "sample_id", all.x = TRUE)

  gene_variance <- aggregate(count ~ gene, data = merged, FUN = var, na.rm = TRUE)
  names(gene_variance)[2] <- "variance"
  top100_genes <- head(gene_variance[order(-gene_variance$variance), "gene"], 100)

  top100_data <- subset(merged, gene %in% top100_genes)
  mean_counts <- aggregate(count ~ gene + cohort + condition,
                           data = top100_data, FUN = mean, na.rm = TRUE)
  names(mean_counts)[4] <- "mean_count"
}


# ----------------------------------------------------------
# FINAL REVISION TASKS (1.1 -> 5.1)
# ----------------------------------------------------------

#TASK 1
#' Combine integration and clustering data
#'
#' Reads and merges integration and clustering tables by cell ID, ensuring consistent ID formatting.
#' @param integration_path Path to CSV file with columns: cell, integration_cluster
#' @param clustering_path Path to CSV file with columns: cell, cell_type, sample_type
#' @return none
#' @export
final_revision_df <- function(integration_path, clustering_path) {
  # Funzione di normalizzazione ID cellule
  normalize_cell_id <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    x <- gsub("_X_", "_", x)
    x <- gsub("_Y_", "_", x)
    x <- gsub("_\\.", ".", x)
    return(x)
  }

  # ==========================================
  # TASK 1.1
  # ==========================================

  # Leggi i file in data.frame
  integration_df <- read.csv(integration_path)
  clustering_df  <- read.csv(clustering_path)

  # Normalizza gli ID
  integration_df$cell_clean <- normalize_cell_id(integration_df$cell)
  clustering_df$cell_clean  <- normalize_cell_id(clustering_df$cell)

  # Join sulle cellule
  combined_df <- merge(integration_df, clustering_df, by = "cell_clean", all = FALSE)

  # Salva file Task 1
  write.csv(combined_df, "Outputs/combined_cell_cluster_celltype.csv", row.names = FALSE)

  # ==========================================
  # TASK 2.1
  # ==========================================

  counts_cluster_celltype <- aggregate(cell_clean ~ integration_cluster + cell_type,
                                       data = combined_df, FUN = length)
  colnames(counts_cluster_celltype)[3] <- "cell_count"

  write.csv(counts_cluster_celltype, "Outputs/counts_per_cluster_celltype.csv", row.names = FALSE)

  # ==========================================
  # TASK 3.1
  # ==========================================

  # Conta per cluster, cell_type, sample_type
  tab_ct_st <- aggregate(cell_clean ~ integration_cluster + cell_type + sample_type,
                         data = combined_df, FUN = length)
  colnames(tab_ct_st)[4] <- "cell_count"

  # Totale per cluster
  totals_cluster <- aggregate(cell_clean ~ integration_cluster,
                              data = combined_df, FUN = length)
  colnames(totals_cluster)[2] <- "cluster_total_cells"

  # Merge per aggiungere totale cluster
  tab_ct_st <- merge(tab_ct_st, totals_cluster, by = "integration_cluster", all.x = TRUE)

  # Percentuale dentro cluster
  tab_ct_st$pct_within_cluster <- (tab_ct_st$cell_count / tab_ct_st$cluster_total_cells) * 100

  # Totale per (cluster, cell_type)
  totals_cluster_celltype <- aggregate(cell_clean ~ integration_cluster + cell_type,
                                       data = combined_df, FUN = length)
  colnames(totals_cluster_celltype)[3] <- "cluster_celltype_total"

  # Merge + percentuale sample_type dentro cell_type del cluster
  tab_ct_st <- merge(tab_ct_st, totals_cluster_celltype,
                     by = c("integration_cluster", "cell_type"), all.x = TRUE)

  tab_ct_st$pct_within_cluster_celltype <- (tab_ct_st$cell_count / tab_ct_st$cluster_celltype_total) * 100

  # Salva Task 3
  write.csv(tab_ct_st, "Outputs/summary_cluster_celltype_sampletype.csv", row.names = FALSE)

  # ==========================================
  # TASK 4.1
  # ==========================================

  plot_data <- as.data.frame(table(combined_df$integration_cluster,
                                   combined_df$cell_type,
                                   combined_df$sample_type))
  colnames(plot_data) <- c("cluster", "cell_type", "tissue", "count")

  # Totale per cluster e tessuto
  totals <- aggregate(count ~ cluster + tissue, data = plot_data, sum)
  colnames(totals)[3] <- "total"

  # Merge + percentuali
  plot_data <- merge(plot_data, totals, by = c("cluster", "tissue"))
  plot_data$percent <- (plot_data$count / plot_data$total) * 100

  # Plot
  plot4_1_df <- ggplot(plot_data, aes(x = cluster, y = percent, fill = cell_type)) +
    geom_bar(stat = "identity") +
    facet_wrap(~ tissue) +
    labs(title = "Distribution of cell types in clusters (Normal vs Tumor)",
         x = "Integration Cluster",
         y = "Cell Percentage (%)") +
    theme_minimal()

  ggsave( "Outputs/plot4_1_df.jpeg", plot = plot4_1_df, width = 3000, height = 1500, units = "px")


  # ==========================================
  # TASK 5.1
  # ==========================================
  # Conta celle
  counts <- as.data.frame(table(combined_df$integration_cluster,
                                combined_df$cell_type,
                                combined_df$sample_type))
  colnames(counts) = c("cluster", "cell_type", "tissue", "count")

  # Totale per cluster e tessuto
  totals <- aggregate(count ~ cluster + tissue, data = counts, sum)
  colnames(totals)[3] <- "total"

  # Merge e calcolo % normalizzata
  counts_norm <- merge(counts, totals, by = c("cluster", "tissue"))
  counts_norm$percent <- round((counts_norm$count / counts_norm$total) * 100, 2)

  # Salva Task 5
  write.csv(counts_norm, "Outputs/normalized_percent_cluster_celltype_tissue.csv", row.names = FALSE)

}


