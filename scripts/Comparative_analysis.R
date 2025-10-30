# COMPARATIVE ANALYSIS: data.frame vs data.table

library(data.table)
library(microbenchmark)
library(ggplot2)

# carico le funzioni
source("R/df_functions.R")
source("R/dt_functions.R")

# Task 1
bench1 <- microbenchmark(
  df_version = bulk_counts_summary_df("project_oct25/bulk_counts_long.csv", "project_oct25/sample_metadata.csv"),
  dt_version = bulk_counts_summary_dt("project_oct25/bulk_counts_long.csv", "project_oct25/sample_metadata.csv"),
  times = 10
)
autoplot(bench1)

# Task 2
bench2 <- microbenchmark(
  df_version = bulk_counts_qc_df("project_oct25/bulk_counts_long.csv"),
  dt_version = bulk_counts_qc_dt("project_oct25/bulk_counts_long.csv"),
  times = 10
)
autoplot(bench2)

#Task 3 non ha il benchmark perchè ha senso solo in data.table

# Task 4
bench4 <- microbenchmark(
  df_version = annotate_counts_df("project_oct25/bulk_counts_long.csv", "project_oct25/sample_metadata.csv"),
  dt_version = annotate_counts_dt("project_oct25/bulk_counts_long.csv", "project_oct25/sample_metadata.csv"),
  times = 10
)
autoplot(bench4)

# Task 5
bench5 <- microbenchmark(
  df_version = classify_labs_df("project_oct25/clinical_labs.csv", "project_oct25/lab_reference_ranges.csv"),
  dt_version = classify_labs_dt("project_oct25/clinical_labs.csv", "project_oct25/lab_reference_ranges.csv"),
  times = 10
)
autoplot(bench5)

 # Task 6
bench6 <- microbenchmark(
  df_version = match_vitals_df("project_oct25/clinical_labs.csv", "project_oct25/vitals_time_series.csv"),
  dt_version = match_vitals_dt("project_oct25/clinical_labs.csv", "project_oct25/vitals_time_series.csv"),
  times = 10
)
autoplot(bench6)

# Task 7
bench7 <- microbenchmark(
  df_version = top_peaks_df("project_oct25/atac_peaks.bed.csv", "chr2", 2000000, 4000000),
  dt_version = top_peaks_dt("project_oct25/atac_peaks.bed.csv", "chr2", 2000000, 4000000),
  times = 10
)
autoplot(bench7)

# Task 8
bench8 <- microbenchmark(
  df_version = gene_stats_filter_df("project_oct25/bulk_counts_long.csv", "project_oct25/sample_metadata.csv"),
  dt_version = gene_stats_filter_dt("project_oct25/bulk_counts_long.csv", "project_oct25/sample_metadata.csv"),
  times = 10
)
autoplot(bench8)

# Task 9
bench9 <- microbenchmark(
  df_version = wide_long_wide_df("project_oct25/bulk_counts_wide.csv", "project_oct25/sample_metadata.csv"),
  dt_version = wide_long_wide_dt("project_oct25/bulk_counts_wide.csv", "project_oct25/sample_metadata.csv"),
  times = 10
)
autoplot(bench9)

# Task 10
bench10 <- microbenchmark(
  df_version = atac_to_gene_df("project_oct25/atac_peaks.bed.csv", "project_oct25/gene_annotation.bed.csv"),
  dt_version = atac_to_gene_dt("project_oct25/atac_peaks.bed.csv", "project_oct25/gene_annotation.bed.csv"),
  times = 10
)
autoplot(bench10)

# Task 11
bench11 <- microbenchmark(
  df_version = variants_to_genes_df("project_oct25/variants.csv", "project_oct25/gene_annotation.bed.csv"),
  dt_version = variants_to_genes_dt("project_oct25/variants.csv", "project_oct25/gene_annotation.bed.csv"),
  times = 10
)
autoplot(bench11)

# Task 12
bench12 <- microbenchmark(
  df_version = combine_cohorts_df("project_oct25/cohortA_samples.csv", "project_oct25/cohortB_samples.csv","project_oct25/bulk_counts_long.csv"),
  dt_version = combine_cohorts_dt("project_oct25/cohortA_samples.csv", "project_oct25/cohortB_samples.csv","project_oct25/bulk_counts_long.csv"),
  times = 10
)
autoplot(bench12)

#FINAL REVISION
benchfr <- microbenchmark(
  df_version = final_revision_df("project_oct25/annotated_GSM3516673_normal_annotated_GSM3516672_tumor_SeuratIntegration.csv", "project_oct25/nt_combined_clustering.output.csv"),
  dt_version = final_revision_dt ("project_oct25/annotated_GSM3516673_normal_annotated_GSM3516672_tumor_SeuratIntegration.csv", "project_oct25/nt_combined_clustering.output.csv"),
  times = 10
)
autoplot(benchfr)

results <- rbindlist(list(
  task1 = as.data.table(bench1),
  task2 = as.data.table(bench2),
  task4 = as.data.table(bench4),
  task5 = as.data.table(bench5),
  task6 = as.data.table(bench6),
  task7 = as.data.table(bench7),
  task8 = as.data.table(bench8),
  task9 = as.data.table(bench9),
  task10 = as.data.table(bench10),
  task11 = as.data.table(bench11),
  task12 = as.data.table(bench12),
  final_revision = as.data.table(benchfr)
), idcol = "task")[, .(df_mean = mean(time[expr=="df_version"])/1e6,
                       dt_mean = mean(time[expr=="dt_version"])/1e6),
                   by = task]

