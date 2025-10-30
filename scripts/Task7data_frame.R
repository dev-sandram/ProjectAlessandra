# ==========================================================
# TASK 7 – data.frame version
# ==========================================================

peaks <- read.csv("project_oct25/atac_peaks.bed.csv", stringsAsFactors = FALSE)
subset_peaks <- subset(peaks, chr == "chr2" & start >= 2000000 & start <= 4000000)
subset_peaks <- subset_peaks[order(-subset_peaks$score), ]
top50_peaks <- head(subset_peaks, 50)

cat("\n--- Top 50 picchi su chr2 (2–4 Mb) ---\n")
print(top50_peaks)
