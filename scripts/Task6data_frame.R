# ==========================================================
# TASK 6 – data.frame version
# ==========================================================

# 1. Import
labs   <- read.csv("project_oct25/clinical_labs.csv")
vitals <- read.csv("project_oct25/vitals_time_series.csv")

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


