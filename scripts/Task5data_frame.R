# 1. Import
labs <- read.csv("project_oct25/clinical_labs.csv")
ref  <- read.csv("project_oct25/lab_reference_ranges.csv")

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
