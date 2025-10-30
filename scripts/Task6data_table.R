# =====================================================
# Goal: Nearest-time matching of vitals to lab draws. Data: clinical_labs.csv, vitals_time_series.csv Tasks:
#• For each lab measurement, attach the nearest HR and SBP reading using a rolling join on (patient_id, time) 
#and report the time lag in minutes.
#• Summarize correlation between CRP and nearest HR/SBP by patient.
#si lavora su dati di pazienti, si introducono misurazioni reali nel tempo,
#si usano strumenti “temporali” come il rolling join per trovare misure vicine nel tempo.

#TASK 6: Nearest-time matching of vitals to lab draws
# =====================================================
library(data.table)

# 1. Carico i dati
labs <- fread("project_oct25/clinical_labs.csv")
vitals <- fread("project_oct25/vitals_time_series.csv")

# 2. Preparo le tabelle e converte le date in un formato tempo che R può confrontare
labs[, time_iso := as.POSIXct(time_iso)]
vitals[, time_iso := as.POSIXct(time_iso)]
#ordina i dati per paziente e per tempo, serve per i join temporali
setorder(labs, patient_id, time_iso)
setorder(vitals, patient_id, time_iso)

# Salvo il tempo del laboratorio, creando una nuova colonna lab_time con l'orario del prelievo
labs[, lab_time := time_iso]
#imposto le chiavi, importante per il join
setkey(labs, patient_id, time_iso)
setkey(vitals, patient_id, time_iso)

# 3. Trovo l'HR più vicino
vitals_hr <- vitals[vital == "HR", .(patient_id, time_iso, value)] #prendo solo le righe dei HR
setnames(vitals_hr, "value", "nearest_HR") #rinomino la colonna value in nearest HR
# SALVO IL TEMPO DELL'HR PRIMA DEL JOIN!
vitals_hr[, hr_time := time_iso] #salvo l'ora della misura HR
setkey(vitals_hr, patient_id, time_iso)

labs_with_hr <- vitals_hr[labs, roll = "nearest"] #rolling join
#per ogni esame di lab, trova il battito HR più vicino nel tempo x lo stesso paziente
labs_with_hr[, hr_lag_minutes := as.numeric(difftime(lab_time, hr_time, units = "mins"))]
#calcola la differenza temporale in minuti tra prelievo e battito

# 4. Trovo l'SBP più vicino, esattamente lo stesso per pressione
vitals_sbp <- vitals[vital == "SBP", .(patient_id, time_iso, value)]
setnames(vitals_sbp, "value", "nearest_SBP")
# SALVO IL TEMPO DELL'SBP PRIMA DEL JOIN!
vitals_sbp[, sbp_time := time_iso]
setkey(vitals_sbp, patient_id, time_iso)

labs_with_vitals <- vitals_sbp[labs_with_hr, roll = "nearest"]
labs_with_vitals[, sbp_lag_minutes := as.numeric(difftime(lab_time, sbp_time, units = "mins"))]

# 6. Analisi CRP: filtra solo le righe dove il lab è CRP
crp_data <- labs_with_vitals[lab == "CRP"]

#per ogni paziente (by sample id) calcola la correlazione tra CRP e battito e CRP e pressione
#complete.obs: ignora le righe con dati mancanti
cor_crp_hr <- crp_data[!is.na(nearest_HR), .(
  correlation_CRP_HR = cor(value, nearest_HR, use = "complete.obs")
), by = patient_id]

cor_crp_sbp <- crp_data[!is.na(nearest_SBP), .(
  correlation_CRP_SBP = cor(value, nearest_SBP, use = "complete.obs")
), by = patient_id]

