#Goal: Classify values against reference intervals.
#Data: clinical_labs.csv, lab_reference_ranges.csv, plus sample_metadata.csv (for sex proxy if you want to extend)
#Tasks:
#• Treat reference ranges as intervals and label each lab as "normal" vs "out_of_range" using a non-equi join 
#(value >= lower & value <= upper).
#• Count abnormal rates by patient and by lab.

library(data.table)

# 1 Carico i file
labs <- fread("project_oct25/clinical_labs.csv")
ref  <- fread("project_oct25/lab_reference_ranges.csv")

# Unisco le tabelle per aggiungere i range di riferimento
# Dato che i range sono uguali per M e F, li prendo una sola volta
ref_unique <- unique(ref[, .(lab, lower, upper)]) 
#unique returns a data table with duplicated rows removed.
#unique è una funzione che prende la tabella lab_ref e controlla quali sono le differenze nelle colonne lab, 
#lower e upper e le salva in ref_unique

# Uso una join per aggiungere a ogni valore clinico i limiti lower e upper del test corrispondente
merged_labs <- merge(labs, ref_unique, by = "lab")
# -----------------------------------------------------
# PARTE 2: Classifico i valori come "normal" o "out_of_range"
# -----------------------------------------------------
merged_labs[, status := ifelse(value >= lower & value <= upper, "normal", "out_of_range")] 
#all'interno di merged labs, a dx, metti status con questa condizione, 
#rispettivamente, se out o normale e le salva in una colonna che è status

# -----------------------------------------------------
# PARTE 3: Calcolo la percentuale di risultati fuori range per paziente
# -----------------------------------------------------
abnormal_by_patient <- merged_labs[, .(
  total_tests = .N,                                 # numero totale di test per paziente
  out_of_range = sum(status == "out_of_range")      # quanti sono fuori range
), by = patient_id]

# -----------------------------------------------------
# PARTE 4: Calcolo i risultati fuori range per tipo di test
# -----------------------------------------------------
abnormal_by_lab <- merged_labs[, .(
  total_tests = .N,
  out_of_range = sum(status == "out_of_range")
), by = lab]
