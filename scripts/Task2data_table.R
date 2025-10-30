#Goal: Add QC-style derived columns without copying. Data: bulk_counts_long.csv 
#Tasks: 
#• Add la log2 counts column and a binary flag high if count > 100. 
#• Overwrite high to use a gene-wise threshold (e.g., count > median(count) by=gene).

# 1 Carico la libreria
library(data.table)

# 2 Leggo il file come data.table
counts <- fread("project_oct25/bulk_counts_long.csv")

# PARTE 1: Aggiungo la colonna log2 dei conteggi
# -----------------------------------------------------
# := modifica la tabella "in place" (senza crearne una copia)
counts[, log2_count := log2(count)]  #modifica o calcola qualcosa in questa tabella
#:= significa che la modifica è in place
#calcola il log in base due dei valori nella colonna count per stabilizzare la varianza
#rende i dati più confrontabili tra geni
# posso aggiungere +1 dopo count per evitare log2(0), che non è definito)

#Aggiungo la colonna binaria 'high' (count > 100)
counts[, high := count > 100]
#se il count è >100 restituisce TRUE altrimenti FALSE
#quindi crea una nuova colonna high con valori booleani
#per classificare i geni in highly expressed o no usando 100 come soglia

# -----------------------------------------------------
# PARTE 3: Sovrascrivo 'high' in modo gene-wise, con questo criterio := (in place)
#per ogni gene verifica se il valore di count è maggiore della mediana di count per quel gene
# (count > median(count) per ciascun gene)
# -----------------------------------------------------
counts[, high := count > median(count), by = gene]  #x ogni gene separatamente


