#Goal: Slice genomics windows efficiently. Data: atac_peaks.bed.csv
#Tasks:
#• Extract peaks on chr2 with start between 2 and 4 Mb.
#• Among those peaks, return the top 50 by score after setorder() (descending)

# Carico la libreria
library(data.table)

# Leggo il file dei picchi ATAC
peaks <- fread("project_oct25/atac_peaks.bed.csv")

# Filtro i picchi su chr2 e nella finestra 2–4 Mb
# Nota: 1 Mb = 1.000.000 basi
subset_peaks <- peaks[chr == "chr2" & start >= 2000000 & start <= 4000000]

#Ordino i picchi per punteggio decrescente
subset_peaks <- setorder(subset_peaks, -score)

#Prendo i primi 50 picchi
top50_peaks <- head(subset_peaks, 50)

result <- top_peaks_dt("project_oct25/atac_peaks.bed.csv","chr2", 2000000,4000000)
