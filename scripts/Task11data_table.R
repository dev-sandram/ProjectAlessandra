#Goal: Map SNPs to genes.
#Data: variants.csv, gene_annotation.bed.csv 
#Tasks:
#• Convert variant positions to 1-bp intervals (pos, pos) 
#and find overlaps with gene intervals.
#• Summarize counts of HIGH impact variants by gene and by sample.
#• List genes with HIGH-impact variants across all samples.

#Obiettivo: 
#mappare le varianti (SNP) sui geni
#contare le varianti ad alto impatto (HIGH) per gene e per campione
#ottenere la lista dei geni che presentano almeno una variante HIGH in qualsiasi campione

# carico la libreria
library(data.table)

# leggo i file
variants <- fread("project_oct25/variants.csv")
genes    <- fread("project_oct25/gene_annotation.bed.csv")

# Creo intervalli 1-bp per le varianti: start = pos, end = pos
#    (foverlaps lavora con intervalli start-end)
variants[, start := pos]  #creo per ogni variante un intervallo 1bp 
variants[, end   := pos]  #foverlaps () lavora con intervalli start-end quindi gli servono

# 5) imposto le chiavi per foverlaps: (chr, start, end), requisito per foverlaps utile per velocizzare
setkey(variants, chr, start, end)
setkey(genes,    chr, start, end)

# 6) eseguo l'overlap: trovo per ogni variante i geni che si sovrappongono
#    nomatch = 0L rimuove le varianti senza overlap, che non sovrappongono alcun gene
overlaps <- foverlaps(variants, genes, type = "any", nomatch = 0L)
#type any cioè prende qualsiasi sovrapposizione (parziale o completa)

# 7) filtro le varianti di alto impatto (impact == "HIGH")
#    normalizzo il valore di impact a maiuscole per sicurezza
overlaps[, impact_upper := toupper(impact)] #normalizzo il campo impact
#overlaps	È il nome della tabella su cui stai lavorando (una data.table).
#[,]	In data.table, serve per dire “fai qualcosa su questa tabella”.
#impact_upper := ...	Con := stai creando una nuova colonna chiamata impact_upper, oppure sovrascrivendo se già esiste.
#toupper(impact)	È una funzione di R che trasforma tutto il testo in maiuscolo prendendo i valori della colonna impact.

high_overlaps <- overlaps[impact_upper == "HIGH"] #seleziono solo le varianti con impact HIGH

# 8) conto le varianti HIGH per gene e per sample_id
#conta quante varianti ci sono per coppia gene x sample_id
#.N è il conteggio delle righe del gruppo
high_counts_by_gene_sample <- high_overlaps[, .(
  high_variant_count = .N
), by = .(gene, sample_id)]

# 9) conto le varianti HIGH per gene (tutte le sample aggregate)
high_counts_by_gene <- high_overlaps[, .(
  total_high_variants = .N
), by = gene][order(-total_high_variants)]
#conta quanti HIGH per gene complessivamente e ordina decrescente i geni

# 10) lista dei geni che hanno almeno una variante HIGH (unique gene)
#Estrae i geni con almeno una variante HIGH
genes_with_high <- unique(high_counts_by_gene$gene)
print(genes_with_high)

# 11) (Opzionale) salvo i risultati su file CSV nella cartella project_oct25
fwrite(high_counts_by_gene_sample, "project_oct25/high_variants_by_gene_sample.csv")
fwrite(high_counts_by_gene, "project_oct25/high_variants_by_gene_total.csv")
fwrite(data.table(gene = genes_with_high), "project_oct25/genes_with_high_variants.csv")
