library(data.table)
options(stringsAsFactors = F)
# compile a large table of resistance values
# that can be made into a heatmap
resistance.table <- data.table(read.table('~/bhatt_local/meta-mustache/klebsiella_resistance_table.tsv', sep='\t', quote="",header = T))
# make measurement simple
resistance.table$measurement_simple <- sapply(resistance.table$Measurement, function(x) mean(as.numeric(strsplit(x, split='/')[[1]])))
uab <- unique(resistance.table$Antibiotic)
usrr <- unique(resistance.table$sample_accession)
res.matrix.txt <- matrix("", nrow=length(usrr), ncol=length(uab), dimnames = list(usrr, uab))
res.matrix.numeric <- matrix(0,  nrow=length(usrr), ncol=length(uab), dimnames = list(usrr, uab))
for (ab in uab){
  for (s in usrr){
    val <- resistance.table[which(Antibiotic==ab & sample_accession==s), Resistance.phenotype]
    val.numeric <- resistance.table[which(Antibiotic==ab & sample_accession==s), measurement_simple]
    if (length(val)==0){
      val <- NA
    } else if (length(val) >1){
      val <- paste(val, collapse = '_')
    }
    if (length(val.numeric)==0){
      val.numeric <- NA
    } else if (length(val.numeric) >1){
      val.numeric <- mean(val.numeric)
    }
    # print(val)
    res.matrix.txt[s, ab] <- val
    res.matrix.numeric[s, ab] <- val.numeric
  }
}
# fix some text issues 
res.matrix.txt[res.matrix.txt=='not defined'] <- NA
res.matrix.txt[res.matrix.txt=='not defined_not defined'] <- NA
res.matrix.txt[res.matrix.txt=='resistant_resistant'] <- 'resistant'

res.matrix.txt[res.matrix.txt=='susceptible_susceptible'] <- 'susceptible'
res.matrix.txt[res.matrix.txt=='resistant_susceptible'] <- 'intermediate'
res.matrix.txt[res.matrix.txt=='susceptible_resistant'] <- 'intermediate'
res.matrix.txt[res.matrix.txt=='resistant_intermediate'] <- 'intermediate'
res.matrix.txt[res.matrix.txt=='susceptible_intermediate'] <- 'susceptible'
res.matrix.txt[res.matrix.txt=='susceptible-dose dependent'] <- 'intermediate'

table(res.matrix.txt, useNA = 'ifany')
write.table(res.matrix.txt, '~/bhatt_local/meta-mustache/klebsiella_resistance_matrix_text.tsv', sep='\t', quote=F, row.names = T, col.names = T)
write.table(res.matrix.numeric, '~/bhatt_local/meta-mustache/klebsiella_resistance_matrix_numeric.tsv', sep='\t', quote=F, row.names = T, col.names = T)

# fraction of tested antibiotics susceptible
fs <- apply(res.matrix.txt, 1, function(x) sum(x=='susceptible', na.rm = T) / sum(!is.na(x)))

best.s <- names(sort(fs, decreasing = T))
best.s.table <- data.frame(accession=best.s,
                           num_susceptable=sapply(best.s, function(x) sum(res.matrix.txt[x,]=='susceptible', na.rm = T)),
                           num_tested=sapply(best.s, function(x) sum(!(is.na(res.matrix.txt[x, ])))))

kg <- read.table('~/klebsiella_get/kleb_genomes.txt', sep='\t', quote='', header=T)
which(best.s.table$accession %in% kg$BioSample)
head(best.s.table[which(best.s.table$accession %in% kg$BioSample), ], 10)
