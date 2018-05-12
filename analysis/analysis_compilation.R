# getting all the data together for meta-mustache
options(stringsAsFactors = F)
library(data.table)
library(rafalib)
# working directory on the cluster 
wd <- '/home/ben/scg_lab/meta-mustache'
setwd(wd)

srr_list.file <- 'srr_download_list_all_ecoli.txt'
srr_list <- read.table(srr_list.file)[,1]

# data table of resistance phenotypes
resistance.table.file <- '/home/ben/bhatt_local/meta-mustache/resistance_table.tsv'
resistance.table <- data.table(read.table(resistance.table.file, sep='\t', header = T, fill=T))
resistance.table.orig <- data.table(read.table(resistance.table.file, sep='\t', header = T, fill=T))
# mapping from id in table to SRR id
id_to_srr.file <- '/home/ben/bhatt_local/meta-mustache/ncbi_webpages/id_to_srr.txt' 
id_to_srr <- read.table(id_to_srr.file, sep='\t', header=F, fill=T)
rownames(id_to_srr) <- id_to_srr[,1]
colnames(id_to_srr) <- c('id','srr')
id_to_srr <- id_to_srr[(id_to_srr[,2] %in% srr_list),]

# get SRR ids in the resistance table
resistance.table$srr <- id_to_srr[as.character(resistance.table$sample_id),2]
# remove NA rows
resistance.table <- resistance.table[!is.na(resistance.table$srr),]
# set measurement_simple for combined treatments
resistance.table$measurement_simple <- sapply(resistance.table$measurement, function(x) mean(as.numeric(strsplit(x, split='/')[[1]])))
# write.table(resistance.table, '~/bhatt_local/meta-mustache/resistance_table_ecoli.tsv', sep='\t', quote=F, row.names = F, col.names = T)
# write.table(resistance.table, '~/scg_lab/meta-mustache/analysis/resistance_table_ecoli.tsv', sep='\t', quote=F, row.names = F, col.names = T)

# plot of density of resistance values for each antibiotic
unique_ab <- names(sort(table(resistance.table$antibiotic))[sort(table(resistance.table$antibiotic))>1])
pdf('~/scg_lab/meta-mustache/analysis/resistance_density.pdf', height=20, width=10)
mypar(10,4, cex.lab=0.9, cex.main=0.9)
for (ab in rev(unique_ab)){
  ab_table <- table(resistance.table[antibiotic==ab, resistance_phenotype])
  ab_table_text <- sapply(1:length(ab_table), function(i) paste(names(ab_table)[i], ':  ', ab_table[i], sep=''))
  d <- density(resistance.table[antibiotic==ab, measurement_simple])
  plot(d, xlab='Measurement (mg/L)',
       main = paste(ab, 'resistance values', '\n', paste(ab_table_text, collapse='  ')))
  polygon(d, col='grey70')
}
dev.off()

# add some information to the resistance table
# number of reads in each after trimming
multiqc_data.file <- '~/scg_lab/meta-mustache/multiqc_data/multiqc_fastqc.txt'
multiqc_data <- read.table(multiqc_data.file, sep='\t', quote='', header=T)
multiqc_data$srr <- sapply(multiqc_data$Sample, function(x) strsplit(x, split='_')[[1]][1])
srr_reads <- data.frame(multiqc_data[!duplicated(multiqc_data$srr), c('srr', 'Total.Sequences')])
rownames(srr_reads) <- srr_reads$srr
resistance.table$reads <- srr_reads[resistance.table$srr, "Total.Sequences"]

# number of insertion seqs found
is_count.file <- '~/scg_lab/meta-mustache/analysis/insertion_seqs_wc.txt'
is_count <- read.table(is_count.file, sep=' ', header=F, quote='')
colnames(is_count) <- c('filename', 'lines')
rownames(is_count) <- sapply(is_count$filename, function(x) strsplit(x, split='/')[[1]][7])

# number of insertion seqs sent to blast
isb_count.file <- '~/scg_lab/meta-mustache/analysis/insertion_seqs_blast_wc.txt'
isb_count <- read.table(isb_count.file, sep=' ', header=F, quote='')
colnames(isb_count) <- c('filename', 'lines')
rownames(isb_count) <- sapply(isb_count$filename, function(x) strsplit(x, split='/')[[1]][7])

# how do number of reads and number of insertion sequences compare
plot(srr_reads[rownames(is_count),2], is_count[,2],
     xlab='Total reads', ylab='IS found',
     main=paste('Reads vs. IS\nspearman: ',
                round(cor(srr_reads[rownames(is_count),2], is_count[,2], method ='spear'), digits = 3)))
# normalize IS per million reads
norm_is_count <- sapply(rownames(is_count), function(x) is_count[x, "lines"]/
                          (srr_reads[x, "Total.Sequences"]/1000000))

# fraction of tested antibiotics resistant
fraction_resistant <- sapply(unique(resistance.table$srr), function(x) {
  a <- resistance.table[srr==x, resistance_phenotype]
  sum(a=='resistant')/length(a)
})

plot(fraction_resistant[names(norm_is_count)], norm_is_count,
     xlab='resistant fraction', ylab='norm IS found',
     main=paste('resistance vs. IS\nspearman: ',
                round(cor(fraction_resistant[names(norm_is_count)], norm_is_count, method ='spear'), digits = 3)))

# write.table(resistance.table, '~/scg_lab/meta-mustache/analysis/resistance_table_ecoli.tsv', sep='\t', quote=F, row.names = F, col.names = T)
# write.table(resistance.table, '~/bhatt_local/meta-mustache/resistance_table_ecoli.tsv', sep='\t', quote=F, row.names = F, col.names = T)

# compile a large table of resistance values
# that can be made into a heatmap
resistance.table <- data.table(read.table('~/bhatt_local/meta-mustache/resistance_table_ecoli.tsv', sep='\t', quote="",header = T))
uab <- unique(resistance.table$antibiotic)
usrr <- unique(resistance.table$srr)
res.matrix.txt <- matrix("", nrow=length(usrr), ncol=length(uab), dimnames = list(usrr, uab))
res.matrix.numeric <- matrix(0,  nrow=length(usrr), ncol=length(uab), dimnames = list(usrr, uab))
for (ab in uab){
  for (s in usrr){
    val <- resistance.table[which(antibiotic==ab & srr==s), resistance_phenotype]
    val.numeric <- resistance.table[which(antibiotic==ab & srr==s), measurement_simple]
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
res.matrix.txt[res.matrix.txt=='resistant_resistant'] <- 'resistant'
res.matrix.txt[res.matrix.txt=='susceptible_susceptible'] <- 'susceptible'
res.matrix.txt[res.matrix.txt=='resistant_susceptible'] <- 'intermediate'
res.matrix.txt[res.matrix.txt=='susceptible_resistant'] <- 'intermediate'
res.matrix.txt[res.matrix.txt=='resistant_intermediate'] <- 'intermediate'
res.matrix.txt[res.matrix.txt=='susceptible-dose dependent'] <- 'intermediate'

table(res.matrix.txt, useNA = 'ifany')
write.table(res.matrix.txt, '~/bhatt_local/meta-mustache/resistance_matrix_text.tsv', sep='\t', quote=F, row.names = T, col.names = T)
write.table(res.matrix.numeric, '~/bhatt_local/meta-mustache/resistance_matrix_numeric.tsv', sep='\t', quote=F, row.names = T, col.names = T)

# gross heatmap of this data
library(gplots)
library(ggplot2)
heatmap.2(matrix(as.numeric(as.factor(res.matrix)), nrow=nrow(res.matrix), ncol=ncol(res.matrix)), Rowv = NA, Colv=NA, trace='none')

ggplot(data.frame(as.factor(res.matrix))) + geom_tile(aes(fill = value),
                                                colour = "white") + scale_fill_manual(values=c("red", "blue", "black"))

ggplot(resistance.table, aes(antibiotic, srr)) + geom_tile(aes(fill = resistance_phenotype),
                                                colour = "white") + scale_fill_manual(values=c("red", "blue", "black", "green", "purple"))

       