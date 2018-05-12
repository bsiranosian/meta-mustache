# testing for differential resistance in any of these groups that 
# have insertion sequences at the same position
library(data.table)
library(gplots)
library(viridis)
library(grDevices)

resistance.table <- data.table(read.table('~/bhatt_local/meta-mustache/resistance_table_ecoli.tsv', sep='\t', quote="",header = T))
uab <- unique(resistance.table$antibiotic)
usrr <- unique(resistance.table$srr)
fr.table <- data.table(read.table('~/bhatt_local/meta-mustache/find_results_big_table.tsv', sep = '\t', quote='', header=T))
unique.srr <- unique(fr.table$srr)
# best consensus insertions
best.is <- data.table(read.table('~/bhatt_local/meta-mustache/analysis/best_mlr_unique_table.tsv', sep='\t', quote='', header = T))

# limit resistance table to the unique srr
resistance.table <- resistance.table[srr %in% unique.srr, ]


# For each event, get srr that had that event and not
srr.with.list <- list()
srr.without.list <- list()
for (i in 1:nrow(best.is)){
  this.is <- best.is[i, ]
  l.table <- fr.table[left_site==this.is$site, ]
  r.table <- fr.table[right_site==this.is$site, ]
  srr.with.event <- unique(c(l.table$srr, r.table$srr))
  srr.without.event <- unique.srr[!(unique.srr %in% srr.with.event)]
  srr.with.list[[as.character(this.is$site)]] <- list(srr.with.event)
  srr.without.list[[as.character(this.is$site)]] <- list(srr.without.event)
}


ab.ttest.table <- matrix(list(), nrow=nrow(best.is), ncol = length(uab),
                         dimnames=list(best.is$site, uab))
ab.ttest.pv <- matrix(NA, nrow=nrow(best.is), ncol = length(uab),
                         dimnames=list(best.is$site, uab))

for (site in best.is$site){
  site <- as.character(site)
  srr.with.event <- srr.with.list[[site]]
  srr.without.event <- srr.without.list[[site]]
  for (ab in uab){
    ab.table <- resistance.table[antibiotic==ab, ]  
    with.measurement <- na.omit(ab.table[srr.with.event, measurement_simple, on='srr'])
    without.measurement <- na.omit(ab.table[srr.without.event, measurement_simple, on='srr'])
    if (length(with.measurement)>=5 & length(without.measurement)>=5){
      t.res <- t.test(with.measurement, without.measurement)
      t.pv <- t.res$p.value
    } else{
      t.pv =1 
    }
    ab.ttest.table[site, ab] <- list(t.res)
    ab.ttest.pv[site, ab] <- t.pv
  }
}

heatmap.2(ab.ttest.pv, trace='none', Rowv=NA, Colv=NA, dendrogram='none', col=rev(viridis(32)))
heatmap.2(-log(ab.ttest.pv), trace='none', col=(viridis(32)), margins = c(5,5))

# some good results here?
ab <- 'ceftriaxone'
test.sites <- names(head(sort(log(ab.ttest.pv[,ab]))))
best.is[as.character(site) %in% test.sites, ]

site <- test.sites[1]

ab.table <- resistance.table[antibiotic==ab, ]
srr.with.event <- srr.with.list[[site]]
srr.without.event <- srr.without.list[[site]]
with.measurement <- na.omit(ab.table[srr.with.event, measurement_simple, on='srr'])
without.measurement <- na.omit(ab.table[srr.without.event, measurement_simple, on='srr'])
d.all <- density(ab.table$measurement_simple)
d.with <- density(with.measurement)
d.without <- density(without.measurement)

plot(d, ylim=c(0, max(c(d.all$y, d.with$y, d$without$y))))
polygon(d, col=adjustcolor('grey60', alpha=0.5))
lines(d.with)
polygon(d.with, col=adjustcolor('firebrick', alpha=0.5))
lines(d.without)
polygon(d.without, col=adjustcolor('steelblue', alpha=0.5))

# top 20 hits
mat.ord <- order(log(ab.ttest.pv))
pdf('~/bhatt_local/meta-mustache/analysis/atibiotic_IS_ttest.pdf', height=15, width=10)
par(mfrow=c(5,4))
for (i in 1:100){
  k <- arrayInd(mat.ord[i], dim(ab.ttest.pv))
  site <- rownames(ab.ttest.pv)[k[,1]]
  ab <- colnames(ab.ttest.pv)[k[,2]]
  pv <- ab.ttest.pv[k]
  ab.table <- resistance.table[antibiotic==ab, ]
  srr.with.event <- srr.with.list[[site]]
  srr.without.event <- srr.without.list[[site]]
  with.measurement <- na.omit(ab.table[srr.with.event, measurement_simple, on='srr'])
  without.measurement <- na.omit(ab.table[srr.without.event, measurement_simple, on='srr'])
  d.all <- density(ab.table$measurement_simple)
  d.with <- density(with.measurement)
  d.without <- density(without.measurement)
  
  title <- paste(ab, site, signif(pv, 3), sep=' | ')
  
  plot(d, ylim=c(0, max(c(d.all$y, d.with$y, d$without$y))), 
       main=title, col='white',cex.main=0.75)
  # polygon(d, col=adjustcolor('grey60', alpha=0.5))
  lines(d.with)
  polygon(d.with, col=adjustcolor('firebrick', alpha=0.5))
  lines(d.without)
  polygon(d.without, col=adjustcolor('steelblue', alpha=0.5))
  legend('topleft', legend=c('with IS', 'without IS'), fill=c(adjustcolor('firebrick', alpha=0.5), adjustcolor('steelblue', alpha=0.5)), bty = 'n')
}
dev.off()


for (ab in rev(unique_ab)){
  ab_table <- table(resistance_table[antibiotic==ab, resistance_phenotype])
  ab_table_text <- sapply(1:length(ab_table), function(i) paste(names(ab_table)[i], ':  ', ab_table[i], sep=''))
  d <- density(resistance_table[antibiotic==ab, measurement_simple])
  plot(d, xlab='Measurement (mg/L)',
       main = paste(ab, 'resistance values', '\n', paste(ab_table_text, collapse='  ')))
  polygon(d, col='grey70')
}
dev.off()