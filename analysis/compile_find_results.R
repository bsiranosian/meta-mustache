# starting evaluation of the insertion sequence find results
options(stringsAsFactors = F)
library(msa)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(rafalib)
library(gplots)
library(viridis)
# working directory on the cluster 
# only do this to refresh the data
if (F){
  wd <- '/home/ben/scg_lab/meta-mustache'
  setwd(wd)
  
  find.dirs <- list.dirs('find_results',recursive = F)
  find.results.files <- sapply(find.dirs, function(fd) list.files(fd[1], pattern = 'insertions_seqs.tsv', full.names = T))
  find.results.files <- as.character(find.results.files[sapply(find.results.files, length) >0])
  
  find.results.srr <- sapply(find.results.files, function(x) strsplit(x, split='/')[[1]][2])
  names(find.results.files) <- find.results.srr
  
  fr.table <- data.table(read.table(find.results.files[1], sep='\t', quote='', header=T))
  fr.table$srr <- find.results.srr[1]
  
  for (i in 2:length(find.results.files)){
    tmp <- data.table(read.table(find.results.files[i], sep='\t', quote='', header=T))
    tmp$srr <- find.results.srr[i]
    fr.table <- rbind(fr.table,tmp)
  }
  
  write.table(fr.table, '~/bhatt_local/meta-mustache/analysis/find_results_big_table.tsv', 
              sep='\t', quote=F, row.names = F, col.names = T)
}

# normally start here
fr.table <- data.table(read.table('~/bhatt_local/meta-mustache/analysis/find_results_big_table.tsv', sep = '\t', quote='', header=T))
# use this start site
fr.table$use_start <- apply(as.matrix(fr.table[, c('left_site', 'right_site')]), 1, function(x) min(x[x>0]))
m.table <- fr.table[orientation=='M']
ml.table <- fr.table[orientation %in% c('M','L')]
l.table <- fr.table[orientation=='L']
r.table <- fr.table[orientation=='R']

# figure of assembly length distribution
pdf('~/bhatt_local/meta-mustache/analysis/figures/find_results_assembly_length_distribution.pdf', height=4, width=8)
par(mfrow=c(1,3))
hist(m.table$assembly_length, main='Merged sequence assembly length', xlab='Assembly length')
hist(l.table$assembly_length, main='L sequence assembly length', xlab='Assembly length')
hist(r.table$assembly_length, main='R sequence assembly length', xlab='Assembly length')
dev.off()


start.sites <- c(m.table$left_site, l.table$left_site, r.table$right_site)
start.sites.ml <- c(ml.table$left_site)
start.sites.r <- c(r.table$right_site)
leftsite.table <- table(m.table$left_site)
all.site.table <- table(start.sites)
site.table.ml <- table(start.sites.ml)
site.table.r <- table(start.sites.r)
mean.length.ml <- sapply(names(site.table.ml), function(x) mean(ml.table[ml.table$left_site==x, assembly_length]))
mean.length.r <- sapply(names(site.table.r), function(x) mean(r.table[r.table$right_site==x, assembly_length]))

# insertion sequences in 5kb bin
# set up breaks
binsize <- 5000
start <- 0
end <- max(fr.table$use_start)
breaks <- (0:ceiling(end/binsize)) * binsize
bin.starts <- table(cut(fr.table$use_start, breaks=breaks, labels=F))
bin.starts.df <- data.frame(x=as.numeric(names(bin.starts))* binsize, y=as.numeric(bin.starts))
# M only 
bin.starts.m <- table(cut(m.table$use_start, breaks=breaks, labels=F))
bin.starts.m.df <- data.frame(x=as.numeric(names(bin.starts.m))* binsize, y=as.numeric(bin.starts.m))
# L only 
bin.starts.l <- table(cut(m.table$use_start, breaks=breaks, labels=F))
bin.starts.l.df <- data.frame(x=as.numeric(names(bin.starts.l))* binsize, y=as.numeric(bin.starts.l))
# R only 
bin.starts.r <- table(cut(m.table$use_start, breaks=breaks, labels=F))
bin.starts.r.df <- data.frame(x=as.numeric(names(bin.starts.r))* binsize, y=as.numeric(bin.starts.r))
# unique sites only
bin.starts.unique <- table(cut(unique(fr.table$use_start), breaks=breaks, labels=F))
bin.starts.unique.df <- data.frame(x=as.numeric(names(bin.starts.unique))* binsize, y=as.numeric(bin.starts.unique))

ggplot(bin.starts.df, aes(x=x, y=y)) +
  theme_bw() + geom_segment(aes(xend=x), yend=0) + xlab('E. coli genome') + 
  ylab('Number of insertion events') + ggtitle('IS sequences detected in 5kb bin (M+L+R)')
ggsave('~/bhatt_local/meta-mustache/analysis/figures/IS_5kb_bin_MLR.pdf')
ggplot(bin.starts.m.df, aes(x=x, y=y)) +
  theme_bw() + geom_segment(aes(xend=x), yend=0) + xlab('E. coli genome') + 
  ylab('Number of insertion events') + ggtitle('IS sequences detected in 5kb bin (M)')
ggsave('~/bhatt_local/meta-mustache/analysis/figures/IS_5kb_bin_M.pdf')
ggplot(bin.starts.l.df, aes(x=x, y=y)) +
  theme_bw() + geom_segment(aes(xend=x), yend=0) + xlab('E. coli genome') + 
  ylab('Number of insertion events') + ggtitle('IS sequences detected in 5kb bin (L)')
ggsave('~/bhatt_local/meta-mustache/analysis/figures/IS_5kb_bin_L.pdf')
ggplot(bin.starts.r.df, aes(x=x, y=y)) +
  theme_bw() + geom_segment(aes(xend=x), yend=0) + xlab('E. coli genome') + 
  ylab('Number of insertion events') + ggtitle('IS sequences detected in 5kb bin (R)')
ggsave('~/bhatt_local/meta-mustache/analysis/figures/IS_5kb_bin_R.pdf')
# and for unique sites
ggplot(bin.starts.unique.df, aes(x=x, y=y)) +
  theme_bw() + geom_segment(aes(xend=x), yend=0) + xlab('E. coli genome') + 
  ylab('Number of unique insertion events') + ggtitle('IS sequences with unique start detected in 5kb bin (M+L+R)')
ggsave('~/bhatt_local/meta-mustache/analysis/figures/IS_5kb_bin_MLR_unique.pdf')

# consensus sequences for alignments
# minimum of this many 
min.seqs <- 5
start.table <- table(fr.table$use_start)
do.starts <- names(start.table)[start.table >= min.seqs]
site.stats <- as.data.frame(matrix(0, nrow = length(do.starts), ncol=8,
                                   dimnames=list(do.starts,
                                                 c('site', 'mean.length','mean.bits', 'num.seqs', 'consensus',
                                                   'num.m','num.l','num.r'))))
# write out this to be aligned on the cluster
library(seqinr)
outfolder <- '~/scg_lab/meta-mustache/analysis/muscle_is/unaligned_seqs'

if (F){
  for (s in do.starts){
    this.table <- fr.table[use_start==s]
    this.seqs <- this.table[,assembly]
    to.write <- as.list(this.seqs)
    outname <- paste(s, '.fasta', sep='')
    write.fasta(to.write, this.table$srr, file.path(outfolder, outname), as.string = T)
  }
}

# read in aligned seqs, make consensus
infolder <- '~/scg_lab/meta-mustache/analysis/muscle_is/muscle_aligned/'
infiles <- list.files(infolder, full.names = T)
infiles.names <- sapply(list.files(infolder, full.names = F), function(x) strsplit(x, split='\\.')[[1]][1])
infiles.sites <- sapply(infiles.names, function(x) strsplit(x, split='_')[[1]][1])
infiles.orientation <- sapply(infiles.names, function(x) strsplit(x, split='_')[[1]][2])
site.stats.mlr <- as.data.frame(matrix(0, nrow = length(infiles), ncol=9,
                                      dimnames=list(infiles.names,
                                                    c('site', 'mean.length','mean.bits', 'num.seqs', 'consensus',
                                                      'num.m','num.l','num.r','primary.orientation'))))

for (i in 1:length(infiles)){
  print(i)
  seqs <- read.fasta(infiles[i])
  seqs.name <- infiles.names[i]
  seqs.site <- infiles.sites[i]
  seqs.orientation <- infiles.orientation[i]
  
  if (seqs.orientation=='r'){
    this.table <- fr.table[right_site==seqs.site]
  } else {
    this.table <- fr.table[left_site==seqs.site]
  }
  
  seqs.dnastring <- DNAStringSet(sapply(seqs, function(x) DNAString(paste(as.character(x), collapse=''))))
  cm <- consensusMatrix(seqs.dnastring,as.prob = T)[c('A','C','G','T'),]
  cm.all <- consensusMatrix(seqs.dnastring,as.prob = T)
  cs <- consensusString(seqs.dnastring, ambiguityMap='N')
  bits <- apply(cm, 2, function(x) sum(x) * (2+(sum(sapply(x, function(y) y * log(y, base=2)), na.rm = T))))

  site.stats.mlr[seqs.name,'site'] <- seqs.site
  site.stats.mlr[seqs.name,'mean.bits'] <- mean(bits)
  site.stats.mlr[seqs.name,'mean.length'] <- mean(this.table$assembly_length)
  site.stats.mlr[seqs.name,'num.seqs'] <- nrow(this.table)
  site.stats.mlr[seqs.name,'consensus'] <- cs
  site.stats.mlr[seqs.name,'num.m'] <- sum(this.table$orientation=='M')
  site.stats.mlr[seqs.name,'num.l'] <- sum(this.table$orientation=='L')
  site.stats.mlr[seqs.name,'num.r'] <- sum(this.table$orientation=='R')
  orientation.map <- c('M', 'L','R')
  site.stats.mlr[seqs.name,'primary.orientation'] <- orientation.map[which.max(site.stats.mlr[seqs.name, c('num.m','num.l','num.r')])]
}
site.stats.mlr <- as.data.table(site.stats.mlr, keep.rownames=T)

library(ggplot2)
ggplot(site.stats.mlr, aes(y=mean.bits, x=mean.length, col=primary.orientation)) + 
  geom_point() + theme_bw() + xlab('Insertion seqeunce length') + ylab('Mean bitscore') + 
  ggtitle('Insertion sequence consensus alignments')

# keep best alignmnets
bits.thresh <- 1
length.thresh <- 300
best.site.stats <- site.stats.mlr[mean.bits >= bits.thresh & mean.length >= length.thresh]
heatmap.2(t(apply(as.matrix(best.site.stats[,c('num.m','num.l','num.r')]), 1, function(x) x/sum(x))),
          trace='none', Colv=NA, col=rev(redblue(32))[16:32], dendrogram = 'none', cexCol = 1.0,
          main='frequency of orientation')

# write out this to be blasted on the cluster
library(seqinr)
if (F){
  # remove gaps
  to.write <- as.list(sapply(best.site.stats$consensus, function(x) gsub('-', '', x)))
  write.fasta(to.write, names=best.site.stats$site,
              "~/scg_lab/meta-mustache/analysis/blast_consensus_is/best_mlr/best_seqs_mlr.fasta", as.string = T)
}

# read in blast results
blast.results <- data.table(read.table('~/scg_lab/meta-mustache/analysis/blast_consensus_is/best_mlr/diamond_results.tsv', 
                                       sep='\t', quote='', header = T))

name.map <- best.site.stats$site
names(name.map) <- rownames(best.site.stats)
blast.results$qseqid <- name.map[blast.results$qseqid]
# look for some specific sequences
bl.transposase <- blast.results[grep('transposase', blast.results$stitle),]
bl.endonuclease <- blast.results[grep('endonuclease', blast.results$stitle),]
bl.insertion <- blast.results[grep('insertion', blast.results$stitle),]

# add this info into the hist
loc.hist <- hist(rep(as.numeric(site.stats.mlr$site), times=as.numeric(site.stats.mlr$num.seqs)), breaks=1000, plot=F)
loc.hist.best <- hist(rep(as.numeric(best.site.stats$site), times=as.numeric(best.site.stats$num.seqs)), breaks=1000, plot=F)
# loc.hist <- hist(site.stats.mlr$site, breaks=1000, plot = F)
  
col.vector <- rep('black', times=length(loc.hist$counts))
col.vector.best <- rep('black', times=length(loc.hist$counts))
for (i in 1:nrow(bl.transposase)){
  col.vector[which.min(abs(loc.hist$mids - as.numeric(bl.transposase[i,'qseqid'])))] <- 'red'
  col.vector.best[which.min(abs(loc.hist.best$mids - as.numeric(bl.transposase[i,'qseqid'])))] <- 'red'
}
for (i in 1:nrow(bl.endonuclease)){
  if (col.vector[which.min(abs(loc.hist$mids - as.numeric(bl.endonuclease[i,'qseqid'])))]=='black'){
    col.vector[which.min(abs(loc.hist$mids - as.numeric(bl.endonuclease[i,'qseqid'])))] <- 'blue'
  }
  if (col.vector.best[which.min(abs(loc.hist.best$mids - as.numeric(bl.endonuclease[i,'qseqid'])))]=='black'){
    col.vector.best[which.min(abs(loc.hist.best$mids - as.numeric(bl.endonuclease[i,'qseqid'])))] <- 'blue'
  }
}
for (i in 1:nrow(bl.insertion)){
  if (col.vector[which.min(abs(loc.hist$mids - as.numeric(bl.insertion[i,'qseqid'])))]=='black'){
    col.vector[which.min(abs(loc.hist$mids - as.numeric(bl.insertion[i,'qseqid'])))] <- 'green'
  }
  if (col.vector.best[which.min(abs(loc.hist.best$mids - as.numeric(bl.insertion[i,'qseqid'])))]=='black'){
    col.vector.best[which.min(abs(loc.hist.best$mids - as.numeric(bl.insertion[i,'qseqid'])))] <- 'green'
  }
}

mypar(2,1)
barplot(loc.hist$counts, col = col.vector, border = col.vector,
        xlab = 'E. coli genome (4,616 kb)', ylab='IS frequency', main='All detected insertion sequences')
legend(x=700, y=800, legend = c('transposase','endonuclease', 'insertion element','other'),
       fill = c('red', 'blue', 'green','black'), bty = 'n', cex = 0.8)
barplot(loc.hist.best$counts, col = col.vector.best, border = col.vector.best,
        xlab = 'E. coli genome (4,616 kb)', ylab='IS frequency', main='IS > 300bp, >1.0 average bitscore')
legend(x=700, y=400, legend = c('transposase','endonuclease', 'insertion element','other'),
       fill = c('red', 'blue', 'green','black'), bty = 'n', cex = 0.8)

# what are the top results then?
top.counts <- rownames(best.leftsite.stats)[order(best.leftsite.stats$num.seqs, decreasing = T)]
best.leftsite.stats[top.counts,1:3]
blast.results[qseqid==top.counts[1], stitle]

top.bin <- loc.hist.best$mids[which.max(loc.hist.best$counts)]
top.site.info <- best.site.stats[abs(as.numeric(site)-top.bin) <2500,]
# all diguanylate cyclase?
blast.results[qseqid==top.site.info$site[1],]
blast.results[qseqid==top.site.info$site[2],]

bins.sorted <- loc.hist.best$mids[order(loc.hist.best$counts, decreasing = T)][1:sum(loc.hist.best$counts>0)]
# top.site.info <- best.site.stats[abs(as.numeric(site)-top.bin) <2500,]

# instead of working with bins, lets just work from the sequences directly
# theres a few duplicates 
dup.sites <- best.site.stats[best.site.stats$site %in% best.site.stats[duplicated(best.site.stats$site),site],]
dup.sites[,c(1:4, 6:9)]

best.unique <- best.site.stats[!(site %in% dup.sites[,site])]
best.unique$top.hit <- sapply(best.unique$site, function(x) blast.results[qseqid==x, stitle][1])
best.unique$top.hit.evalue <- sapply(best.unique$site, function(x) blast.results[qseqid==x, evalue][1])
best.unique <- best.unique[order(best.unique$num.seqs, decreasing = T), ]


# read in the MERGEM database results
mergm.results <- data.table(read.table('~/scg_lab/meta-mustache/analysis/blast_consensus_is/best_mlr/blast_mergem_results.tsv', 
                                       sep='\t', quote='', header=T))
best.unique$top.hit.mergem <- sapply(best.unique$site, function(x) mergm.results[qseqid==x, stitle][1])
best.unique$top.hit.mergem.evalue <- sapply(best.unique$site, function(x) mergm.results[qseqid==x, evalue][1])
best.unique$top.hit.mergem.family <- sapply(best.unique$top.hit.mergem, function(x) strsplit(x, split=' ')[[1]][3])
write.table(best.unique, '~/bhatt_local/meta-mustache/analysis/best_mlr_unique_table.tsv', sep='\t', quote=F, row.names = F, col.names = T)

# MUCH BETTER PLOT
ggplot(best.unique, aes(x=as.numeric(site), y=as.numeric(num.seqs), col=top.hit.mergem.family)) +
  theme_bw() + geom_segment(aes(xend=as.numeric(site)), yend=0) + xlab('E. coli genome') + 
  ylab('Number of insertion events') + ggtitle('IS sequences, MERGEM database results')

# see where they land (in genes?) 








