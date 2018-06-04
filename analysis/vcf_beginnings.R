library(SNPRelate)
library(SeqArray)
library(ggtree)
# setwd('~/bhatt_local/meta-mustache/vcf/')
# vcf.f <- 'freebayes_all_10k.vcf'
# vcf.o <- 'freebayes_all_10k.gds'
# ntasks <- 4

setwd('/labs/asbhatt/bsiranos/meta-mustache_nomobile/freebayes_simul')
vcf.f <- 'freebayes_all_snps.vcf'
vcf.o <- 'freebayes_all_snps.gds'
vcf.ef.f <- 'e_fergusonii.vcf'
vcf.ef.o <- 'e_fergusonii.gds'
vcf.both.o <- 'ecoli_fergusonii_snps.gds'
ntasks <- 32
# seqVCF2GDS(vcf.f, vcf.o, parallel=ntasks)
# seqVCF2GDS(vcf.ef.f, vcf.ef.o, parallel=ntasks)

seqMerge(c(vcf.o, vcf.ef.o),vcf.both.o)

genofile <- seqOpen(vcf.o)
genofile.ef <- seqOpen(vcf.ef.o)

snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2, autosome.only = F, num.thread=ntasks, slide.max.bp=50000)
snpset.id <- unlist(snpset)
saveRDS(snpset, 'freebayes_LD_snpset_50k.rds')

pca <- snpgdsPCA(genofile, autosome.only = FALSE, num.thread=1, snp.id=snpset.id)
# pca <- snpgdsPCA(genofile, autosome.only = FALSE, num.thread=ntasks, snp.id = NULL)

# variance proportion (%)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
# make a data.frame
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],    
                  stringsAsFactors = FALSE)
head(tab)
pdf('pca.pdf')
plot(tab$EV1, tab$EV2, xlab="eigenvector 1", ylab="eigenvector 2")
dev.off()

dissMatrix  <-  snpgdsDiss(genofile, sample.id=NULL, snp.id=snpset.id, 
                           autosome.only=FALSE,remove.monosnp=TRUE, maf=NaN, 
                             missing.rate=NaN, num.thread=1, verbose=TRUE)

snpHCluster <-  snpgdsHCluster(dissMatrix, sample.id=NULL, need.mat=TRUE, hang=0.25)

cutTree <- snpgdsCutTree(snpHCluster, z.threshold=15, n.perm = 5000,
                         samp.group=NULL,outlier.n = 1,
                         col.list=NULL, pch.outlier=4,
                         pch.list=NULL,label.H=FALSE, label.Z=FALSE, verbose=TRUE)

pdf('tree.pdf', height=8, width=40)
snpgdsDrawTree(cutTree, main = "E. coli isolates", edgePar=list(col=rgb(0.5,0.5,0.5,0.75),t.col="black"),
               y.label.kinship=F,leaflab="perpendicular", cex.lab=0.5)
dev.off()


dm <- as.dist(dissMatrix$diss)
hc <- hclust(dm,method='average')
p <- ggtree(as.phylo(hc), branch.length = 'none', layout='circular')
p + theme(legend.position = 'right')  
