library(SNPRelate)
library(SeqArray)
library(ggtree)

vcf.f <- '~/bhatt_local/freebayes_all_10k.vcf'
vcf.o <- '~/bhatt_local/freebayes_all_10k.gds'
seqVCF2GDS(vcf.f, vcf.o)

genofile <- seqOpen(vcf.o)

snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2, autosome.only = F)
snpset.id <- unlist(snpset)

pca <- snpgdsPCA(genofile, autosome.only = FALSE, num.thread=1, snp.id=snpset.id)
pca <- snpgdsPCA(genofile, autosome.only = FALSE, num.thread=1, snp.id = NULL)

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
plot(tab$EV1, tab$EV2, xlab="eigenvector 1", ylab="eigenvector 2")


dissMatrix  <-  snpgdsDiss(genofile, sample.id=NULL, snp.id=NULL, 
                           autosome.only=FALSE,remove.monosnp=TRUE, maf=NaN, 
                           missing.rate=NaN, num.thread=4, verbose=TRUE)

snpHCluster <-  snpgdsHCluster(dissMatrix, sample.id=NULL, need.mat=TRUE, hang=0.25)

cutTree <- snpgdsCutTree(snpHCluster, z.threshold=15, n.perm = 5000,
                         samp.group=NULL,outlier.n = 1,
                         col.list=NULL, pch.outlier=4,
                         pch.list=NULL,label.H=FALSE, label.Z=FALSE, verbose=TRUE)


snpgdsDrawTree(cutTree, main = "E. coli isolates", edgePar=list(col=rgb(0.5,0.5,0.5,0.75),t.col="black"),
               y.label.kinship=F,leaflab="perpendicular", cex=0.5)



dm <- as.dist(dissMatrix$diss)
hc <- hclust(dm,method='average')
p <- ggtree(as.phylo(hc), branch.length = 'none', layout='circular')
p + theme(legend.position = 'right')  
