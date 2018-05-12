library(ggtree)
library(treeio)
library(data.table)
options(stringsAsFactors = F)

tree.file <- "~/scg/software/nsegata-phylophlan-1d174e34b2ae/output/multiprotein/multiprotein.tree.nwk"
nwk <- system.file("extdata", "~/scg/software/nsegata-phylophlan-1d174e34b2ae/output/multiprotein/multiprotein.tree.nwk", package="treeio")
tree <- read.tree(tree.file)


resistance.table <- data.table(read.table('~/bhatt_local/meta-mustache/resistance_table_ecoli.tsv', sep='\t', quote="",header = T))
res.matrix.txt <- read.table('~/bhatt_local/meta-mustache/resistance_matrix_text.tsv', sep='\t', quote="",header = T, na.strings = 'NA')
res.matrix.txt <- cbind(rownames(res.matrix.txt), res.matrix.txt)
colnames(res.matrix.txt)[1] <- 'taxa'
res.matrix.numeric <- as.matrix(read.table('~/bhatt_local/meta-mustache/resistance_matrix_numeric.tsv', sep='\t', quote="",header = T))
uab <- colnames(res.matrix.txt)
usrr <- rownames(res.matrix.txt)

for (ab in uab){
  
  p <- ggtree(tree, branch.length = 'none', layout='circular')
  p <- p %<+% res.matrix.txt + geom_tippoint(aes_string(color=ab))
  p + theme(legend.position = 'right')
  ggsave(paste('~/bhatt_local/meta-mustache/analysis/figures/phylo_tree_colored_ab/',ab,'.png', sep=''))
}

  p <- p %<+% res.matrix.txt + geom_tippoint(aes(color=amikacin))
p + theme(legend.position = 'right')

p <- ggtree(tree, branch.length = 'none', layout='circular')
p <- p %<+% as.data.frame(res.matrix.numeric) + geom_tippoint(aes(size=ampicillin))
p + theme(legend.position = 'right')

