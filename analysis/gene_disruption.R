# now that we have the compiled list of where the common insertions are
# what their consensus is, where they land, and what they blast to
# need to find out what's actually being disrupted by these events
    # is there a promoter in them?

library(biomartr)
cds.file <- getCDS(db='genbank', organism = 'Escherichia coli str. K-12 substr. MG1655')
ecoli.cds <- read_cds(file=cds.file)


parse.names <- function(cds){
  n.list <- list()
  for (n in names(cds)){
    n.split <- strsplit(n, split=' \\[')[[1]]
    n.split <- sapply(n.split, function(x) gsub(']','',x))
    names(n.split) <- NULL
    
    gene.name <- n.split[1]
    gene.num <- tail(strsplit(gene.name, '_')[[1]],1)
    n.split <- n.split[2:length(n.split)]
    keys <- sapply(n.split, function(x) strsplit(x, '=')[[1]][1])
    vals <- sapply(n.split, function(x) strsplit(x, '=')[[1]][2])
    n.data <- c(gene.name, gene.num, vals)
    names(n.data) <- c('name', 'num', keys)
    n.list[[gene.name]] <- as.list(n.data)
    }
  
  # bind this all into a big df
  n.df <- bind_rows(n.list)
  # keep only the first 10 cols
  n.df <- n.df[,1:10]
  return(n.df)
}

n.df <- parse.names(ecoli.cds)
write.table(n.df, '~/bhatt_local/meta-mustache/ecoli_cds_table.txt', sep='\t', quote=F, row.names = F, col.names = T)

# for a 
parse.location <- function(l, is.complement=F){
  if (substr(l, 1,11)=='complement('){
    is.complement <- T
    l.sub <- substr(l, 12, nchar(l)-1)
    parsed <- parse.location(l.sub)
    
      
  } else if (substr(l, 1,5)=='join('){
    is.complement = (is.complement & F)
    l.sub <- substr(l, 6, nchar(l)-1)
    l.intervals <- strsplit(l.sub, ',')[[1]]  
  } else {
    is.complement <- (is.complement & F)
    l.split <- strsplit(l, '\\.\\.')[[1]]
  }
}