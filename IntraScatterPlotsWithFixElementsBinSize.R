library(ggplot2)
library(dplyr)
library(tibble)

start <- Sys.time()
cat("Starting at, ", format(start, "%H:%M:%OS3"), "\n" )

setwd("/labs/csbig/subtipos_mama")

cat("Loading data\n")
load(file="genes.annot.RData")
load(file="adjmatrix_basal.RData")

ch <- as.character(1)
all.genes.annot <- genes.annot
genes.annot <- all.genes.annot[all.genes.annot$Chr == ch, ]
genes <- genes.annot$symbol
g <- parallel::mclapply(X = genes, mc.cores = 6,  mc.cleanup = FALSE, FUN = function(gene1) {
  #cat("Gene: ", gene1, "\n")
  mivals <- lapply(X = genes, 
                   FUN = function(gene2) {
                     distance <- max(max(genes.annot[gene2, "start"], genes.annot[gene1, "start"]) - min(genes.annot[gene1, "end"], genes.annot[gene2, "end"]), 0)
                   })
  names(mivals) <- genes
  unlist(mivals)
})

distancias  <- plyr::ldply(g)
rownames(distancias) <- colnames(distancias)

distancias <- lista_distmatrices[[ch]]
genes <- rownames(distancias)
r <- rep(genes, each = length(genes))
c <- rep(genes, times = length(genes))
dist.df <- data.frame(source = r, target = c, distance = unlist(distancias),  row.names = NULL)
dist.df <- dist.df[!is.na(dist.df$distance), ]

mis <- adjmatrix[rownames(adjmatrix) %in% genes, colnames(adjmatrix) %in% genes]
genes <- rownames(mis)
r <- rep(rownames(dist.df), each = length(genes))
c <- rep(rownames(dist.df), times = length(genes))
mi.df <- data.frame(source = r, target = c, mi = unlist(mis), row.names = NULL)
mi.df <- mi.df[!is.na(mi.df$mi), ]

dist.mi.df <- merge(dist.df, mi.df)
dist.mi.df <- dist.mi.df[with(dist.mi.df, order(distance, mi)), ]
