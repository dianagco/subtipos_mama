library(ggplot2)

start <- Sys.time()
cat("Starting at, ", format(start, "%H:%M:%OS3"), "\n" )

setwd("~/Workspace/INMEGEN/subtipos_mama/")

cat("Loading data\n")
load(file="lista_distancias_intracromosomales.RData")
conds <- c("basal")

for (cond in conds) {
  
  cat("Working with condition ", cond, "\n")
  load(file=paste("adjmatrix_", cond, ".RData", sep=""))
  
  ch <- as.character(1)
  
  distancias <- lista_distmatrices[[ch]]
  genes <- rownames(distancias)
  r <- rep(genes, each = length(genes))
  c <- rep(genes, times = length(genes))
  dist.df <- data.frame(source = r, target = c, distance = unlist(distancias),  row.names = NULL)
  dist.df <- dist.df[!is.na(dist.df$distance), ]
  dist.df <- dist.df[dist.df$distance > 0, ]
  
  mis <- adjmatrix[rownames(adjmatrix) %in% genes, colnames(adjmatrix) %in% genes]
  genes <- rownames(mis)
  r <- rep(genes, each = length(genes))
  c <- rep(genes, times = length(genes))
  mi.df <- data.frame(source = r, target = c, mi = unlist(mis), row.names = NULL)
  mi.df <- mi.df[!is.na(mi.df$mi), ]
  
  dist.mi.df <- merge(dist.df, mi.df)
  dist.mi.df <- dist.mi.df[with(dist.mi.df, order(distance, mi)), ]
  
  cat("Valor máximo duplicados: ", max(table(dist.mi.df[duplicated(dist.mi.df$distance), "distance"])), "\t")
  
  write.table(dist.mi.df, file = paste(cond,"_mi_distance_", ch, ".tsv", sep=""), sep="\t", 
              col.names = T, row.names = F, quote = F)
  
}