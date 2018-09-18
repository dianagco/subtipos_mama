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
  
  #write.table(dist.mi.df, file = paste(cond,"_mi_distance_", ch, ".tsv", sep=""), sep="\t", 
  #            col.names = T, row.names = F, quote = F)
  
  bins <- split(dist.mi.df, seq(from =1, to = nrow(dist.mi.df)-1)%/%50)
  testy <- parallel::mclapply(X = bins, mc.cores = 3,  mc.cleanup = FALSE, FUN = function(b){
    myplot <- ggplot(b, aes(mi))+geom_density()
    p <- layer_data(myplot)
    return(matrix(c(p[which.max(p$y), "x"], mean(b$mi), median(b$mi), max(b$mi)), ncol = 4))
  })
  stats <- plyr::ldply(testy)
  colnames(stats) <- c("id", "moda", "media", "mediana", "max")
  write.table(stats, file = paste("stats_", cond,"_100k_", ch, ".tsv", sep=""), sep="\t", 
              col.names = T, row.names = F, quote = F)
  myplot <- ggplot(data=stats, aes(x=id+1, y=max, group=1)) +
    geom_point(size = 1.5)
  print(myplot)
  
  plot(stats$id +1 , stats$media, main="Scatterplot Example", 
       xlab="id", ylab="media", pch=19)
  abline(lm(stats$id +1~stats$media), col="red") # regression line (y~x) 
  lines(lowess(stats$media,stats$id), col="blue") 
}

