library(ggplot2)

start <- Sys.time()
cat("Starting at, ", format(start, "%H:%M:%OS3"), "\n" )

setwd("~/Workspace/subtipos_mama")

cat("Loading data\n")
load(file="genes.annot.RData")
load(file="adjmatrix_basal.RData")


all.genes.annot <- genes.annot
chrs <- c(as.character(1:22), "X")
binlength <- 100000

## Extract intra-chromosmal interactions
for (ch in chrs) {
  cat("Working with chromosome", ch, "\n")
  genes.annot <- all.genes.annot[all.genes.annot$Chr == ch, ]
  genes <- genes.annot$symbol
  genes <- genes[1:10]
  g <- parallel::mclapply(X = genes, mc.cores = 4,  mc.cleanup = FALSE, FUN = function(gene1) {
    #cat("Gene: ", gene1, "\n")
    other.genes <- genes[which(genes != gene1)]
    
    mivals <- lapply(X = other.genes, 
                                FUN = function(gene2) {
                                  distance <- max(genes.annot[gene2, "start"], genes.annot[gene1, "start"]) - min(genes.annot[gene1, "end"], genes.annot[gene2, "end"])
                                  bin <- max(floor(distance/binlength), 1)
                                  list(bin, max(adjmatrix[gene1, gene2], adjmatrix[gene2, gene1], na.rm = T))
                                })
    matrix(unlist(mivals), ncol = 2, byrow = T)
  })
  cat("Getting bin means.\n")
  df <- plyr::ldply(g)
  colnames(df) <- c("bin", "mi")
  addtodf <- data.frame("bin" = 2, "mi"= c(0.7, 0.75, 0.5) )
  df <- rbind(df, addtodf)
  
  testy <- lapply(unique(df$bin), function(b) {
    l <- df[df$bin == b, ]
    l$order <- order(l$mi)
    l
  })

  colnames(dfmeans) <- c("distance", "avg.mi") 
  dfmeans$distance <- dfmeans$distance*binlength
  cat("Saving plot.\n")
  png(paste("Basal_ch", ch, ".png", sep=""), width =800 , height = 400)
  myplot <- ggplot(data=dfmeans, aes(x=distance, y=avg.mi, group=1)) +
    geom_point(size = 1.5)
  print(myplot)
  dev.off() 
  cat("Saving data.\n")
  save(dfmeans, file=paste("basal_binmeans_ch", ch, ".RData", sep=""))
 
}
end <- Sys.time()
cat("Whole thing took ", (end - start), "\n")
cat("Ending at, ", format(end, "%H:%M:%OS3"), "\n" )

