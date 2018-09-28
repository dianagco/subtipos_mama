library(ggplot2)
start <- Sys.time()
cat("Starting at, ", format(start, "%H:%M:%OS3"), "\n" )

setwd("~/Workspace/subtipos_mama/")

cat("Loading data\n")
cytobands <- read.delim(file="cytoBand.txt", header = F)
colnames(cytobands) <- c("chr", "start", "end", "cytoname", "gpos")
load("rdata/genes.annot.RData")
load(file="rdata/adjmatrix_basal.RData")

chrs <- c(as.character(1:22), "X", "Y")
chrs <- c(as.character(1))
conds <- c("basal", "healthy")

all.genes.annot <- genes.annot
for (ch in chrs) {
  genes.annot <- all.genes.annot[all.genes.annot$Chr == ch, ]
  genes <- genes.annot$symbol
  g <- parallel::mclapply(X = genes, mc.cores = 6,  mc.cleanup = FALSE, FUN = function(gene1) {
    #cat("Gene: ", gene1, "\n")
    other.genes <- genes[which(genes != gene1)]
    
    mivals <- lapply(X = other.genes, 
                     FUN = function(gene2) {
                       cito.gene1 <- cytobands[cytobands$chr == paste("chr", ch, sep = "") &  
                                                 (genes.annot[gene1, "start"] >= cytobands$start 
                                                  & genes.annot[gene1, "start"] <= cytobands$end), ]
                       cito.gene2 <- cytobands[cytobands$chr == paste("chr", ch, sep = "") &  
                                                 (genes.annot[gene2, "start"] >= cytobands$start 
                                                  & genes.annot[gene2, "start"] <= cytobands$end), ]
                       
                       cito.distance <- abs(as.numeric(rownames(cito.gene1)) - as.numeric(rownames(cito.gene2)))
                       list(cito.distance, max(adjmatrix[gene1, gene2], adjmatrix[gene2, gene1], na.rm = T))
                      })
    matrix(unlist(mivals), ncol = 2, byrow = T)
  })
  cat("Getting bin means.\n")
  df <- plyr::ldply(g)
  colnames(df) <- c("cito.distance", "mi")
  dfmeans <- aggregate(mi ~ cito.distance, data=df, mean, na.rm=TRUE)
  colnames(dfmeans) <- c("cito.distance", "avg.mi") 
  cat("Saving plot.\n")
  png(paste("Basal_ch", ch, ".png", sep=""), width =800 , height = 400)
  myplot <- ggplot(data=df, aes(x=cito.distance, y=mi, group=1)) +
    geom_point(size = 1.5)
  print(myplot)
  dev.off() 
  cat("Saving data.\n")
}

