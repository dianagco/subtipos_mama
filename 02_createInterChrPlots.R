library(ggplot2)
library(reshape2)
start <- Sys.time()
cat("Starting at, ", format(start, "%H:%M:%OS3"), "\n" )

setwd("/labs/csbig/subtipos_mama")

cat("Loading data\n")
load(file="genes.annot.RData")
load(file="adjmatrix_basal.RData")

all.genes.annot <- genes.annot
chrs <- c(as.character(1:22), "X", "Y")

## Extract inter-chromosmal interactions
for (ch1 in chrs) {
  other.chrs <- chrs[which(chrs != ch1)]
  for (ch2 in other.chrs) {
    cat("Working with chromosome ", ch1, " and chromosome ", ch2, "\n")
    genes.ch1 <- all.genes.annot[all.genes.annot$Chr == ch1, ]
    genes.ch2 <- all.genes.annot[all.genes.annot$Chr == ch2, ]
    m1 <- adjmatrix[rownames(adjmatrix) %in% genes.ch1$symbol, 
                    colnames(adjmatrix) %in% genes.ch2$symbol]
    m2 <-  adjmatrix[rownames(adjmatrix) %in% genes.ch2$symbol, 
                     colnames(adjmatrix) %in% genes.ch1$symbol]
    m1[is.na(m1)] <- 0
    m2[is.na(m2)] <- 0
    m3 <- m1 + t(m2)
    cat("Saving plot.\n")
    png(paste("Basal_ch", ch1, "_ch", ch2, ".png", sep=""), width =800 , height = 400)
    myplot <- ggplot(data=melt(m3, value.name = "MI"), aes(MI))+geom_density()
    print(myplot)
    dev.off() 
    cat("Saving data.\n")
    save(m3, file=paste("basal_ch", ch1, "_ch", ch2, ".RData", sep=""))
    write.table(m3, file=paste("basal_ch", ch1, "_ch", ch2, ".tsv", sep=""), 
                quote = F, row.names = T, col.names = T, sep = "\t")
    
  }
}
end <- Sys.time()
cat("Whole thing took ", (end - start), ".\n")
cat("Ending at, ", format(end, "%H:%M:%OS3"), "\n" )

