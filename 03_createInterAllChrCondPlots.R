library(ggplot2)
library(RColorBrewer)

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
start <- Sys.time()
cat("Starting at, ", format(start, "%H:%M:%OS3"), "\n" )

setwd("/labs/csbig/subtipos_mama")
#setwd("~/Workspace/subtipos_mama")

conditions <- c("healthy", "luma", "lumb", "her2", "basal")

g <- parallel::mclapply(X = conditions, mc.cores = 5,  mc.cleanup = FALSE, FUN = function(con) {
  cat("Working with condition ", con, "\n")
  cat("Loading data\n")
  
  load(file=paste("adjmatrix_", con, ".RData", sep="")) 
  all.mis <- data.frame(mi = unlist(adjmatrix[!is.na(adjmatrix) & adjmatrix != 1]))
  all.mis$condition <- con
  colnames(all.mis) <- c("mi", "condition")
  png(paste("All_", con, ".png", sep=""), width = 1200, height = 600)
  myplot <- ggplot(data=all.mis) + 
    geom_density(aes(x = mi, colour = myPalette(1)), show.legend = F) +
    coord_cartesian(xlim = c(0, 0.25)) 
  print(myplot)
  dev.off() 
})
end <- Sys.time()
cat("Whole thing took ", (end - start), ".\n")
cat("Ending at, ", format(end, "%H:%M:%OS3"), "\n" )

    
