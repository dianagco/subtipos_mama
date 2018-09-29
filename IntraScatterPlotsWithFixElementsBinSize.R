library(ggplot2)
library(car)
library(reshape2)
start <- Sys.time()
cat("Starting at, ", format(start, "%H:%M:%OS3"), "\n" )

setwd("~/Workspace/subtipos_mama/")

cat("Loading data\n")
conds <- c("basal", "luma", "lumb", "her2", "healthy")

generateDistanciaMIByChromosome <- function () {
  load(file="rdata/lista_distancias_intracromosomales.RData")
  for (cond in conds) {
    cat("Working with condition ", cond, "\n")
    load(file=paste("rdata/adjmatrix_", cond, ".RData", sep=""))
    chrs <- c(as.character(1:22), "X")
    
    parallel::mclapply(X = chrs, mc.cores = 6,  mc.cleanup = FALSE, FUN = function(ch){
      cat("\t Working with chromosome ", ch, "\n")
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
      dist.mi.df$cond <- cond
      dist.mi.df$ch <- ch
      
      write.table(dist.mi.df, file = paste("rdata/pares-distancia-mi/", cond,"_distance_mi", ch, ".tsv", sep=""), sep="\t", 
                             col.names = T, row.names = F, quote = F)
    })
  }
}

generateDistanciaMIByChromosome()

testDifferentBinSizes <- function() {
  bs <- c(seq(from = 10, to = 200, by = 10))
  #bs <- seq(from = 10, to = 20, by = 10)
  chrs <- c(as.character(1:22), "X")
  for (cond in conds) {
    
    cat("Working with condition ", cond, "\n")
    load(file=paste("rdata/adjmatrix_", cond, ".RData", sep=""))
    chrs <- c(as.character(1:22), "X")
    
    c <- lapply(X = chrs, FUN = function(ch){
      dist.mi.df <- read.delim(file=paste("rdata/pares-distancia-mi/", cond,"_distance_mi", ch, ".tsv", sep=""), header = T)
      
      testy <- parallel::mclapply(X = bs, mc.cores = 6,  mc.cleanup = FALSE, FUN = function(b){
        dist.mi.df$bin <- ((as.numeric(rownames(dist.mi.df)) - 1)%/%b) + 1
        dfmeans <- aggregate(mi ~ bin, data=dist.mi.df, FUN=mean, na.rm=TRUE)
        recta <- lm(log10(mi) ~ 1 + log10(bin) + I(log10(bin)^2) + I(log10(bin)^3), data = dfmeans) 
        srecta <- summary(recta)
        return(structure(srecta$adj.r.squared, names = b))
      })
      df <- plyr::ldply(testy)
      rownames(df) <- bs
      colnames(df) <- "r2"
      return(df)
    })
    df <-do.call(cbind, c)
    colnames(df) <- chrs
    write.table(df, file = paste("rdata/",  cond, "_binsize_adjr2_thirdorder. tsv", sep=""), sep="\t", 
                col.names = T, row.names = F, quote = F)
    
    melted <- melt(df)
    colnames(melted) <- c("chr", "adjr2" )
    melted$binsize <- rep(bs, times = 23)
    
    png(paste("plots/", cond, "_binsize_adjr2_thirdorder_200.png", sep=""), width =1200 , height = 800)
    myplot <- ggplot(data=melted, aes(x=binsize, y=adjr2, group=chr)) +
      geom_line(aes(color=chr))+
      geom_point(aes(color=chr))
    print(myplot)
    dev.off()
    
  }
}

#testDifferentBinSizes()

testSecondAndThirdOrderPol <- function(binsize) {
  conds <- c("basal", "healthy")
  chrs <- c(as.character(1:22), "X")
  
  for (cond in conds) {
   
    cat("Working with condition ", cond, "\n")
   
    c <- lapply(X = chrs, FUN = function(ch){
      dist.mi.df <- read.delim(file=paste("rdata/pares-distancia-mi/", cond,"_distance_mi", ch, ".tsv", sep=""), header = T)
      dist.mi.df <- dist.mi.df[with(dist.mi.df, order(distance, -mi)), ]
      
      dist.mi.df$bin <- ((as.numeric(rownames(dist.mi.df)) - 1)%/%binsize) + 1
      dfmeans <- aggregate(cbind(distance, mi)~bin, data=dist.mi.df, FUN=mean, na.rm=TRUE)
      
      recta3 <- lm(log10(mi) ~ 1 + log10(distance) + I(log10(distance)^2) + I(log10(distance)^3), data = dfmeans)
      
      png(paste("plots/intra-fixed-bins-size/", cond, "_ch", ch, "_secondorder.png", sep=""), width =1200 , height = 800)
      plot(mi~bin, data = dfmeans)
      y <- coefficients(recta2)["(Intercept)"] + coefficients(recta2)["log10(bin)"]*log10(dfmeans$bin) + coefficients(recta2)["I(log10(bin)^2)"]*(log10(dfmeans$bin)^2)
      y <- 10^y
      lines(
        x = dfmeans$bin,
        y = y,
        col = "blue"
      )
      dev.off()
      
      recta3 <- lm(log10(mi) ~ 1 + log10(bin) + I(log10(bin)^2) + I(log10(bin)^3), data = dfmeans)
      png(paste("plots/intra-fixed-bins-size/", cond, "_ch", ch, "_thirdorder.png", sep=""), width =1200 , height = 800)
      plot(mi~bin, data = dfmeans)
      y <- coefficients(recta3)["(Intercept)"] + coefficients(recta3)["log10(bin)"]*log10(dfmeans$bin) + coefficients(recta3)["I(log10(bin)^2)"]*(log10(dfmeans$bin)^2) + coefficients(recta3)["I(log10(bin)^3)"]*(log10(dfmeans$bin)^3)
      y <- 10^y
      lines(
        x = dfmeans$bin,
        y = y,
        col = "green"
      )
      dev.off()
     
      t <- anova(recta2, recta3)
      sink(paste("plots/intra-fixed-bins-size/", cond, "_ch", ch, "_secondvsthirdanova.txt", sep=""))
      print(t)
      sink()
      
    })
    closeAllConnections()
    
  }
}
 
#testSecondAndThirdOrderPol(100)

getThirdOrderModel <- function(binsize) {
  conds <- c("basal", "luma", "lumb", "her2", "healthy")
  chrs <- c(as.character(1:22), "X")
  #chrs <- c(as.character(1:2))
  
  for (cond in conds) {
    
    cat("Working with condition ", cond, "\n")
    
    models <- lapply(X = chrs, FUN = function(ch){
      dist.mi.df <- read.delim(file=paste("rdata/pares-distancia-mi/", cond,"/", cond, "_distance_mi", ch, ".tsv", sep=""), header = T)
      dist.mi.df <- dist.mi.df[with(dist.mi.df, order(distance, -mi)), ]
      
      dist.mi.df$bin <- ((as.numeric(rownames(dist.mi.df)) - 1)%/%binsize) + 1
      dfmeans <- aggregate(cbind(distance, mi)~bin, data=dist.mi.df, FUN=mean, na.rm=TRUE)
      
      recta3 <- lm(log10(mi) ~ 1 + log10(distance) + I(log10(distance)^2) + I(log10(distance)^3), data = dfmeans)
      
      png(paste("plots/intra-fixed-bins-size/", cond, "/ch_", ch, ".png", sep=""), width =1200 , height = 800)
      plot(mi~distance, data = dfmeans)
      y <- coefficients(recta3)["(Intercept)"] + coefficients(recta3)["log10(distance)"]*log10(dfmeans$distance) + coefficients(recta3)["I(log10(distance)^2)"]*(log10(dfmeans$distance)^2) + coefficients(recta3)["I(log10(distance)^3)"]*(log10(dfmeans$distance)^3)
      y <- 10^y
      lines(
        x = dfmeans$distance,
        y = y,
        col = "green"
      )
      dev.off()
      
      return(recta3)
      
    })
    
    names(models) <- chrs
    save(models, file=paste("rdata/intra-fixed-bins-size/", cond, "_models.RData", sep=""), compress="xz")
    
  }
}

#getThirdOrderModel(100)
