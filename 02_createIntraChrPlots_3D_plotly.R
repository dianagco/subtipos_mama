library(plotly)

Sys.setenv("plotly_username"="dianagco")
Sys.setenv("plotly_api_key"="jh47zwssSzbvn6gUsfAE")

start <- Sys.time()
cat("Starting at, ", format(start, "%H:%M:%OS3"), "\n" )

setwd("~/Workspace/subtipos_mama")

#conditions <- c("healthy", "luma", "lumb", "her2", "basal")
conditions <- c("basal")
chrs <- c(as.character(1:22), "X")
cond <- conditions[1]
ch <- chrs[1]

#for (ch in chrs) {
  mi.df <- read.delim(paste(cond, "_intra_order_", ch, ".tsv", sep=""),
                        stringsAsFactors = T, header = T, sep="\t" )
  mi.df <- mi.df[c(1,3,2)]
  
  mi.agg <- aggregate(mi ~ order, data=mi.df, FUN = cbind)
  
  max_length <- max(unlist(lapply (mi.agg$mi, FUN = length)))
  
  mi.matrix <- sapply (mi.agg$mi, function (x) {
                  length (x) <- max_length; 
                  return (x)
                })
  mi.matrix[is.na(mi.matrix)] <- 0 
  min.matrix <- mi.matrix[1:50, 1:50]
  p <- plot_ly(z = ~min.matrix,  type = 'surface')
  p
#}
