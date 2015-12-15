#!/usr/bin/Rscript

load('sim-study-checkpoint3.rda')

MakeDynPlot <- function(i){
  o <- res[[i]]$mantel.observed
  par(mfrow=c(2, 1))
  p <- unlist(df[i, c('seasonal.amplitude', 'tprob.sp', 'tprob.net', 'starting.grid.x')])
  dotchart(p)
  plot.new()
  text(x=0.5, y=0.5, labels=paste('raster columns:', df[i, 'raster.ncol']))
  heatmap(log(o + 1), Rowv=NA)
}
lapply(1:10, MakeDynPlot)
dev.off()
