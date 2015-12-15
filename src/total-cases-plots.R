#!/usr/bin/Rscript

load('sim-study-checkpoint3.rda')

tmpf <- function(x) sum(x$mantel.observed)
total.cases <- sapply(res, tmpf)

pdf('total-cases-plots.pdf')
plot(df$tprob.net, total.cases)
plot(df$tprob.sp, total.cases)
plot(df$raster.ncol, total.cases)
plot(df$seasonal.amplitude, total.cases)
plot(df$tprob.outside, total.cases)
dev.off()
