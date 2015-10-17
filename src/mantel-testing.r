#!/usr/bin/Rscript

library(Hmisc)
library(pander)
library(plyr)
library(reshape2)
library(vegan)

source('mantel-testing-functions.R')

#' Create distance matrices

unwanted <- c('week', 'totalNumberSwineAccessions', 'Unk')
ind <- which(!colnames(real.case.data) %in% unwanted)
observed <- real.case.data[, ind]

CheckCrossCorLagSensitivity <- function(observed){
  lags <- 0:10
  tmpf <- function(x, y){
    sx <- observed[, x]
    sy <- observed[, y]
    ret <- ccf(sx, sy, plot=FALSE, demean=TRUE)
    ret <- ret[lags, ]
    ret$acf[, ,1]
  }
  grd <- expand.grid(x=1:ncol(observed), y=1:ncol(observed))
  grd <- grd[grd$x != grd$y, ]
  cross.cors <- mapply(x=grd$x, y=grd$y, tmpf)

  png('cross-cors-by-lag-boxplot.png')
  boxplot(t(cross.cors), names=lags, xlab='lag', ylab='cross-correlation')
  dev.off()

  lag.max.ind <- apply(cross.cors, 2, which.max)
  lag.max <- lags[lag.max.ind]
  png('counts-of-times-each-lag-was-max-for-a-pair.png')
  plot(table(lag.max))
  dev.off()
}
CheckCrossCorLagSensitivity(observed)

pop.struct.mats <- MakePopStructMats(observed)
pop.dyn.mats <- MakePopDynMats(observed, lags.sel=0:2)

#' Test for association between distance matrices

struct.v.struct.tests <- DoMantelTests(pop.struct.mats, pop.struct.mats,
                                       permutations=1e4)
dyn.v.struct.tests <- DoMantelTests(pop.dyn.mats, pop.struct.mats,
                                    permutation=1e4)

sink('mantel-table-population-structure-matrices.txt')
pander(struct.v.struct.tests)
sink()
sink('mantel-table-population-structure-vs-dynamics-matrices.txt')
pander(dyn.v.struct.tests)
sink()

#' Get descriptive stats

getDesc <- function(M) {
  M[upper.tri(M)]
}
desc <- sapply(c(pop.dyn.mats['lag1'], pop.struct.mats), getDesc)
invisible(latex(describe(desc), file='mantel-describe.tex'))

#' Run partial tests

partial.tests <- DoPartialMantelTests(pop.dyn.mats, pop.struct.mats['shipment'],
                               pop.struct.mats['gcd'], permutation=1e3)

sink('mantel-partial-table.txt')
pander(partial.tests)
sink()

save.image(file='mantel-testing-checkpoint1.RData')
