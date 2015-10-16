#!/usr/bin/Rscript

library(Hmisc)
library(pander)
library(plyr)
library(reshape2)
library(vegan)

load('common-data.RData')

unwanted <- c('week', 'totalNumberSwineAccessions', 'Unk')
ind <- which(!colnames(real.case.data) %in% unwanted)
observed <- real.case.data[, ind]

nms <- colnames(observed)
nhood <- shared.border.adjacency[nms, nms]
ep <- flows.matrix[nms, nms]
epl <- log10(ep + 1)
centerDists <- state.to.state.dists[nms, nms]

## Check cross correlations at different lags

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

GetCrossCorrs <- function(lag, obs){
    n <- ncol(obs)
    CC <- matrix(nrow=n, ncol=n)
    getcc <- function(x,y, lag){
        foo <- ccf(x, y, plot=FALSE)
        ind <- which(foo$lag == lag)
        foo$acf[ind]
    }
    for(i in seq_len(n)){
        for(j in seq_len(n)){
            ## CC[i,j] will be high if deviations from the mean in series i
            ## are shifted to the left of similar deviations to the mean in series j
            ## i.e. i's deviations are indicative of j's future deviations
            CC[i, j] <- getcc(obs[, j], obs[, i], lag=lag)
        }
      }

    colnames(CC) <- rownames(CC) <- colnames(obs)
    CC
}
lags.sel <- 0:2
CC <- lapply(lags.sel, GetCrossCorrs, obs=observed)
names(CC) <- paste0('lag', lags.sel)
stopifnot(colnames(CC$lag1) == nms)


## Testing for association between distance matrices based on structure

pop.struct.mats <- list('shipment'=epl, 'gcd'=-centerDists, 'sharedBord'=nhood)
pop.dyn.mats <- CC

doTest <- function(x, y, symmetrize=FALSE, ...){
    if(symmetrize){
        x <- (x + t(x)) * 0.5
        y <- (y + t(y)) * 0.5
    }
    mantel(x, y, ...)
}

methods <- c('spearman', 'pearson')
symmetrize <- c(TRUE, FALSE)
des <- expand.grid(M1=names(pop.struct.mats), M2=names(pop.struct.mats),
                   method=methods, symmetrize=symmetrize,
                   stringsAsFactors=FALSE)
des <- des[des$M1 != des$M2, ]
des$permutations <- 10000

res <- mapply(doTest, x=pop.struct.mats[des$M1], y=pop.struct.mats[des$M2],
              symmetrize=des$symmetrize, permutations=des$permutations,
              SIMPLIFY=FALSE)

des$r <- sapply(res, '[[', 'statistic')
des$pValues <- sapply(res, '[[', 'signif')

sink('mantel-table-population-structure-matrices.txt')
pander(des)
sink()

## Testing for associaton between distance matrices based on dynamics and structure

des.dyn.vs.struct <- expand.grid(M1=names(pop.dyn.mats),
                                 M2=names(pop.struct.mats),
                                 method=methods, symmetrize=symmetrize,
                                 stringsAsFactors=FALSE)
des.dyn.vs.struct$permutations <- 10000

res <- mapply(doTest, x=pop.dyn.mats[des.dyn.vs.struct$M1],
              y=pop.struct.mats[des.dyn.vs.struct$M2],
              symmetrize=des.dyn.vs.struct$symmetrize,
              permutations=des.dyn.vs.struct$permutations,
              SIMPLIFY=FALSE)

des.dyn.vs.struct$r <- sapply(res, '[[', 'statistic')
des.dyn.vs.struct$pValues <- sapply(res, '[[', 'signif')

sink('mantel-table-population-structure-vs-dynamics-matrices.txt')
pander(des.dyn.vs.struct)
sink()

## get descriptive stats

getDesc <- function(M) {
  M[upper.tri(M)]
}
desc <- sapply(c(pop.dyn.mats['lag1'], pop.struct.mats), getDesc)
invisible(latex(describe(desc), file='mantel-describe.tex'))

## partial tests

doPartialTest <- function(x, y, z, symmetrize=FALSE, ...){
    if(symmetrize){
        x <- (x + t(x)) * 0.5
        y <- (y + t(y)) * 0.5
    }
    mantel.partial(x, y, z, ...)
}

desp <- data.frame(M1='shipment', M2='lag0', M3='gcd',
                   symmetrize=c(TRUE, FALSE),
                   permutations=1000, stringsAsFactors=FALSE)

resp <- mapply(doPartialTest, x=pop.struct.mats[desp$M1],
               y=pop.dyn.mats[desp$M2], z=pop.struct.mats[desp$M3],
               symmetrize=desp$symmetrize,
               permutations=desp$permutations,
               SIMPLIFY=FALSE)

desp$r <- sapply(resp, '[[', 'statistic')
desp$pValues <- sapply(resp, '[[', 'signif')

sink('mantel-partial-table.txt')
pander(desp)
sink()

save.image(file='mantel-testing-checkpoint1.RData')
