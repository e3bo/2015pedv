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

GetCrossCorrs <- function(){
    n <- ncol(observed)
    CC <- matrix(nrow=n, ncol=n)
    getcc <- function(x,y, lag=1){
        foo <- ccf(x, y, plot=FALSE)
        ind <- which(foo$lag == lag)
        foo$acf[ind]
    }
    for(i in seq_len(n)){
        for(j in seq_len(n)){
            ## CC[i,j] will be high if deviations from the mean in series i
            ## are shifted to the left of similar deviations to the mean in series j
            ## i.e. i's deviations are indicative of j's future deviations
            CC[i,j] <- getcc(observed[,j], observed[,i])
        }
    }
    colnames(CC) <- rownames(CC) <- colnames(observed)
    CC
}

CC <- GetCrossCorrs()
CC <- CC[nms, nms]

check_cc_directionality <- function(){
    t <- 1:100 * .25
    x <- sin(t)
    ## y is shifted to the right
    y <- sin(t-1)

    xylag1 <- getcc(x, y)
    yxlag1 <- getcc(y,x)
    main <- paste('getcc(x, y) ==', round(xylag1,2),
                  '; getcc(y,x) == ', round(yxlag1,2))
    if(xylag1 > yxlag1){
        conc <- 'left arg delayed by lag'
    } else{
        conc <- 'right arg delayed by lag'
    }
    plot(x~t, type='l', ylab='f(t)', main=main, sub=conc)
    lines(y~t, col=2)
    legend('topright', col=1:2, legend=c('x', 'y'), lty=1)
}

png('direction-check.png')
check_cc_directionality()
dev.off()

## Hypothesis tesing

getMat <- function(x) switch(x, 'shipment'=epl, 'cor'=CC,
                             'gcd'=-centerDists, 'sharedBord'=nhood)

doTest <- function(M1, M2, symmetrize=FALSE, ...){
    x <- getMat(M1)
    y <- getMat(M2)
    if(symmetrize){
        x <- (x + t(x)) * 0.5
        y <- (y + t(y)) * 0.5
    }
    mantel(x, y, ...)
}

mats <- c('shipment', 'cor', 'gcd', 'sharedBord')
methods <- c('spearman', 'pearson')
symmetrize <- c(TRUE, FALSE)
des <- expand.grid(M1=mats, M2=mats, method=methods, symmetrize=symmetrize,
                   stringsAsFactors=FALSE)
des <- des[des$M1 != des$M2, ]
des$permutations <- 10000

res <- list()
for(i in seq_len(nrow(des))){
    print(i)
    res[[i]] <- do.call(doTest, as.list(des[i,]))
}

des$r <- sapply(res, '[[', 'statistic')
des$pValues <- sapply(res, '[[', 'signif')

sink('mantel-table.txt')
pander(des)
sink()

## get descriptive stats

getDesc <- function(x) {
  M <- getMat(x)
  M[upper.tri(M)]
}
desc <- sapply(c('shipment', 'cor', 'gcd'), getDesc)
invisible(latex(describe(desc), file='mantel-describe.tex'))

## partial tests

doPartialTest <- function(M1, M2, M3, symmetrize=FALSE, ...){
    x <- getMat(M1)
    y <- getMat(M2)
    z <- getMat(M3)
    if(symmetrize){
        x <- (x + t(x)) * 0.5
        y <- (y + t(y)) * 0.5
    }
    mantel.partial(x, y, z, ...)
}

desp <- data.frame(M1='shipment', M2='cor', M3='gcd',
                   symmetrize=c(TRUE, FALSE),
                   permutations=1000, stringsAsFactors=FALSE)

resp <- list()
for(i in seq_len(nrow(desp))){
    print(i)
    resp[[i]] <- do.call(doPartialTest, as.list(desp[i, ]))
}

desp$r <- sapply(resp, '[[', 'statistic')
desp$pValues <- sapply(resp, '[[', 'signif')

sink('mantel-partial-table.txt')
pander(desp)
sink()

save.image(file='mantel-testing-checkpoint1.RData')
