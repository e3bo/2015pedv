#!/usr/bin/Rscript

library(Hmisc)
library(igraph)
library(mapproj)
library(maps)
library(pander)
library(plyr)
library(reshape2)
library(vegan)

Sys.setlocale("LC_TIME", "C") #Needed for identical()
Sys.setlocale("LC_COLLATE", "C")

## Get case data

dataDir <- '.'

tmpf <- function(){
  fn <- file.path(dataDir, 'PEDvweeklyreport-state-ts-01-08-14.csv')
  ret <- read.csv(fn)
  ret[1:9, c('CA', 'MD', 'NE', 'WY')] <- 0
  key <- order(colnames(ret)[-c(1:2)])
  ret <- cbind(ret[, c(1:2)], ret[, -c(1:2)][, key])
  ret$Unk <- NULL
  target <- structure(list(week = structure(c(20L, 29L, 32L, 8L, 12L),
.Label = c("10/13/2013", "10/20/2013", "10/27/2013", "10/6/2013",
"11/10/2013", "11/17/2013", "11/24/2013", "11/3/2013", "12/1/2013",
"12/15/2013", "12/22/2013", "12/29/2013", "12/8/2013", "4/15/2013",
"4/22/2013", "4/29/2013", "5/13/2013", "5/20/2013", "5/27/2013",
"5/6/2013", "6/10/2013", "6/16/2013", "6/23/2013", "6/3/2013",
"6/30/2013", "7/14/2013", "7/21/2013", "7/28/2013", "7/7/2013",
"8/11/2013", "8/18/2013", "8/25/2013", "8/4/2013", "9/1/2013",
"9/15/2013", "9/22/2013", "9/29/2013", "9/8/2013"), class = "factor"),
totalNumberSwineAccessions = c(17L, 34L, 26L, 90L, 134L), CA = c(0, 0,
0, 0, 1), CO = c(1L, 1L, 0L, 0L, 0L), IA = c(8L, 6L, 2L, 38L, 54L), IL
= c(0L, 1L, 0L, 1L, 14L), IN = c(3L, 2L, 1L, 1L, 3L), KS = c(0L, 4L,
3L, 6L, 4L), KY = c(0L, 0L, 0L, 0L, 0L), MD = c(0, 0, 0, 0, 0), MI =
c(0L, 0L, 0L, 2L, 0L), MN = c(1L, 2L, 2L, 7L, 20L), MO = c(0L, 0L, 0L,
2L, 4L), NC = c(0L, 3L, 4L, 14L, 18L), NE = c(0, 0, 0, 0, 2), NY =
c(0L, 0L, 0L, 0L, 0L), OH = c(0L, 2L, 1L, 5L, 5L), OK = c(0L, 11L,
10L, 10L, 2L), PA = c(1L, 0L, 3L, 1L, 0L), SD = c(0L, 0L, 0L, 0L, 3L),
TN = c(0L, 0L, 0L, 1L, 1L), TX = c(0L, 0L, 0L, 1L, 1L), WI = c(0L, 0L,
0L, 1L, 0L ), WY = c(0, 0, 0, 0, 1)), .Names = c("week",
"totalNumberSwineAccessions", "CA", "CO", "IA", "IL", "IN", "KS",
"KY", "MD", "MI", "MN", "MO", "NC", "NE", "NY", "OH", "OK", "PA",
"SD", "TN", "TX", "WI", "WY" ), row.names = c(4L, 13L, 20L, 30L, 38L),
class = "data.frame")
  stopifnot(identical(target, ret[c(4,13,20,30,38),]))
  ret
}
caseData <- tmpf()

unwanted <- c('week', 'totalNumberSwineAccessions', 'Unk')
ind <- which(!colnames(caseData) %in% unwanted)
observed <- caseData[, ind]

## Get cross correlations

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

## Get shared-border neighborhoods

nbEdgelist <- read.csv(file.path(dataDir, 'state_neighbors_fips.txt'), header=FALSE)
data(state.fips)
g <- graph.data.frame(nbEdgelist, directed=FALSE)
g <- simplify(g)
key <- match(V(g)$name, state.fips$fips)
abb <- state.fips$abb[key]
V(g)$name <- as.character(abb)
vids <- which(V(g)$name %in% colnames(caseData))
g2 <- induced.subgraph(g, vids)
nms <- colnames(observed)
nhood <- as.matrix(g2[nms,nms])

## Get shipment flows

epO <- ep <- read.csv(file.path(dataDir, 'shipment-flows-origins-on-rows-dests-on-columns.csv'),row.names=1)
key <- match(colnames(observed), colnames(ep))
ep <- ep[key, key]
rownames(ep) <- colnames(ep)

### assertions based on inspection of original .xls file
stopifnot(ep['CA', 'IL'] == 9415)
stopifnot(ep['MI', 'KS'] == 6)
stopifnot(ep['IL', 'IA'] == 1424813)
stopifnot(ep['KY', 'IA'] == 17658)
stopifnot(ep['MO', 'IA'] == 2389932)
stopifnot(ep['OK', 'CA'] == 16762)
stopifnot(ep['MO', 'CA'] == 645)

epl <- log10(ep +1)
epl <- data.matrix(epl)
## We discard the within-state flows because they do not enter into
## the analysis. Thus it is mostly a matter of how the plots will
## look, and because the original spreadsheet had no data on
## within-state flows (they were added for some other analyses based
## on an assumption of 90% of total flow), there's no good reason to
## include them in the plot.
diag(epl) <- NA

key <- match(colnames(ep), colnames(CC))
CC <- CC[key,key]

## Get great circle distance

key <- match(colnames(CC), state.abb)

cx <- state.center$x[key]
cy <- state.center$y[key]

n <- ncol(CC)
centerDists <- matrix(nrow=n, ncol=n)

# Calculates the geodesic distance between two points specified by radian latitude/longitude using the
# Haversine formula (hf)
# source: http://www.r-bloggers.com/great-circle-distance-calculations-in-r/
gcd.hf <- function(long1, lat1, long2, lat2) {
      R <- 6371 # Earth mean radius [km]
        delta.long <- (long2 - long1)
        delta.lat <- (lat2 - lat1)
        a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
        c <- 2 * asin(min(1,sqrt(a)))
        d = R * c
        return(d) # Distance in km
  }

deg2rad <- function(deg) return(deg*pi/180)

getCenterDist <- function(state1,state2){
    long1 <- deg2rad(cx[state1])
    long2 <- deg2rad(cx[state2])
    lat1 <- deg2rad(cy[state1])
    lat2 <- deg2rad(cy[state2])
    gcd.hf(long1, lat1, long2, lat2)
}

for(i in seq_len(n)){
    for(j in seq_len(n)){
        centerDists[i,j] <- getCenterDist(i, j)
    }
}

colnames(centerDists) <- rownames(centerDists) <- colnames(CC)

getDesc <- function(x) {
  M <- getMat(x)
  M[upper.tri(M)]
}

## hypothesis tesing

getMat <- function(x) switch(x, 'shipment'=epl, 'cor'=CC,
                             'gcd'=-centerDists, 'sharedBord'=nhood)

doTest <- function(M1, M2, symmetrize=FALSE, ...){
    x <- getMat(M1)
    y <- getMat(M2)
    if(symmetrize){
        x <- (x + t(x))*0.5
        y <- (y + t(y))*0.5
    }
    mantel(x, y, ...)
}

mats <- c('shipment', 'cor', 'gcd', 'sharedBord')
methods <- c('spearman', 'pearson')
symmetrize <- c(TRUE, FALSE)
des <- expand.grid(M1=mats, M2=mats, method=methods, symmetrize=symmetrize, stringsAsFactors=FALSE)
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

## get Descriptive stats

desc <- sapply(c('shipment', 'cor', 'gcd'), getDesc)
invisible(latex(describe(desc), file='mantel-describe.tex'))

## partial tests


doPartialTest <- function(M1, M2, M3, symmetrize=FALSE, ...){
    x <- getMat(M1)
    y <- getMat(M2)
    z <- getMat(M3)
    if(symmetrize){
        x <- (x + t(x))*0.5
        y <- (y + t(y))*0.5
    }
    mantel.partial(x, y, z, ...)
}

desp <- data.frame(M1='shipment', M2='cor', M3='gcd', symmetrize=c(TRUE, FALSE),
                   permutations=1000, stringsAsFactors=FALSE)

resp <- list()
for(i in seq_len(nrow(desp))){
    print(i)
    resp[[i]] <- do.call(doPartialTest, as.list(desp[i,]))
}

desp$r <- sapply(resp, '[[', 'statistic')
desp$pValues <- sapply(resp, '[[', 'signif')

sink('mantel-partial-table.txt')
pander(desp)
sink()

save.image(file='mantel-testing-checkpoint1.RData')
