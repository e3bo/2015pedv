#+setup, include=FALSE, cache=FALSE
library(knitr)
opts_chunk$set(fig.path='simulator-figures/', fig.align='center', fig.show='hold')
options(replace.assign=TRUE,width=80)
Sys.setlocale("LC_TIME", "C") #Needed for identical()
Sys.setlocale("LC_COLLATE", "C")

#' ## Load Packages

set.seed(4253)
library(maptools)
library(MASS)
library(raster)
library(reshape2)
library(rgdal)
library(igraph)
library(vegan)
print(sessionInfo())

#' ## Source R scripts

source('mantel-testing-functions.R')
source('regression-modeling-functions.R')

#' ## Load data

ea.proj <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs" ## NAD83 Lambert Azimuthal Equal Area

## cb.proj <- "+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs" ## from .proj file
## cb <- readShapeSpatial("cb_2014_us_county_500k/cb_2014_us_county_500k")
##                         proj4string=CRS(cb.proj))
cb <- readOGR('cb_2014_us_county_500k', 'cb_2014_us_county_500k')

#' ## Sample random coordinates

cb.ea <- spTransform(cb, CRS(ea.proj))

tmpf <- function(geo, countyData){
  sf <- countyData$STFIPS
  cf <- countyData$COFIPS
  sf <- formatC(sf, flag="0", format="d", width=2)
  cf <- formatC(cf, flag="0", format="d", width=3)
  test1 <- geo@data[, 'STATEFP'] %in% sf
  test2 <- geo@data[, 'COUNTYFP'] %in% cf
  test <- test1 & test2
  geo[test, ]
}
cb2 <- tmpf(geo=cb.ea, countyData=county.hogs.pigs)

plotAllCountiesWithFarmData <- function(){
    plot(cb2)
}

getSamp <- function(spdf=cb2, cfps, sfps, n, plot.samp=FALSE){
    if(n == 0) {
        samp <- NULL
    } else {
        cfps <- formatC(cfps, flag="0", format="d", width=3)
        sfps <- formatC(sfps, flag="0", format="d", width=2)
        ind <- which(spdf@data[, 'COUNTYFP']==cfps & spdf@data[, 'STATEFP']==sfps)
        cty.poly <- spdf[ind, ]
        samp <- spsample(cty.poly, type='random', n=n, iter=10)
        if(plot.samp){
            plot(cty.poly)
            points(samp)
        }
    }
    samp
}

testSamplingOfMultiPolygonCounties <- function(){
   getSamp(cfps=21, sfps=25, n=100, plot.samp=TRUE)
}

testSamplingOfCountiesWithHoles <- function(){
   getSamp(cfps=15, sfps=51, n=100, plot.samp=TRUE)
}

coord.samps <- mapply(getSamp, cfps=county.hogs.pigs$COFIPS,
                      sfps=county.hogs.pigs$STFIPS,
                      n=county.hogs.pigs$DATA, SIMPLIFY=FALSE)

#' #Make raster grid

r <- raster(cb2)
res(r) <- 1600 * 10
#adj <- adjacent(r, 1:ncell(r), directions=8) ## this could be quickly calculated as needed

tmpf <- function(xy) {
  if(is.null(xy)) {
    NULL
  } else {
    cellFromXY(object=r, xy=xy)
  }
}
cell.samps <- lapply(coord.samps, tmpf)

plotCellCounts <- function(){
  tab <- table(unlist(cell.samps))
  r[as.integer(names(tab))] <- tab
  plot(r)
}

## convert fips to abbreviation
stfips.rle <- list(lengths=county.hogs.pigs$DATA, values=county.hogs.pigs$STFIPS)
data(state.fips, package='maps')
key <- match(county.hogs.pigs$STFIPS, state.fips$fips)
county.hogs.pigs$abb <- as.character(state.fips$abb[key])

## create agent data structures

adf <- data.frame(cell=unlist(cell.samps),
                  abb=rep(county.hogs.pigs$abb,
                      times=county.hogs.pigs$DATA),
                  infection.time=NA,
                  recovery.time=NA)

sp.nbs <- vector('list', length=nrow(adf))

occupied.cells <- unique(adf$cell)
cell2id <- lapply(occupied.cells, function(x) which(adf$cell == x))
names(cell2id) <- occupied.cells

#' calc weight for sbm

tmpf <- function(rel='directed', flowMat, flows, fiAll) {
    fmo <- t(flowMat[state.abb, state.abb])
    diag(fmo) <- flows[state.abb, 'impInternalFlow']
    Fund <- fmo + t(fmo)
    Fdir <- fmo
   ## Fdir[i, j] == head sent to i from j
    if(rel == 'directed') {
        F <- fmo
    } else {
        F <- fmo + t(fmo)
    }
    Fsum <- rowSums(F)
    Fnorm <- F/Fsum
    n <- fiAll[state.abb, ]$X2002.Farms
    W <- Fnorm * n
    W <- t(t(W) / n)
    list(weights=W, Fsum=Fsum)
}
wcDir <- tmpf('directed', flowMat=flows.matrix, flows=internal.flows,
              fiAll=farms.and.inventory)

adf <- adf[order(adf$abb), ]
state.tots <- rle(as.character(adf$abb))
wcDir$weights <- wcDir$weights[state.tots$values, state.tots$values]

g <- sample_sbm(sum(state.tots$lengths), pref.matrix=wcDir$weights / 1000,
                block.sizes=state.tots$lengths, directed=TRUE)

## initialize infections

                                        #cases <- which(adf$cell %in% c(27799))
adf$infection.time <- NA
adf$recovery.time <- NA
cases <- sample.int(nrow(adf), 100)

get.nbs <- function(cell, rast=r) {
  nb <- adjacent(x=rast, cell, directions=8)
  nb.cells <- nb[, 'to']
  cell.names <- as.character(c(cell, nb.cells))
  unlist(cell2id[cell.names])
}
sp.nbs[cases] <- lapply(adf$cell[cases], get.nbs)
adf$infection.time[cases] <- 0

nsteps <- 38
step <- 1
tprob <- 0.001
tprob.net <- 0.8
rprob <- 0.5
seasonal.factor <- function(x) -sinpi((x + 3)/ 52 * 2)
seasonal.amplitude <- 0

while(step < nsteps){
  new.cases <- repeat.cases <- integer(0)
  for (case in cases){
    ## spatial transmission
    contacts <- sp.nbs[[case]]
    is.susceptible <- is.na(adf$infection.time[contacts])
    contacts <- contacts[is.susceptible]
    rand <- runif(n=length(contacts))
    test <- rand < tprob * (1 + seasonal.amplitude * seasonal.factor(step))
    if(any(test)){
      new.cases <- c(new.cases, contacts[test])
    }
    ## network transmission
    contacts <- neighbors(g, v=case, mode='out')
    is.susceptible <- is.na(adf$infection.time[contacts])
    contacts <- contacts[is.susceptible]
    rand <- runif(n=length(contacts))
    test <- rand < tprob.net * (1 + seasonal.factor(step))
    if(any(test)){
      new.cases <- c(new.cases, contacts[test])
    }
    rand2 <- runif(n=1)
    if(rand2 > rprob){
      repeat.cases <- c(repeat.cases, case)
    } else {
      adf$recovery.time[case] <- step
    }
  }
  print(system.time(sp.nbs[new.cases] <- lapply(adf$cell[new.cases], get.nbs)))
  adf$infection.time[new.cases] <- step
  cases <- c(new.cases, repeat.cases)
  cat('step: ', step, '\n')
  step <- step + 1
}

## Tabulation of case counts

step <- seq(1, to=nsteps)

events.by.state <- function(x, what) {
  tapply(adf[[what]] < x, adf$abb, sum, na.rm=TRUE)
}
cum.infections <- sapply(step, events.by.state, what='infection.time')
cum.recoveries <- sapply(step, events.by.state, what='recovery.time')
no.infected <- cum.infections - cum.recoveries
new.cases <- t(apply(cum.infections, 1, diff))
new.cases <- cbind(cum.infections[, 1], new.cases)

tmpf <- function(x, size=1.75, prep=.51) {
  xrep <- rbinom(x, size=x, prob=prep)
  rnbinom(n=xrep, size=size, mu=.76 + 1.92 *xrep)
}
reports <- t(apply(new.cases, 1, tmpf))

observed <- t(reports)
observed <- observed[, colSums(observed) > 0]

## Mantel tests

pop.struct.mats <- MakePopStructMats(observed)
pop.dyn.mats <- MakePopDynMats(observed)
mantel.tests <- DoMantelTests(pop.dyn.mats, pop.struct.mats, permutations=1e2)

## Regression model

om <- GetRegressionData(case.data=as.data.frame(observed), flows=internal.flows,
                        fiAll=farms.and.inventory, stateCty=state.summaries)

m <- list()
f <- list()
f$lm <- as.formula(clog(cases) ~ clogCases1wa + logInternalFlowScaled + logCmedDenseScaled + weekCent + offset(log(nSusc1WksAgo) - 2*log(nFarms)))
m$lm <- lm(f$lm, data=om)
summary(m$lm)
