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
cb2 <- tmpf(geo=cb.ea, countyData=county.hogs.pigs.02)

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

n <- county.hogs.pigs.02$DATA
coord.samps <- mapply(getSamp, cfps=county.hogs.pigs.02$COFIPS,
                      sfps=county.hogs.pigs.02$STFIPS,
                      n=n, SIMPLIFY=FALSE)

#' #Make raster grid

r <- raster(cb2)
res(r) <- 1600 * 10

my.8.nbs <- function(r, cell){
  nc <- ncol(r)
  nr <- nrow(r)
  row.num <- (cell - 1) %/% nc + 1
  col.num <- (cell - 1) %% nc + 1
  nbs <- list()
  i <- 1
  ## First row
  if (row.num > 1){
      nbs[[i]] <- cell - nc
      i <- i + 1
      if (col.num > 1){
          nbs[[i]] <- cell - nc - 1
          i <- i + 1
      }
      if (col.num < nc){
          nbs[[i]] <- cell - nc + 1
          i <- i + 1
      }
  }
  ## Second row
  if (col.num > 1){
      nbs[[i]] <- cell - 1
      i <- i + 1
  }
  nbs[[i]] <- cell
  i <- i + 1
  if (col.num < nc){
      nbs[[i]] <- cell + 1
      i <- i + 1
  }
  ## Last row
  if (row.num < nr){
      nbs[[i]] <- cell + nc
      i <- i + 1
      if (col.num > 1){
          nbs[[i]] <- cell + nc - 1
          i <- i + 1
      }
      if (col.num < nc){
          nbs[[i]] <- cell + nc + 1
          i <- i + 1
      }
  }
  unlist(nbs)
}

tmpf <- function(x) my.8.nbs(r, x)
adj <- lapply(1:ncell(r), tmpf)

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
stfips.rle <- list(lengths=county.hogs.pigs.02$DATA,
                   values=county.hogs.pigs.02$STFIPS)
data(state.fips, package='maps')
key <- match(county.hogs.pigs.02$STFIPS, state.fips$fips)
county.hogs.pigs.02$abb <- as.character(state.fips$abb[key])

## create agent data structures

adf <- data.frame(cell=unlist(cell.samps),
                  abb=rep(county.hogs.pigs.02$abb,
                      times=n),
                  infection.time=NA,
                  recovery.time=NA)

occupied.cells <- unique(adf$cell)
cell2id <- lapply(occupied.cells, function(x) which(adf$cell == x))
names(cell2id) <- occupied.cells

#' calc weight for sbm

tmpf <- function(rel='directed', flowMat, flows, farms.by.state) {
    fmo <- t(flowMat[state.abb, state.abb])
    diag(fmo) <- flows[state.abb, 'impInternalFlow']
    ## fmo[i, j] == head sent to i from j
    if(rel == 'directed') {
        F <- fmo
    } else {
        F <- fmo + t(fmo)
    }
    Fsum <- rowSums(F)
    Fnorm <- F / Fsum
    n <- farms.by.state[state.abb]
    W <- Fnorm * n
    W <- t(t(W) / n)
    list(weights=W, Fsum=Fsum)
}

adf <- adf[order(adf$abb), ]
state.tots <- rle(as.character(adf$abb))
names(state.tots$lengths) <- state.tots$values

wcDir <- tmpf('directed', flowMat=flows.matrix, flows=internal.flows,
              farms.by.state=state.tots$lengths)

wcDir$weights <- wcDir$weights[state.tots$values, state.tots$values]

g <- sample_sbm(sum(state.tots$lengths), pref.matrix=wcDir$weights / 1000,
                block.sizes=state.tots$lengths, directed=TRUE)

get.nbs <- function(cell, nb.adj=adj) {
  #nb <- adjacent(x=rast, cell, directions=8, include=TRUE)
  #nb.cells <- nb[, 'to']
  #ind <- nb.adj[, 'from'] == cell
  nb.cells <- adj[[cell]]
  cell.names <- as.character(nb.cells)
  unlist(cell2id[cell.names])
}
sp.nbs <- lapply(adf$cell, get.nbs)

net.nbs <- adjacent_vertices(g, v=V(g), mode='out')
net.nbs <- sapply(net.nbs, as.integer)

run.sims <- function(adf, verbose=TRUE) {
  while(step < nsteps){
    sf <- (1 + seasonal.amplitude * seasonal.factor(step))
    net.contact.dist <- table(unlist(net.nbs[cases]))
    sp.contact.dist <- table(unlist(sp.nbs[cases]))
    tprob.sp.t <- tprob.sp * sf
    tprob.net.t <- tprob.net * sf

    avoidance.probs.sp <- (1 - tprob.sp.t)^sp.contact.dist
    avoidance.probs.net <- (1 - tprob.net.t)^net.contact.dist

    rand <- runif(n=length(avoidance.probs.sp))
    test <- rand > avoidance.probs.sp
    ids <- as.integer(names(sp.contact.dist))
    new.cases.sp <- ids[test]

    rand <- runif(n=length(avoidance.probs.net))
    test <- rand > avoidance.probs.net
    ids <- as.integer(names(net.contact.dist))
    new.cases.net <- ids[test]

    new.cases <- unique(c(new.cases.net, new.cases.sp))
    is.susceptible <- is.na(adf$infection.time[new.cases])
    new.cases <- new.cases[is.susceptible]
    adf$infection.time[new.cases] <- step

    rand <- runif(n=length(cases))
    test <- rand > rprob
    repeat.cases <- cases[test]
    adf$recovery.time[cases[!test]] <- step
    cases <- c(new.cases, repeat.cases)
    if(verbose) cat('step: ', step, '\n')
    step <- step + 1
  }
  adf
}

adf$infection.time <- NA
adf$recovery.time <- NA
cases <- sample.int(nrow(adf), 10)
adf$infection.time[cases] <- 0

nsteps <- 38
step <- 1
tprob.sp <- 0.01
tprob.net <- 0.0
rprob <- 0.5
seasonal.factor <- function(x) -sinpi((x + 3)/ 52 * 2)
seasonal.amplitude <- .2

adf.out <- run.sims(adf)
step <- seq(1, to=nsteps)

events.by.state <- function(x, what) {
  tapply(adf.out[[what]] < x, adf.out$abb, sum, na.rm=TRUE)
}

cum.infections <- sapply(step, events.by.state, what='infection.time')
cum.recoveries <- sapply(step, events.by.state, what='recovery.time')
no.infected <- cum.infections - cum.recoveries
new.cases <- t(apply(cum.infections, 1, diff))
new.cases <- cbind(cum.infections[, 1], new.cases)

tmpf <- function(x, size=1.75, prep=.02) {
  xrep <- rbinom(x, size=x, prob=prep)
  rnbinom(n=xrep, size=size, mu=.76 + 1.92 *xrep)
}
reports <- t(apply(new.cases, 1, tmpf))

observed <- t(reports)
#observed <- t(new.cases)


## Mantel tests

observed <- observed[, colSums(observed) > 0]
pop.struct.mats <- MakePopStructMats(observed)
pop.dyn.mats <- MakePopDynMats(observed)
mantel.tests <- DoMantelTests(pop.dyn.mats, pop.struct.mats, permutations=1e3)

## Regression model

om <- GetRegressionData(case.data=as.data.frame(observed), flows=internal.flows,
                        fiAll=farms.and.inventory, stateCty=state.summaries)

m <- list()
f <- list()
f$lm <- as.formula(clog(cases) ~ clogCases1wa + logInternalFlowScaled + logCmedDenseScaled + weekCent + offset(log(nSusc1WksAgo) - 2*log(nFarms)))

f$lm <- as.formula(clog(cases) ~ clogCases1wa + logCmedDenseScaled + logInternalFlowScaled + weekCent + offset(log(nSusc1WksAgo) - 2*log(nFarms)))


#f$lm <- as.formula(clog(cases) ~ clogCases1wa + logInternalFlowScaled + logCmedDenseScaled + weekCent + offset(-log(nFarms)))

f$ct <- update(f$lm, cases ~ .)

m$lm <- lm(f$lm, data=om)
m$ct <- glm.nb(f$ct, data=om)
summary(m$lm)
summary(m$ct)
