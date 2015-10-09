#+setup, include=FALSE, cache=FALSE
library(knitr)
opts_chunk$set(fig.path='simulator-figures/', fig.align='center', fig.show='hold')
options(replace.assign=TRUE,width=80)
Sys.setlocale("LC_TIME", "C") #Needed for identical()
Sys.setlocale("LC_COLLATE", "C")

#' ## Load Packages

set.seed(4253)
library(plyr)
library(maptools)
library(raster)
library(rgdal)
print(sessionInfo())

#' ## Load data

load('flows-checkpoint1.RData')

ea.proj <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs" ## NAD83 Lambert Azimuthal Equal Area

## cb.proj <- "+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs" ## from .proj file
## cb <- readShapeSpatial("cb_2014_us_county_500k/cb_2014_us_county_500k")
##                         proj4string=CRS(cb.proj))
cb <- readOGR('cb_2014_us_county_500k', 'cb_2014_us_county_500k')

#' ## Sample random coordinates

cb.ea <- spTransform(cb, ea.proj)

tmpf <- function(geo=cb.ea){
  sf <- countyData$STFIPS
  cf <- countyData$COFIPS
  sf <- formatC(sf, flag="0", format="d", width=2)
  cf <- formatC(cf, flag="0", format="d", width=3)
  test1 <- geo@data[, 'STATEFP'] %in% sf
  test2 <- geo@data[, 'COUNTYFP'] %in% cf
  test <- test1 & test2
  geo[test, ]
}
cb2 <- tmpf()

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

coord.samps <- mapply(getSamp, cfps=countyData$COFIPS, sfps=countyData$STFIPS,
                      n=countyData$DATA, SIMPLIFY=FALSE)

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
stfips.rle <- list(lengths=countyData$DATA, values=countyData$STFIPS)
data(state.fips, package='maps')
key <- match(countyData$STFIPS, state.fips$fips)
countyData$abb <- as.character(state.fips$abb[key])

## create agent data structures

adf <- data.frame(cell=unlist(cell.samps),
                  abb=rep(countyData$abb, times=countyData$DATA),
                  infection.time=NA,
                  recovery.time=NA)
sp.nbs <- vector('list', length=nrow(adf))

occupied.cells <- unique(adf$cell)
cell2id <- lapply(occupied.cells, function(x) which(adf$cell == x))
names(cell2id) <- occupied.cells

## initialize infections

                                        #cases <- which(adf$cell %in% c(27799))
adf$infection.time <- NA
adf$recovery.time <- NA
cases <- sample.int(nrow(adf), 10)

get.nbs <- function(cell, rast=r) {
  nb <- adjacent(x=rast, cell, directions=8)
  nb.cells <- nb[, 'to']
  cell.names <- as.character(c(cell, nb.cells))
  unlist(cell2id[cell.names])
}
sp.nbs[cases] <- lapply(adf$cell[cases], get.nbs)
adf$infection.time[cases] <- 0

nsteps <- 20
step <- 1
tprob <- 0.01
rprob <- 0.5

while(step < nsteps){
  new.cases <- repeat.cases <- integer(0)
  for (case in cases){
    contacts <- sp.nbs[[case]]
    is.susceptible <- is.na(adf$infection.time[contacts])
    contacts <- contacts[is.susceptible]
    rand <- runif(n=length(contacts))
    test <- rand < tprob
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

## tabulation of case counts

step <- seq(1, to=nsteps)

events.by.state <- function(x, what) {
  tapply(adf[[what]] < x, adf$abb, sum, na.rm=TRUE)
}
cum.infections <- sapply(step, events.by.state, what='infection.time')
cum.recoveries <- sapply(step, events.by.state, what='recovery.time')
no.infected <- cum.infections - cum.recoveries

