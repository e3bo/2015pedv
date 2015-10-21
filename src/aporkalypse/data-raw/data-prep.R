#+setup, include=FALSE, cache=FALSE
library(knitr)
opts_chunk$set(fig.path='data-prep-figures/', fig.align='center', fig.show='hold')
options(replace.assign=TRUE,width=80)
Sys.setlocale("LC_TIME", "C") #Needed for identical()
Sys.setlocale("LC_COLLATE", "C")

#+include=TRUE
print(sessionInfo())
dataDir <- file.path('.')
obs <- list()

GetCaseData <- function(){
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
  wk <- as.character(ret$week)
  ret$week <- as.Date(wk, format = "%m/%d/%Y")
  ret
}
obs$real.case.data <- GetCaseData()

GetSharedBorderAdjacencyMatrix <- function(){
    nbEdgelist <- read.csv(file.path(dataDir, 'state_neighbors_fips.txt'),
                           header=FALSE)
    data(state.fips, package='maps')
    g <- igraph::graph.data.frame(nbEdgelist, directed=FALSE)
    g <- igraph::simplify(g)
    key <- match(igraph::V(g)$name, state.fips$fips)
    abb <- state.fips$abb[key]
    igraph::V(g)$name <- as.character(abb)
    as.matrix(g[, ])
}
obs$shared.border.adjacency <- GetSharedBorderAdjacencyMatrix()

GetShipmentFlows <- function(){
    file <- file.path(dataDir, 'shipment-flows-origins-on-rows-dests-on-columns.csv')
    flows.matrix <- read.csv(file, row.names=1)
    flows.matrix <- data.matrix(flows.matrix)
    names(dimnames(flows.matrix)) <- c('origin', 'destination')

    stopifnot(rownames(flows.matrix) == colnames(flows.matrix))
    diag(flows.matrix) <- NA ## these are not part of original data

    ### assertions based on inspection of original .xls file
    stopifnot(flows.matrix['CA', 'IL'] == 9415)
    stopifnot(flows.matrix['MI', 'KS'] == 6)
    stopifnot(flows.matrix['IL', 'IA'] == 1424813)
    stopifnot(flows.matrix['KY', 'IA'] == 17658)
    stopifnot(flows.matrix['MO', 'IA'] == 2389932)
    stopifnot(flows.matrix['OK', 'CA'] == 16762)
    stopifnot(flows.matrix['MO', 'CA'] == 645)

    flows.matrix
  }
obs$flows.matrix <- GetShipmentFlows()

GetStateDists <- function(){
  x <- do.call(cbind, state.center)
  M <- sp::spDists(x, longlat=TRUE)
  colnames(M) <- rownames(M) <- state.abb
  M
}
obs$state.to.state.dists <- GetStateDists()

GetFarmsInventory <- function(){
    censPath <- file.path(dataDir, 'table19-2002.csv')
    cens <- read.csv(censPath, strip.white=TRUE, na.strings='(D)',
                     stringsAsFactors=FALSE)
    ## The introduction to the reports says '-' represents 0
    ## and (D) means deleted for privacy
    tmpf <- function(x) {
        if(is.character(x)){
            ret <- ifelse(x=='-', 0, x)
            type.convert(ret)
        }else{
            x
        }
    }
    cens <- lapply(as.list(cens), tmpf)
    cens <- data.frame(cens)
    test <- cens$GEO %in% state.name
    stateCounts <- cens[test, ]
    key <- match(stateCounts$GEO, state.name)
    rownames(stateCounts) <- state.abb[key]
    stateCounts
}
obs$farms.and.inventory <- GetFarmsInventory()

obs$balance.sheet.01.02 <- read.csv('state-hogBalanceSheetDec2000Dec2001.csv',
                                na.strings='-', row.names=1,
                                colClasses=c(state2='NULL'))

GetFarmsSales <- function(){
    censPath <- file.path(dataDir, 'table26-2002.csv')
    cens <- read.csv(censPath,strip.white=TRUE, na.strings='(D)',
                     stringsAsFactors=FALSE)
    ## The introduction to the reports says '-' represents 0
    ## and (D) means deleted for privacy
    tmpf <- function(x) {
        if(is.character(x)){
            ret <- ifelse(x=='-', 0, x)
            type.convert(ret)
        }else{
            x
        }
    }
    cens <- lapply(as.list(cens), tmpf)
    cens <- data.frame(cens)
    test <- cens$GEO %in% state.name
    stateCounts <- cens[test, ]
    key <- match(stateCounts$GEO, state.name)
    stateCounts$abb <- state.abb[key]
    stateCounts
}
obs$farms.and.sales <- GetFarmsSales()

GetInternalFlows <- function(balanceSheet, flowMat, farmsSales) {
  plotPredDemandSales <- function(feederDemand, deaths, estFeederProduction,
                                  imp, exp){
    par(mfrow=c(3,2))
    internalSupply <- estFeederProduction - exp
    demand <- feederDemand + deaths
    internalDemand <- demand - imp
    plot(internalSupply, internalDemand)
    abline(0,1)
    plot(internalSupply, internalDemand, log='xy')
    text(internalSupply, internalDemand, labels=names(exp))
    abline(0,1)
    res <- log(internalSupply) - log(internalDemand)
    barplot(res, names.arg=names(exp), las=2)
    qqnorm(res)
    predInternalFlow <- exp(log(internalDemand) + res/2)
    plot(internalSupply, predInternalFlow, log='xy')
    text(internalSupply, predInternalFlow, labels=names(exp))
    abline(0,1)
    plot(internalDemand, predInternalFlow, log='xy')
    text(internalDemand, predInternalFlow, labels=names(exp))
    abline(0,1)
    par(mfrow=c(1,1))
    cbind(internalSupply, internalDemand, predInternalFlow, exports=exp,
          imports=imp)
  }
  deaths <- balanceSheet[-51,]$deaths*1000
  labs <- rownames(balanceSheet)[1:50]
  diag(flowMat) <- 0
  imports <- colSums(flowMat)[labs]
  exports <- rowSums(flowMat)[labs]
  feederDemand <- rowSums(cbind(farmsSales$Finish.only.Number,
                                farmsSales$Nursery.number),
                          na.rm=FALSE)
  deaths <- balanceSheet[-51, ]$deaths * 1000
  estFeederProduction <- rowSums(cbind(farmsSales$Farrow.to.wean.Number,
                                       farmsSales$Farrow.to.feeder.Number,
                                       farmsSales$Nursery.Number), na.rm=FALSE)
  flows <- plotPredDemandSales(feederDemand=feederDemand, deaths=deaths,
                             estFeederProduction=estFeederProduction, imp=imports,
                               exp=exports)
  imp <- ifelse(is.na(flows[ ,'internalSupply']) | flows[, 'internalSupply'] <= 0,
                flows[, 'internalDemand'], flows[, 'internalSupply'])
  imp <- ifelse(is.na(flows[, 'predInternalFlow']), imp, flows[, 'predInternalFlow'])
  cbind(flows, impInternalFlow=imp)
}
obs$internal.flows <- GetInternalFlows(balanceSheet=obs$balance.sheet.01.02,
                                       flowMat=obs$flows.matrix,
                                       farmsSales=obs$farms.and.sales)

PlotFlows <- function(flows) {
  totalFlows <- flows[, c('exports', 'imports', 'impInternalFlow')] %*% c(1, 1, 2)
  ## Multiply internal flows by 2 since they are both imports and exports
  barplot(sort(totalFlows[,1]), log='y', ylab='Head', xlab='State', las=2)
    ord <- order(totalFlows)
    sel <- rownames(totalFlows)[ord]
    cols <- c('impInternalFlow', 'exports', 'imports')
    h <- flows[sel, cols]
    h <- h %*% diag(c(2, 1, 1))
    h <- h[complete.cases(h), ]
    colnames(h) <- cols
  barplot(t(h), las=2, legend.text=T, args.legend=list(x="topleft"), ylab='Head',
          xlab='State')
}
PlotFlows(flows=obs$internal.flows)

obs$county.areas <- read.delim('2013_Gaz_counties_national.txt')
obs$county.areas <- obs$county.areas[, c('GEOID', 'ALAND')]

GetCountyHogs <- function(ITEM.pattern){
  data(county.fips, package='maps')
  cens <- read.csv('table-12-hogs-and-pigs-by-county.csv',strip.white=TRUE,
                   na.strings='(D)', stringsAsFactors=FALSE)
  ## The introduction to the reports says '-' represents 0
  ## and (D) means deleted for privacy
  tmpf <- function(x) {
      if(is.character(x)){
          ret <- ifelse(x=='-', 0, x)
          type.convert(ret)
      }else{
          x
      }
  }
  cens <- lapply(as.list(cens), tmpf)
  cens <- data.frame(cens)
  test <- cens$STCOFIPS %in% county.fips$fips
  cens <- cens[test,]
  ind <- grep(ITEM.pattern, cens$ITEM)
  cens <- cens[ind,]
  cens
}
obs$county.hogs.pigs <- GetCountyHogs(ITEM.pattern='farms, 2007)$')
obs$county.hogs.pigs.02 <- GetCountyHogs(ITEM.pattern='farms, 2002)$')

obs$resource.regs <- read.csv('reglink.csv', skip=2,
                              colClasses=c(NA, NA, 'NULL', 'NULL'))

GetStateSummaries <- function(regs, rd, countyData) {
  ctyTots <- plyr::ddply(countyData, 'STCOFIPS', plyr::summarize,
                         totalFarms=sum(DATA), STFIPS=STFIPS[1],
     smallFarms = sum(DATA["Inventory \\ Total hogs and pigs \\ Farms by inventory \\ 1 to 24 (farms, 2007)" == ITEM]))
    ctyTots$nonSmallFarms <- with(ctyTots, totalFarms - smallFarms)
    mg <- merge(ctyTots, rd, by.x='STCOFIPS', by.y='GEOID')
    mg <- merge(regs, mg, by.x='Fips', by.y='STCOFIPS')
    mg$kmsq <- mg$ALAND/1e6
    mg$nonSmallDense <- mg$nonSmallFarms / mg$kmsq
    mg$farmDense <- mg$totalFarms / mg$kmsq
  stateCty <- plyr::ddply(mg, 'STFIPS', plyr::summarize,
                          totNonSmall=sum(nonSmallFarms),
                            mDense=mean(farmDense),
                            cmDense=mean(farmDense[totalFarms > 0]),
                            medDense=median(farmDense),
                            cmedDense=median(farmDense[totalFarms > 0]),
                            reg1=weighted.mean(ERS.resource.region==1, w=totalFarms),
                            reg2=weighted.mean(ERS.resource.region==2, w=totalFarms),
                            reg3=weighted.mean(ERS.resource.region==3, w=totalFarms),
                            reg4=weighted.mean(ERS.resource.region==4, w=totalFarms),
                            reg5=weighted.mean(ERS.resource.region==5, w=totalFarms),
                            reg6=weighted.mean(ERS.resource.region==6, w=totalFarms),
                            reg7=weighted.mean(ERS.resource.region==7, w=totalFarms),
                            reg8=weighted.mean(ERS.resource.region==8, w=totalFarms),
                            reg9=weighted.mean(ERS.resource.region==9, w=totalFarms))
    data(state.fips, package='maps')
    key <- match(stateCty$STFIPS, state.fips$fips)
    stateCty$state <- state.fips$abb[key]
    stateCty
}
obs$state.summaries <- GetStateSummaries(regs=obs$resource.regs, rd=obs$county.areas,
                                         countyData=obs$county.hogs.pigs)

GetMap <- function(countyData){
  ea.proj <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs"
  ## NAD83 Lambert Azimuthal Equal Area
  cb <- rgdal::readOGR('cb_2014_us_county_500k', 'cb_2014_us_county_500k')
  cb.ea <- sp::spTransform(cb, sp::CRS(ea.proj))

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
  tmpf(cb.ea, countyData)
}
#obs$county.hogs.pigs.02.map <- GetMap(countyData=obs$county.hogs.pigs.02)

if(!file.exists('data'))
  dir.create('data')
setwd('data')
nms <- names(obs)
files <- paste0(nms, '.rda')
env <- as.environment(obs)
tmpf <- function(x, y) save(list=x, file=y, envir=env)
mapply(tmpf, nms, files)
setwd('..')
