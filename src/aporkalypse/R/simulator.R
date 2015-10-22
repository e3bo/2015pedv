#' Get neighbors of a raster cell using 8-cell neighborhood
#'
#' @param r A raster object.
#' @param cell A cell in the raster object.
#' @return A vector of the neighboring cells of \code{cell} as well as \code{cell}.
#'
Get8nbs <- function(r, cell){
  if (cell > raster::ncell(r)) {
    return(NA)
  }
  nc <- raster::ncol(r)
  nr <- raster::nrow(r)
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

CalcSeasonalFactor <- function(x) -sinpi((x + 3)/ 52 * 2)

FormatFips <- function(x, type) {
  if (type=='county'){
    formatC(x, flag="0", format="d", width=3)
  } else if (type=='state'){
    formatC(x, flag="0", format="d", width=2)
  } else NA
}

GetCountySPDF <- function(sfps, cfps){
  spdf <- county.hogs.pigs.02.map
  cfps <- FormatFips(cfps, 'county')
  sfps <- FormatFips(sfps, 'state')
  ind <- which(spdf@data[, 'COUNTYFP']==cfps & spdf@data[, 'STATEFP']==sfps)
  loadNamespace('sp') ## needed for next line to work 1st time if not already loaded
  spdf[ind, ]
}

SubsetFlows <- function(nms, rel='directed'){
  fmo <- flows.matrix[nms, nms]
  diag(fmo) <- internal.flows[nms, 'impInternalFlow']
  ## fmo[i, j] == head sent to i from j
  if(rel == 'directed') {
    F <- fmo
  } else {
    F <- fmo + t(fmo)
  }
  F
}

GetPrefMat <- function(state.tots, target.mean.deg) {
  nms <- state.tots$values
  F <- SubsetFlows(nms)
  n <- state.tots$lengths
  block.sizes <- outer(n, n)
  pm <- F / block.sizes
  pm <- pm / max(pm) # Probs must be <=1
  mean.deg <- sum(block.sizes * pm) / sum(n)
  if(mean.deg >= target.mean.deg){
    pm <- pm * target.mean.deg / mean.deg
  } else {
    stop("Target mean degree not possible")
  }
  pm
}

GetNetNbs <- function(block.labels, target.mean.deg){
  state.tots <- rle(as.character(block.labels))
  stopifnot(anyDuplicated(state.tots$values)==0)
  pref.matrix <- GetPrefMat(state.tots=state.tots,
                            target.mean.deg=target.mean.deg)
  n <- state.tots$lengths
  trans.net <- igraph::sample_sbm(sum(n), pref.matrix=pref.matrix,
                                  block.sizes=n, directed=TRUE)
  net.nbs <- igraph::adjacent_vertices(trans.net, v=igraph::V(trans.net),
                                       mode='out')
  sapply(net.nbs, as.integer)
}

CreateAgents <- function(job, static,
                         raster.cell.side.meters=16000,
                         census.dilation=1,
                         target.mean.deg=1, ...){

  data('state.fips', package='maps', envir=environment())
  key <- match(county.hogs.pigs.02$STFIPS, state.fips$fips)
  county.hogs.pigs.02$abb <- as.character(state.fips$abb[key])

  ## Sample coordinates of farms-------------------------------------
  GetSampCoord <- function(cfps, sfps, n,
                           plot.samp=FALSE){
    if(n == 0) {
        samp <- NULL
      } else {
        cty.poly <- GetCountySPDF(sfps, cfps)
       len <- 0
       while(len != n){
         samp <- sp::spsample(cty.poly, type='random', n=n, iter=10)
         len <- length(samp)
       }
       if(plot.samp){
          plot(cty.poly)
          points(samp)
       }
    }
    samp
  }
  n <- floor(county.hogs.pigs.02$DATA * census.dilation)
  coord.samps <- mapply(GetSampCoord, cfps=county.hogs.pigs.02$COFIPS,
                        sfps=county.hogs.pigs.02$STFIPS,
                        n=n, SIMPLIFY=FALSE)
  coord.mats <- lapply(coord.samps, function(x) if(!is.null(x)) sp::coordinates(x))
  coord.df <- do.call(rbind, coord.mats)

  ## Convert coordinates to cell membership in raster layer------------
  raster.map <- raster::raster(county.hogs.pigs.02.map)
  raster::res(raster.map) <- raster.cell.side.meters

  GetCell <- function(xy, r=raster.map) {
    if(is.null(xy)) {
      NULL
    } else {
      raster::cellFromXY(object=r, xy=xy)
    }
  }
  cell.samps <- lapply(coord.samps, GetCell)

  ## Generate agent data frame
  adf <- data.frame(cell=unlist(cell.samps),
                    x=coord.df[ ,1],
                    y=coord.df[, 2],
                    cofips=rep(county.hogs.pigs.02$COFIPS, times=n),
                    stfips=rep(county.hogs.pigs.02$STFIPS, times=n),
                    place.name=rep(county.hogs.pigs.02$GEO, times=n),
                    abb=rep(county.hogs.pigs.02$abb, times=n),
                    infection.time=NA,
                    recovery.time=NA)
  adf <- adf[order(adf$stfips), ]
  adf$id <- 1:nrow(adf)

  ## Generate lookup table of neighbors by space-------------------------
  occupied.cells <- unique(adf$cell)
  cell2id <- lapply(occupied.cells, function(x) adf$id[adf$cell == x])
  names(cell2id) <- occupied.cells

  cell2nb.cells <- lapply(1:raster::ncell(raster.map), Get8nbs, r=raster.map)

  GetSpNbs <- function(cell, adj=cell2nb.cells, c2i=cell2id) {
    nb.cells <- adj[[cell]]
    cell.names <- as.character(nb.cells)
    unlist(c2i[cell.names])
  }
  sp.nbs <- lapply(adf$cell, GetSpNbs)

  net.nbs <- GetNetNbs(adf$abb, target.mean.deg)

  list(adf=adf, net.nbs=net.nbs, sp.nbs=sp.nbs)
}

RunSim <- function(adf, net.nbs, sp.nbs, nsteps=38, tprob.sp=0.01,
                    tprob.net=0.01, rprob=1, seasonal.amplitude=0,
                    verbose=FALSE, cases=1) {
  adf$infection.time <- NA
  adf$recovery.time <- NA
  adf$infection.time[cases] <- 0
  step <- 1
  adf <- adf[order(adf$id), ] ## Assumed for looking up neighbors
  while(step < nsteps){
    sf <- (1 + seasonal.amplitude * CalcSeasonalFactor(step))
    net.contact.dist <- table(unlist(net.nbs[cases]))
    sp.contact.dist <- table(unlist(sp.nbs[cases]))
    tprob.sp.t <- max(min(tprob.sp * sf, 1), 0)
    tprob.net.t <- max(min(tprob.net * sf, 1), 0)

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

events.by.state <- function(x, df, what) {
  tapply(df[[what]] < x, df$abb, sum, na.rm=TRUE)
}

GetTimeSeries <- function(adf) {
  nsteps <- max(adf[, c('infection.time', 'recovery.time')], na.rm=TRUE)
  step <- seq(1, nsteps + 1)
  cum.infections <- sapply(step, events.by.state, df=adf, what='infection.time')
  cum.recoveries <- sapply(step, events.by.state, df=adf, what='recovery.time')
  no.infected <- cum.infections - cum.recoveries
  foo <- cbind(0, cum.infections)
  new.cases <- t(apply(foo, 1, diff))
  list(new.cases=new.cases, no.infected=no.infected)
}

SimulateAndSummarize <- function(job, static, dynamic,
                                 starting.state='OH',
                                 nstarters=1, verbose=TRUE, ...){
  possible.starters <- which(dynamic$abb == starting.state)
  inds <- sample.int(length(possible.starters), size=nstarters)
  cases <- possible.starters[inds]
  adf.out <- RunSim(adf=dynamic$adf, net.nbs=dynamic$net.nbs,
                    sp.nbs=dynamic$sp.nbs, cases=cases, ...)
  tsl <- GetTimeSeries(adf.out)
  observed <- t(tsl$new.cases)
  isStateAffected <- colSums(observed) > 0
  if(sum(isStateAffected) == 1 | nrow(observed) < 3){
    return(NA)
  } else {
    observed <- observed[, isStateAffected]
    pop.struct.mats <- MakePopStructMats(observed)
    pop.dyn.mats <- MakePopDynMats(observed)
    mantel.tests <- DoMantelTests(pop.dyn.mats, pop.struct.mats, permutations=1e3)
    list(mantel.tests=mantel.tests, observed=mantel.observed)
  }
}

