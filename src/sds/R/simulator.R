#' Get neighbors of a raster cell using 8-cell neighborhood
#'
#' @param r A raster object.
#' @param cell A cell in the raster object.
#' @param ncells Number of cells in raster object.
#' @param nr Number of rows in raster object.
#' @param nc Number of columns in raster object.
#' @return A vector of the neighboring cells of \code{cell} as well as \code{cell}.
#'
Get8nbs <- function(r, cell, ncells, nc, nr){
  if (cell > ncells) {
    return(NA)
  }
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

#' Get spatial polygons data frame corresponding to a county
#' @param sfps State fips code
#' @param cfps County fips code
#' @param spdf Spatial polygons data frame with all counties
#'
#' @import methods
GetCountySPDF <- function(sfps, cfps, spdf){
  cfps <- FormatFips(cfps, 'county')
  sfps <- FormatFips(sfps, 'state')
  ind <- which(spdf@data[, 'COUNTYFP']==cfps & spdf@data[, 'STATEFP']==sfps)
  loadNamespace('sp') ## needed for next line to work 1st time if not already loaded
  spdf[ind, ]
}

SubsetFlows <- function(nms, rel='directed'){
  data('flows.matrix', envir=environment(), package='sds')
  flows.matrix <- get('flows.matrix')
  data('internal.flows', envir=environment(), package='sds')
  internal.flows <- get('internal.flows')
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

#' @export
GetNetNbs <- function(block.labels, target.mean.deg){
  state.tots <- rle(as.character(block.labels))
  stopifnot(anyDuplicated(state.tots$values)==0)
  pref.matrix <- GetPrefMat(state.tots=state.tots,
                            target.mean.deg=target.mean.deg)
  n <- state.tots$lengths
  trans.net <- igraph::sbm.game(sum(n), pref.matrix=pref.matrix,
                                block.sizes=n, directed=TRUE)
  net.nbs <- igraph::get.adjlist(trans.net, mode='out')
  sapply(net.nbs, as.integer)
}

#' Create agent data structures
#' @export
CreateAgents <- function(raster.cell.side.meters=16000, census.dilation=1,
                         target.mean.deg=1){

  state.fips <- maps::state.fips
  data('county.hogs.pigs.02', package='sds', envir=environment())
  chp <- get('county.hogs.pigs.02')
  key <- match(chp$STFIPS, state.fips$fips)
  chp$abb <- as.character(state.fips$abb[key])

  ## Sample coordinates of farms-------------------------------------
  data('county.hogs.pigs.02.map', package='sds', envir=environment())
  map <- get('county.hogs.pigs.02.map')
  GetSampCoord <- function(cfps, sfps, n, spdf=map,
                           plot.samp=FALSE){
    if(n == 0) {
        samp <- NULL
      } else {
        cty.poly <- GetCountySPDF(sfps, cfps, spdf=spdf)
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
  n <- floor(chp$DATA * census.dilation)
  coord.samps <- mapply(GetSampCoord, cfps=chp$COFIPS, sfps=chp$STFIPS,
                        n=n, SIMPLIFY=FALSE)
  coord.mats <- lapply(coord.samps, function(x) if(!is.null(x)) sp::coordinates(x))
  coord.df <- do.call(rbind, coord.mats)

  ## Convert coordinates to cell membership in raster layer------------
  raster.map <- raster::raster(map)
  raster::res(raster.map) <- raster.cell.side.meters

  GetCell <- function(xy, r=raster.map) {
    if(is.null(xy)) {
      NULL
    } else {
      samps <- raster::cellFromXY(object=r, xy=xy)
      if(any(is.na(samps))) stop("Some farms outside of raster cells")
      samps
    }
  }
  cell.samps <- lapply(coord.samps, GetCell)

  ## Generate agent data frame----------------------------------------
  adf <- data.frame(cell=unlist(cell.samps),
                    x=coord.df[ ,1],
                    y=coord.df[, 2],
                    cofips=rep(chp$COFIPS, times=n),
                    stfips=rep(chp$STFIPS, times=n),
                    place.name=rep(chp$GEO, times=n),
                    abb=rep(chp$abb, times=n),
                    infection.time=NA,
                    recovery.time=NA)
  adf <- adf[order(adf$stfips), ]
  adf$id <- 1:nrow(adf)

  ## Generate lookup table of neighbors by space-------------------------
  occupied.cells <- unique(adf$cell)
  cell2id <- lapply(occupied.cells, function(x) adf$id[adf$cell == x])
  names(cell2id) <- occupied.cells

  cell2nb.cells <- lapply(1:raster::ncell(raster.map), Get8nbs, r=raster.map,
                          ncells=raster::ncell(raster.map),
                          nc=raster::ncol(raster.map),
                          nr=raster::nrow(raster.map))

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
                   tprob.net=0.01, tprob.outside=0, rprob=1,
                   seasonal.amplitude=0, verbose=FALSE, cases=1) {
  adf$infection.time <- NA
  adf$recovery.time <- NA
  adf$infection.time[cases] <- 0
  step <- 1
  adf <- adf[order(adf$id), ] ## Assumed for looking up neighbors
  nagents <- nrow(adf)
  while(step < nsteps){
    sf <- (1 + seasonal.amplitude * CalcSeasonalFactor(step))
    net.contact.dist <- table(unlist(net.nbs[cases]))
    sp.contact.dist <- table(unlist(sp.nbs[cases]))
    tprob.sp.t <- max(min(tprob.sp * sf, 1), 0)
    tprob.net.t <- max(min(tprob.net * sf, 1), 0)
    tprob.outside.t <- max(min(tprob.outside * sf, 1), 0)

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

    n.new.cases.outside <- rbinom(1, size=nagents, prob=tprob.outside.t)
    new.cases.outside <- sample.int(n=nagents, size=n.new.cases.outside)

    new.cases <- unique(c(new.cases.net, new.cases.sp, new.cases.outside))
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

GetTimeSeries <- function(adf, nsteps) {
  step <- seq(1, nsteps + 1)
  cum.infections <- sapply(step, events.by.state, df=adf, what='infection.time')
  cum.recoveries <- sapply(step, events.by.state, df=adf, what='recovery.time')
  no.infected <- cum.infections - cum.recoveries
  foo <- cbind(0, cum.infections)
  new.cases <- t(apply(foo, 1, diff))
  list(new.cases=new.cases, no.infected=no.infected)
}

NewCasesToReports <- function(x, size=1 / 1.75, prep=.5){
  xrep <- rbinom(length(x), size=x, prob=prep)
  rnbinom(n=length(xrep), size=size, mu=0.76 + 1.92 * xrep)
}

#' @export
SimulateAndSummarize <- function(agent.data,
                                 lags.sel=c(0, 1),
                                 net.nbs,
                                 nstarters=1,
                                 nsteps=38,
                                 permutations=1e3,
                                 prep=0.5,
                                 starting.state=NULL,
                                 starting.grid.nx=10,
                                 starting.grid.ny=2,
                                 starting.grid.x=NULL,
                                 starting.grid.y=NULL,
                                 size=0.75,
                                 verbose=TRUE,
                                  ...){
  if (!is.null(starting.state)){
    possible.starters <- which(agent.data$adf$abb == starting.state)
  } else if (!is.null(starting.grid.x) & !is.null(starting.grid.y)){
    cx <- as.integer(cut(agent.data$adf$x, breaks=starting.grid.nx))
    cy <- as.integer(cut(agent.data$adf$y, breaks=starting.grid.ny))
    thex <- ceiling(starting.grid.x * starting.grid.nx)
    they <- ceiling(starting.grid.y * starting.grid.ny)
    test <- cx == thex & cy == they
    possible.starters <- which(test)
  } else {
    possible.starters <- seq(1, nrow(agent.data$adf))
  }
  if(length(possible.starters) == 0){
    stop("No farms that meet requirements to start epizootic")
  }
  inds <- sample.int(length(possible.starters), size=nstarters)
  cases <- possible.starters[inds]
  adf.out <- RunSim(adf=agent.data$adf, net.nbs=agent.data$net.nbs,
                    sp.nbs=agent.data$sp.nbs, cases=cases, nsteps=nsteps, ...)
  tsl <- GetTimeSeries(adf.out, nsteps=nsteps)
  observed <- apply(tsl$new.cases, 1, NewCasesToReports, prep=prep, size=size)
  is.state.count.constant <- apply(observed, 2, var) < sqrt(.Machine$double.eps)
  n.constant.positive <- sum(observed[1, is.state.count.constant] > 0)
  observed <- observed[, !is.state.count.constant, drop=FALSE]
  if(ncol(observed) == 1 | nsteps < 3 | n.constant.positive > 0){
    des <- MakeDesign(mat.list1, mat.list2, ...)
    list(mantel.tests=NA, mantel.observed=observed, record=NA)
  } else {
    pop.struct.mats <- MakePopStructMats(observed)
    pop.dyn.mats <- MakePopDynMats(observed, lags.sel=lags.sel)
    mantel.tests <- DoMantelTests(pop.dyn.mats, pop.struct.mats,
                                  permutations=permutations)
    record <- SpreadMantelTestsDf(mantel.tests)
    list(mantel.tests=mantel.tests, mantel.observed=observed, record=record)
  }
}

SpreadMantelTestsDf <- function(df){
  r <- as.list(df$r)
  p.values <- as.list(df$pValues)
  dep.vars <- !grepl('^r$|^pValues$', colnames(df))
  nms <- as.character(interaction(df[, dep.vars]))
  names(r) <- names(p.values) <- nms
  cbind(mantel.r=as.data.frame(r), mantel.p=as.data.frame(p.values))
}

SSRealCases <- function(mts){
  r <- real.case.data[, -c(1,2)]
  nrow <- max(nrow(r), nrow(mts))
  vars <- union(colnames(r), colnames(mts))
  M2 <- M1 <- matrix(0, nrow=nrow, ncol=length(vars))
  colnames(M2) <- colnames(M1) <- vars
  pad <- integer(nrow - nrow(r))
  for(v in colnames(r)){
    M1[, v] <- c(r[, v], pad)
  }
  pad <- integer(nrow - nrow(mts))
  for(v in colnames(mts)){
    M2[, v] <- c(mts[, v], pad)
  }
  resss <- sum(as.numeric(M1 - M2)^2)
  list(resss=resss, observed=M1, modeled=M2)
}
