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
  spdf[ind, ]
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
                    coord.df,
                    cofips=rep(county.hogs.pigs.02$COFIPS, times=n),
                    stfips=rep(county.hogs.pigs.02$STFIPS, times=n),
                    place.name=rep(county.hogs.pigs.02$GEO, times=n),
                    abb=rep(county.hogs.pigs.02$abb, times=n),
                    infection.time=NA,
                    recovery.time=NA)

  ## Generate lookup table of neighbors by space-------------------------
  occupied.cells <- unique(adf$cell)
  cell2id <- lapply(occupied.cells, function(x) which(adf$cell == x))
  names(cell2id) <- occupied.cells

  cell2nb.cells <- lapply(1:raster::ncell(raster.map), Get8nbs, r=raster.map)

  GetSpNbs <- function(cell, adj=cell2nb.cells, c2i=cell2id) {
    nb.cells <- adj[[cell]]
    cell.names <- as.character(nb.cells)
    unlist(c2i[cell.names])
  }
  sp.nbs <- lapply(adf$cell, GetSpNbs)

  ## Generate lookup table of neighbors by transport network---------------
  GetPrefMat <- function(rel='directed', flowMat, flows, state.tots) {
    nms <- state.tots$values
    fmo <- t(flowMat[nms, nms])
    diag(fmo) <- flows[nms, 'impInternalFlow']
    ## fmo[i, j] == head sent to i from j
    if(rel == 'directed') {
      F <- fmo
    } else {
      F <- fmo + t(fmo)
    }
    n <- state.tots$lengths
    pm <- F / n
    pm <- t(t(pm) / n)
    pm
  }

  adf <- adf[order(adf$abb), ]
  state.tots <- rle(as.character(adf$abb))
  pref.matrix <- GetPrefMat(flowMat=flows.matrix, flows=internal.flows,
                            state.tots=state.tots)
  pref.matrix <- pref.matrix / max(pref.matrix)
  block.sizes <- outer(state.tots$lengths, state.tots$lengths)
  mean.deg <- sum(block.sizes * pref.matrix) / sum(n)
  if(mean.deg >= target.mean.deg){
    pref.matrix <- pref.matrix * target.mean.deg / mean.deg
  } else {
    stop("target mean degree not possible")
  }
  trans.net <- igraph::sample_sbm(sum(state.tots$lengths), pref.matrix=pref.matrix,
                                  block.sizes=state.tots$lengths, directed=TRUE)
  net.nbs <- igraph::adjacent_vertices(trans.net, v=igraph::V(trans.net), mode='out')
  net.nbs <- sapply(net.nbs, as.integer)

  list(adf=adf, net.nbs=net.nbs, sp.nbs=sp.nbs)
}

RunSim <- function(adf, net.nbs, sp.nbs, nsteps=38, tprob.sp=0.01,
                    tprob.net=0.01, rprob=1, seasonal.amplitude=0,
                    verbose=FALSE, cases=1) {
  adf$infection.time <- NA
  adf$recovery.time <- NA
  adf$infection.time[cases] <- 0
  step <- 1
  while(step < nsteps){
    sf <- (1 + seasonal.amplitude * CalcSeasonalFactor(step))
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

SimulateAndSummarize <- function(job, static, dynamic, ...){
  adf.out <- RunSim(adf=dynamic$adf, net.nbs=dynamic$net.nbs,
                     sp.nbs=dynamic$sp.nbs, ...)
  events.by.state <- function(x, df=adf.out, what) {
    tapply(df[[what]] < x, df$abb, sum, na.rm=TRUE)
  }
  nsteps <- max(adf.out[, c('infection.time', 'recovery.time')], na.rm=TRUE)
  step <- seq(1, nsteps + 1)
  cum.infections <- sapply(step, events.by.state, what='infection.time')
  cum.recoveries <- sapply(step, events.by.state, what='recovery.time')
  no.infected <- cum.infections - cum.recoveries
  cum.infections <- cbind(0, cum.infections)
  new.cases <- t(apply(cum.infections, 1, diff))
  new.cases
}
