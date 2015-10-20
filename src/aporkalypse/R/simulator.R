
Get8nbs <- function(r, cell){
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

## TODO DATA to generate cb2 -> county.hogs.pigs.02.map
## review pref matrix derivation for sensibility

CreateAgents <- function(job, static,
                         raster.cell.side.meters=16000, ...){
  data(county.hogs.pigs.02)
  data(county.hogs.pigs.map)
  data(state.fips, package='maps')

  key <- match(county.hogs.pigs.02$STFIPS, state.fips$fips)
  county.hogs.pigs.02$abb <- as.character(state.fips$abb[key])

  ## Sample coordinates of farms
  GetSampCoord <- function(spdf=county.hogs.pigs.02.map, cfps, sfps, n,
                           plot.samp=FALSE){
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
  n <- county.hogs.pigs.02$DATA
  coord.samps <- mapply(getSamp, cfps=county.hogs.pigs.02$COFIPS,
                        sfps=county.hogs.pigs.02$STFIPS,
                        n=n, SIMPLIFY=FALSE)

  ## Convert coordinates to cell membership in raster layer

  raster.map <- raster(county.hogs.pigs.02.map)
  res(raster.map) <- raster.cell.side.meters

  tmpf <- function(xy, r=raster.map) {
    if(is.null(xy)) {
      NULL
    } else {
      cellFromXY(object=r, xy=xy)
    }
  }
  cell.samps <- lapply(coord.samps, tmpf)

  ## Generate agent data frame
  adf <- data.frame(cell=unlist(cell.samps),
                    abb=rep(county.hogs.pigs.02$abb, times=n),
                    infection.time=NA,
                    recovery.time=NA)

  ## Generate lookup table of neighbors by space
  occupied.cells <- unique(adf$cell)
  cell2id <- lapply(occupied.cells, function(x) which(adf$cell == x))
  names(cell2id) <- occupied.cells

  tmpf <- function(x, r=raster.map) my.8.nbs(r=r, x)
  cell2nb.cells <- lapply(1:ncell(r), tmpf)

  GetSpNbs <- function(cell, adj=cell2nb.cells, c2i=cell2id) {
    nb.cells <- adj[[cell]]
    cell.names <- as.character(nb.cells)
    unlist(c2i[cell.names])
  }
  sp.nbs <- lapply(adf$cell, GetSpNbs)

  ## Generate lookup table of neighbors by transport network





  res
}

RunSim <- function(adf, net.nbs, sp.nbs, nsteps=38, trpob.sp=0.01,
                    tprob.net=0.01, rprob=1, seasonal.amplitude=0,
                    verbose=FALSE, cases=1) {
  adf$infection.time <- NA
  adf$recovery.time <- NA
  adf$infection.time[cases] <- 0
  step <- 1
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

SimulateAndSummarize <- function(job, static, dynamic, ...){
  adf.out <- RunSim(adf=dynamic$adf, net.nbs=dynamic$net.nbs,
                     sp.nbs=dynamic$sp.nbs, ...)
  events.by.state <- function(x, df=adf.out, what) {
    tapply(df[[what]] < x, df$abb, sum, na.rm=TRUE)
  }
  cum.infections <- sapply(step, events.by.state, what='infection.time')
  cum.recoveries <- sapply(step, events.by.state, what='recovery.time')
  no.infected <- cum.infections - cum.recoveries
  new.cases <- t(apply(cum.infections, 1, diff))
  new.cases <- cbind(cum.infections[, 1], new.cases)
  new.cases
}
