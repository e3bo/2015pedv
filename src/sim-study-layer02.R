#!/usr/bin/Rscript
library(methods) #for raster

set.seed(123, "L'Ecuyer")

mc.cores <- ifelse(Sys.info()['sysname'] == "Linux",
                   parallel::detectCores() - 1, 1)
mc.cores <- ifelse(mc.cores > 20, 20, mc.cores)
options('mc.cores'=mc.cores)

load('sim-study-checkpoint2.rda')

df <- do.call(cbind, c(list(des), as.list(ag.data.inds)))
df <- as.data.frame(df)

df$raster.cell.side <- raster.cell.side.grid[df$ag.ind]
df$target.mean.deg <- target.mean.deg.grid[df$net.nbs.ind]
df$lags.sel <- 1
df$nstarters <- 1
df$permutations <- 1000
df$starting.grid.nx <- 10
df$starting.grid.ny <- 2
Wrapper <- function(...) try(sds::SimulateAndSummarize(...))

system.time(res <- parallel::mcMap(Wrapper,
                                   agent.data=ag[df$ag.ind],
                                   lags.sel=df$lags.sel,
                                   net.nbs=net.nbs[df$net.nbs.ind],
                                   nstarters=df$nstarters,
                                   permutations=df$permutations,
                                   prep=df$prep,
                                   rprob=df$rprob,
                                   size=df$size,
                                   seasonal.amplitude=df$seasonal.amplitude,
                                   starting.grid.nx=df$starting.grid.nx,
                                   starting.grid.ny=df$starting.grid.ny,
                                   starting.grid.x=df$starting.grid.x,
                                   starting.grid.y=df$starting.grid.y,
                                   tprob.net=df$tprob.net,
                                   tprob.outside=df$tprob.outside,
                                   tprob.sp=df$tprob.sp))
recs <- lapply(res, '[[', 'record')
recs <- do.call(rbind, recs)
resall <- cbind(df, recs)
save.image('sim-study-checkpoint3.rda')

sub <- resall
sub$r <- sub[ , paste0('mantel.r.lag1.shipment.spearman.TRUE.', df$permutations[1])]

GetVaryingInputs <- function(df){
  rngs <- apply(df, 2, range)
  is.varying <- apply(rngs, 2, diff) > sqrt(.Machine$double.eps)
  ret <- names(is.varying[is.varying])
  is.ind <- grepl('\\.ind$', ret)
  ret[!is.ind]
}
RunKriging <- function(df, design){
  km.vars <- GetVaryingInputs(df)
  rhs <- do.call(paste, c(as.list(km.vars), list(sep="+")))
  formula <- as.formula(paste('~', rhs))
  DiceKriging::km(formula=formula, design=design[, km.vars], response=design$r,
                  nugget.estim=TRUE, covtype='matern3_2',
                  control=list(maxit=1000))
}

m <- RunKriging(df=df, design=sub)
save.image('sim-study-checkpoint4.rda')

nmeta <- 1e5
extra.par.ranges <- list(target.mean.deg=range(target.mean.deg.grid),
                         raster.cell.side=range(raster.cell.side.grid))
all.par.ranges <- c(par.ranges, extra.par.ranges)

GetRandLHSDes <- function(n, ranges){
  X <- lhs::randomLHS(n, length(ranges))
  tmpf <- function(samp, range){
    d <- diff(sort(range))
    samp <- samp * d
    samp + range[1]
  }
  for(i in 1:ncol(X)){
    X[, i] <- tmpf(X[ ,i], ranges[[i]])
  }
  colnames(X) <- names(ranges)
  X
}

RunSobol <- function(nmeta, km, all.par.ranges, order=1){
  vars <- km@covariance@var.names
  rngs <- all.par.ranges[vars]
  X1 <- GetRandLHSDes(nmeta, rngs)
  X2 <- GetRandLHSDes(nmeta, rngs)
  colnames(X1) <- colnames(X2) <- vars
  Wrapper <- function(X){
    predict(km, newdata=X, type='UK')$m
  }
  sob <- sensitivity::sobol(model=Wrapper, X1=X1, X2=X2, order=order, nboot=100)
  vym <- sob$V['global', 'original']
  vyd <- DiceKriging::coef(m)$sd2 + DiceKriging::coef(m)$nugget
  vy <- vym + vyd
  sob.inds <- sob$S[, 'original'] * vym / vy
  names(sob.inds) <- rownames(sob$S)
  sob.ind.rand <- vyd / vy
  list(sob, sob.inds, sob.ind.rand)
}

(sob.out <- RunSobol(nmeta, km=m, all.par.ranges=all.par.ranges, order=1))

save.image('sim-study-checkpoint5.rda')

MakePlots <- function(km, all.par.ranges, npoints=1000, sub, df, v1='tprob.net',
                      v2='tprob.sp'){
  vars <- km@covariance@var.names
  rngs <- all.par.ranges[vars]
  X1 <- GetRandLHSDes(npoints, rngs)
  colnames(X1) <- vars
  p <- predict(m, newdata=X1, type='UK')
  x <- X1[, v1]
  plot(x, p$mean, ylim=c(0, 0.5), xlab=v1)
  points(x, p$lower95, col='grey')
  points(x, p$upper95, col='grey')
  points(sub$tprob.net, sub$r, col=2)
  plotdes <- sensitivity::parameterSets(par.ranges=all.par.ranges[c(v1, v2)],
                                      samples=sqrt(npoints), method='grid')
  des <- X1
  des[, v1] <- plotdes[1:nrow(des), v1]
  des[, v2] <- plotdes[1:nrow(des), v2]
  pmean <- predict(m, newdata=des, type='UK')$m
  f <- as.formula(paste('pmean~', paste(v1, v2, sep="*")))
  lattice::levelplot(f, data=as.data.frame(des))
}

MakePlots(km=m, all.par.ranges=all.par.ranges, sub=sub, df=df)

save.image('sim-study-checkpoint6.rda')