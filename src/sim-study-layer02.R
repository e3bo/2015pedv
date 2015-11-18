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
df$permutations <- 2
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

X <- sub[, GetVaryingInputs(df)]
Y <- sub$r
ntest <- 0.5 * nrow(X)
test.ind <- seq_len(ntest)
X.test <- X[test.ind, ]
Y.test <- Y[test.ind]
X.train <- X[-test.ind, ]
Y.train <- Y[-test.ind]
testdf <- cbind(X.test, Y.test)

pdf('modelComparisonc-plots.pdf')
mc <- modelComparison(X.train, Y.train, K=10, type='all', test=testdf, penalty=2,
                      degree=2, gcv=4, covtype='matern3_2', formula=Y~.)
dev.off()

#' The modelComparison function does not allow some important options
#' to the kriging model to be controlled, thus we'll train that model
#' separately and compare it's performance with other models.

cl <- parallel::makeForkCluster()
doParallel::registerDoParallel(cl)

km.train <- modelFit(X=X.train, Y=Y.train, type='Kriging', formula=Y~.,
                     covtype='matern3_2', control=list(maxit=1e3, trace=FALSE),
                     nugget.estim=TRUE, multistart=mc.cores)

Y.test.km <- modelPredict(km.train, X.test)
val.km <- c(R2=R2(Y.test, Y.test.km), RMSE=RMSE(Y.test, Y.test.km))
cbind(mc$Test, val.km)

#' Both kriging models have superior predictive performance on the
#' test set.

mc$CV

#' The untuned kriging model also did better than any other model in
#' terms of the cross validation.

plot(km.train)

#' The mean prediction errors in leave-one-out cross validation are
#' somewhat overdispersed relative to a normal distribution, but no
#' outliers are apparent.

noise.var <- rep(val.km['RMSE']^2, length(Y))

km.m1 <- modelFit(X=X, Y=Y, type='Kriging', formula=Y~.,
                  covtype='matern3_2', control=list(maxit=1e3, trace=FALSE),
                  noise.var=noise.var, multistart=mc.cores)

Y.km.m1 <- modelPredict(km.m1, newdata=X)
Yres.m1 <- (Y - Y.km.m1)^2
km.v1 <- modelFit(X=X, Y=Yres.m1, type='Kriging', formula=Y~1, covtype='matern3_2',
                  control=list(maxit=1e3, trace=TRUE), multistart=mc.cores)

center <- apply(X, 2, median)
pdf('section-views-km.v1.pdf', width=5, height=30)
sectionview.km(model=km.v1$model, center=center, mfrow=c(9, 1))
dev.off()

#contourview.km(model=km.v1$model, center=center, axis=matrix(c(3, 8), nrow=1))

plot(km.v1$model)

#' The leave-one-out errors are by no means gaussian, but still we
#' should be capturing peaks and valleys in the variance.

noise.var.km.v1 <- modelPredict(km.v1, newdata=X)

km.m2 <- modelFit(X=X, Y=Y, type='Kriging', formula=Y~.,
                  covtype='matern3_2', control=list(maxit=1e3, trace=FALSE),
                  noise.var=noise.var.km.v1, multistart=mc.cores)

Y.km.m2 <- modelPredict(km.m2, newdata=X)
Yres.m2 <- (Y - Y.km.m2)^2
km.v2 <- modelFit(X=X, Y=Yres.m2, type='Kriging', formula=Y~1, covtype='matern3_2',
                  control=list(maxit=1e3, trace=TRUE), multistart=mc.cores)

parallel::stopCluster(cl)
save.image('sim-study-checkpoint4.rda')

nmeta <- 1e4
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

RunSobol <- function(nmeta, kmm2, kmv2, all.par.ranges, order=1){
  vars <- km@covariance@var.names
  rngs <- all.par.ranges[vars]
  X1 <- GetRandLHSDes(nmeta, rngs)
  X2 <- GetRandLHSDes(nmeta, rngs)
  colnames(X1) <- colnames(X2) <- vars
  Wrapper <- function(X){
    predict(kmm2, newdata=X, type='UK')$m
  }
  sob <- sensitivity::sobol(model=Wrapper, X1=X1, X2=X2, order=order, nboot=100)
  vym <- sob$V['global', 'original']
  vyd <- DiceKriging::coef(kmv2)$trend
  vy <- vym + vyd
  sob.inds <- sob$S[, 'original'] * vym / vy
  names(sob.inds) <- rownames(sob$S)
  sob.ind.rand <- vyd / vy
  list(sob, sob.inds, sob.ind.rand)
}

(sob.out <- RunSobol(nmeta, kmm2=km.m2$model, kmv2=km.v2$model,
                     all.par.ranges=all.par.ranges, order=1))

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
