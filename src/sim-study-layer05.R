#!/usr/bin/Rscript
library(methods) #for raster

set.seed(125, "L'Ecuyer")

mc.cores <- ifelse(Sys.info()['sysname'] == "Linux",
                   parallel::detectCores() - 1, 1)
mc.cores <- ifelse(mc.cores > 20, 20, mc.cores)
mc.cores <- ifelse(mc.cores == 0, 1, mc.cores)
options('mc.cores'=mc.cores)

load('sim-study-checkpoint4.rda')

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

RunSobol <- function(nmeta, kmm2, kmv2, all.par.ranges, order=1){
  vars <- kmm2@covariance@var.names
  rngs <- all.par.ranges[vars]
  X1 <- GetRandLHSDes(nmeta, rngs)
  X2 <- GetRandLHSDes(nmeta, rngs)
  colnames(X1) <- colnames(X2) <- vars
  sob <- sensitivity::sobol(model=NULL, X1=X1, X2=X2, order=order, nboot=100)
  X <- sob$X
  max.chunksize <- 1e3 ## approximate max, based on limitation in predict.km and RAM available
  chunksize <- min(ceiling(nrow(X) / getOption('mc.cores')), max.chunksize)
  nchunks <- ceiling(nrow(X) / chunksize)

  breaks <- chunksize * seq(0, nchunks)
  breaks[length(breaks)] <- nrow(X)
  chunks <- list()
  for(i in seq(2, nchunks + 1)){
    l <- breaks[i - 1] + 1
    u <- breaks[i]
    chunks[[i - 1]] <- X[l:u, ]
  }
  Wrapper <- function(X){
    predict(kmm2, newdata=X, type='UK')$m
  }
  y <- unlist(parallel::mclapply(chunks, Wrapper))
  sensitivity::tell(sob, y=y)

  vym <- sob$V['global', 'original']
  vyd <- DiceKriging::coef(kmv2)$trend
  vy <- vym + vyd
  sob.inds <- sob$S[, 'original'] * vym / vy
  names(sob.inds) <- rownames(sob$S)
  sob.ind.rand <- vyd / vy
  list(sob, sob.inds, sob.ind.rand)
}

sob.out <- RunSobol(nmeta, kmm2=kms$m2$model, kmv2=kms$v2$model,
                     all.par.ranges=all.par.ranges, order=1)
save.image('sim-study-checkpoint5.rda')

sob.out
tp.ind <- which(names(kms$center) == 'tprob.net')
sa.ind <- which(names(kms$center) == 'seasonal.amplitude')

pdf('km-m2-views.pdf')
DiceView::sectionview(kms$m2, axis=tp.ind, center=kms$center, mfrow=c(1,1))
DiceView::contourview(kms$m2, axis=matrix(c(sa.ind, tp.ind), nrow=1), center=kms$center, mfrow=c(1, 1))
dev.off()
