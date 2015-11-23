#!/usr/bin/Rscript
library(methods) #for raster

set.seed(124, "L'Ecuyer")

mc.cores <- ifelse(Sys.info()['sysname'] == "Linux",
                   parallel::detectCores() - 1, 1)
mc.cores <- ifelse(mc.cores > 20, 20, mc.cores)
mc.cores <- ifelse(mc.cores == 0, 1, mc.cores)
options('mc.cores'=mc.cores)

load('sim-study-checkpoint3.rda')

GetMetaModels <- function(resall, df){
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

  pdf('modelComparison-plots.pdf')
  mc <- DiceEval::modelComparison(X.train, Y.train, K=10, type='all', test=testdf, penalty=2,
                                  degree=2, gcv=4, covtype='matern3_2', formula=Y~.)
  dev.off()

  #' The modelComparison function does not allow some important options
  #' to the kriging model to be controlled, thus we'll train that model
  #' separately and compare it's performance with other models.

  cl <- parallel::makeForkCluster()
  on.exit(parallel::stopCluster(cl))
  doParallel::registerDoParallel(cl)

  km.train <- DiceEval::modelFit(X=X.train, Y=Y.train, type='Kriging', formula=Y~.,
                                 covtype='matern3_2', control=list(maxit=1e3, trace=FALSE),
                                 nugget.estim=TRUE, multistart=max(mc.cores, 5))

  Y.test.km <- DiceEval::modelPredict(km.train, X.test)
  val.km <- c(R2=DiceEval::R2(Y.test, Y.test.km), RMSE=DiceEval::RMSE(Y.test, Y.test.km))
  print(cbind(mc$Test, val.km))

  #' Both kriging models have superior predictive performance on the
  #' test set.

  print(mc$CV)

  #' The untuned kriging model also did better than any other model in
  #' terms of the cross validation.

  DiceKriging::plot(km.train$model)

  #' The mean prediction errors in leave-one-out cross validation are
  #' somewhat overdispersed relative to a normal distribution, but no
  #' outliers are apparent.

  noise.var <- rep(val.km['RMSE']^2, length(Y))

  km.m1 <- DiceEval::modelFit(X=X, Y=Y, type='Kriging', formula=Y~.,
                              covtype='matern3_2', control=list(maxit=1e3, trace=FALSE),
                              noise.var=noise.var, multistart=max(mc.cores, 5))

  Y.km.m1 <- DiceEval::modelPredict(km.m1, newdata=X)
  Yres.m1 <- (Y - Y.km.m1)^2
  km.v1 <- DiceEval::modelFit(X=X, Y=Yres.m1, type='Kriging', formula=Y~1, covtype='matern3_2',
                              control=list(maxit=1e3, trace=TRUE), multistart=max(mc.cores, 5))

  center <- apply(X, 2, median)
  pdf('section-views-km.v1.pdf', width=5, height=30)
  DiceView::sectionview.km(model=km.v1$model, center=center, mfrow=c(9, 1))
  dev.off()

  #contourview.km(model=km.v1$model, center=center, axis=matrix(c(3, 8), nrow=1))

  DiceKriging::plot(km.v1$model)

  #' The leave-one-out errors are by no means gaussian, but still we
  #' should be capturing peaks and valleys in the variance.

  noise.var.km.v1 <- DiceEval::modelPredict(km.v1, newdata=X)

  km.m2 <- DiceEval::modelFit(X=X, Y=Y, type='Kriging', formula=Y~.,
                              covtype='matern3_2', control=list(maxit=1e3, trace=FALSE),
                              noise.var=noise.var.km.v1, multistart=max(mc.cores, 5))

  Y.km.m2 <- DiceEval::modelPredict(km.m2, newdata=X)
  Yres.m2 <- (Y - Y.km.m2)^2
  km.v2 <- DiceEval::modelFit(X=X, Y=Yres.m2, type='Kriging', formula=Y~1, covtype='matern3_2',
                              control=list(maxit=1e3, trace=TRUE), multistart=max(mc.cores, 5))

  list(m1=km.m1, m2=km.m2, v1=km.v1, v2=km.v2, comparisons=mc, center=center)
}

kms <- GetMetaModels(resall, df)
save.image('sim-study-checkpoint4.rda')
