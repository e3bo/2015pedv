#!/usr/bin/Rscript
library(methods) #for raster

set.seed(124, "L'Ecuyer")

mc.cores <- ifelse(Sys.info()['sysname'] == "Linux",
                   parallel::detectCores() - 1, 1)
mc.cores <- ifelse(mc.cores > 20, 20, mc.cores)
mc.cores <- ifelse(mc.cores == 0, 1, mc.cores)
options('mc.cores'=mc.cores)

load('sim-study-checkpoint3.rda')

GetMetaModels <- function(resall, df, covtype='matern3_2'){
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
  mc <- DiceEval::modelComparison(X.train, Y.train, type='all', test=testdf, penalty=2,
                                  degree=2, gcv=4, covtype=covtype, formula=Y~.,
                                  nugget.estim=TRUE)
  dev.off()

  #' The modelComparison function does not allow some important options
  #' to the kriging model to be controlled, thus we'll train that model
  #' separately and compare it's performance with other models.

  if (getOption('mc.cores') > 1){
    cl <- parallel::makeForkCluster()
    on.exit(parallel::stopCluster(cl))
    doParallel::registerDoParallel(cl)
  }

  km.train <- DiceEval::modelFit(X=X.train, Y=Y.train, type='Kriging', formula=Y~.,
                                 covtype=covtype, control=list(maxit=1e3, trace=FALSE),
                                 nugget.estim=TRUE, multistart=max(mc.cores, 5))

  Y.test.km <- DiceEval::modelPredict(km.train, X.test)
  val.km <- c(R2=DiceEval::R2(Y.test, Y.test.km), RMSE=DiceEval::RMSE(Y.test, Y.test.km))
  print(cbind(mc$Test, val.km))

  #' Both kriging models have superior predictive performance on the
  #' test set.

  print(mc$CV)

  #' The untuned kriging model also did better than any other model in
  #' terms of the cross validation.
  pdf('km-training-model-diagnostics.pdf')
  DiceKriging::plot(km.train$model)
  dev.off()

  #' The mean prediction errors in leave-one-out cross validation are
  #' somewhat overdispersed relative to a normal distribution, but no
  #' outliers are apparent.

  km.m1 <- DiceEval::modelFit(X=X, Y=Y, type='Kriging', formula=~.,
                              covtype=covtype, control=list(maxit=1e3, trace=TRUE),
                              nugget.estim=TRUE, multistart=max(mc.cores, 5),
                              nugget=1e-7)

  GetPredNuggetAsNoise <- function(mod){
    noise.var <- rep(mod@covariance@nugget, len=nrow(mod@X))
    mpred <- DiceKriging::km(mod@trend.formula, design=mod@X, response=mod@y,
                             covtype=mod@covariance@name,
                             coef.cov=DiceKriging::covparam2vect(mod@covariance),
                             coef.trend=mod@trend.coef,
                             coef.var=mod@covariance@sd2,
                             noise.var=noise.var)
    predict(mpred, newdata=mod@X, type='UK', se.fit=FALSE, light.return=TRUE)$mean
  }
  Y.km.m1 <- GetPredNuggetAsNoise(km.m1$model)
  Yres2.m1 <- (Y - Y.km.m1)^2
  km.v1 <- DiceEval::modelFit(X=X, Y=Yres2.m1, type='Kriging', formula=Y~1, covtype=covtype,
                              control=list(maxit=1e3, trace=TRUE), multistart=max(mc.cores, 5),
                              nugget.estim=TRUE)

  center <- apply(X, 2, median)
  pdf('section-views-km.v1.pdf', width=5, height=30)
  DiceView::sectionview.km(model=km.v1$model, center=center, mfrow=c(9, 1))
  dev.off()

  Y.km.v1 <- GetPredNuggetAsNoise(km.v1$model)
  Yres.v1 <- Yres2.m1 - Y.km.v1

  pdf('km-v1-model-residuals.pdf')
  par(mfrow=c(3,1))
  plot(Y.km.v1~Yres2.m1)
  plot(Yres.v1/sd(Yres.v1))
  qqnorm(Yres.v1)
  dev.off()

  noise.var <- ifelse(Y.km.v1 < 0, 0, Y.km.v1)
  pdf('noise-var-distribution.pdf')
  par(mfrow=c(3, 1))
  hist(Y.km.v1)
  hist(noise.var)
  qqnorm(noise.var)
  dev.off()

  km.m2 <- DiceEval::modelFit(X=X, Y=Y, type='Kriging', formula=Y~.,
                              covtype=covtype, control=list(maxit=1e3, trace=FALSE),
                              noise.var=noise.var, multistart=max(mc.cores, 5))

  Y.km.m2 <- DiceEval::modelPredict(km.m2, newdata=X)
  Yres2.m2 <- (Y - Y.km.m2)^2
  km.v2 <- DiceEval::modelFit(X=X, Y=Yres2.m2, type='Kriging', formula=Y~1, covtype=covtype,
                              control=list(maxit=1e3, trace=TRUE), nugget.estim=TRUE,
                              multistart=max(mc.cores, 5))

  list(m1=km.m1, m2=km.m2, v1=km.v1, v2=km.v2, comparisons=mc, center=center)
}

kms <- GetMetaModels(resall, df)
save.image('sim-study-checkpoint4.rda')
saveRDS(kms$m2$model, "kmm2.rds")
saveRDS(kms$v2$model, "kmv2.rds")
saveRDS(kms$center, "center.rds")

tp.ind <- which(names(kms$center) == 'tprob.net')
sa.ind <- which(names(kms$center) == 'seasonal.amplitude')
pdf('km-m2-views.pdf')
DiceView::sectionview(kms$m2$model, axis=tp.ind, center=kms$center, mfrow=c(1,1))
DiceView::sectionview(kms$m2, axis=tp.ind, center=kms$center, mfrow=c(1,1))
DiceView::contourview(kms$m2, axis=matrix(c(sa.ind, tp.ind), nrow=1), center=kms$center, mfrow=c(1, 1))
dev.off()
