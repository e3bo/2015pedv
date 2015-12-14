
#' @export
GetCores <- function(){
  mc.cores <- ifelse(Sys.info()['sysname'] == "Linux", parallel::detectCores() - 1, 1)
  mc.cores <- ifelse(mc.cores > 20, 20, mc.cores)
  ifelse(mc.cores == 0, 1, mc.cores)
}

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

#' @export
RunSobol <- function(nmeta, kmm2, kmv2, all.par.ranges, order=1){
  vars <- kmm2@covariance@var.names
  rngs <- all.par.ranges[vars]
  X1 <- GetRandLHSDes(nmeta, rngs)
  X2 <- GetRandLHSDes(nmeta, rngs)
  colnames(X1) <- colnames(X2) <- vars
  sob <- sensitivity::sobol(model=NULL, X1=X1, X2=X2, order=order, nboot=1000)
  X <- ff::as.ff(sob$X)
  sob$X <- NULL        ## Conserve memory here, X can easily be regenerated from X1 and X2
  max.chunksize <- 1e3 ## approximate max, based on complexity of predict.km and RAM available
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
    predict(kmm2, newdata=X, type='UK', se.compute=FALSE, light.return=TRUE)$m
  }
  y <- unlist(parallel::mclapply(chunks, Wrapper))
  sensitivity::tell(sob, y=y)

  vym <- sob$V['global', 'original']
  vyd <- DiceKriging::coef(kmv2)$trend
  vy <- vym + vyd
  sob.inds <- sob$S[, c("original", "min. c.i.", "max. c.i.")] * vym / vy
  rownames(sob.inds) <- rownames(sob$S)
  sob.ind.rand <- vyd / vy
  list(sob, sob.inds, sob.ind.rand)
}

#' @export
GetMetaModels <- function(resall, df, covtype='matern3_2', var1='lag1',
                          var2='shipment', cortype='spearman',
                          is.symmetric='TRUE', statistic='r', make.plots=FALSE){
  sub <- resall
  Y <- sub[ , paste('mantel', statistic, var1, var2, cortype, is.symmetric,
                        df$permutations[1], sep='.')]

  GetVaryingInputs <- function(df){
    rngs <- apply(df, 2, range)
    is.varying <- apply(rngs, 2, diff) > sqrt(.Machine$double.eps)
    ret <- names(is.varying[is.varying])
    is.ind <- grepl('\\.ind$', ret)
    ret[!is.ind]
  }

  X <- sub[, GetVaryingInputs(df)]

  ntest <- 0.5 * nrow(X)
  test.ind <- seq_len(ntest)
  X.test <- X[test.ind, ]
  Y.test <- Y[test.ind]
  X.train <- X[-test.ind, ]
  Y.train <- Y[-test.ind]
  testdf <- cbind(X.test, Y.test)

  file <- paste0('modelComparison-plots-', statistic, '-', var2, '.pdf')
  pdf(file)
  mc <- DiceEval::modelComparison(X.train, Y.train,
                                  type=c('StepLinear', 'Additive', 'PolyMARS', 'MARS'),
                                  test=testdf, penalty=2, degree=2,
                                  gcv=4, formula=Y~.)
  dev.off()

  # The modelComparison function does not allow some important options
  # to the kriging model to be controlled, thus we'll train that model
  # separately and compare it's performance with other models.

  GetKm <- function(X, Y, ...){
    control <- list(max.generations=100, control=list(maxit=1e3))
    DiceEval::modelFit(X=X, Y=Y, type='Kriging', ..., covtype=covtype,
                       control=control, optim.method='gen')
  }
  km.train <- GetKm(X=X.train, Y=Y.train, formula=~.,
                    nugget.estim=TRUE, nugget=1e-7)

  Y.test.km <- DiceEval::modelPredict(km.train, X.test)
  val.km <- c(R2=DiceEval::R2(Y.test, Y.test.km), RMSE=DiceEval::RMSE(Y.test, Y.test.km))
  print(cbind(mc$Test, val.km))

  print(mc$CV)

  km.m1 <- GetKm(X=X, Y=Y, formula=~., nugget.estim=TRUE, nugget=1e-7)

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
  km.v1 <- GetKm(X=X, Y=Yres2.m1, formula=~1, nugget.estim=TRUE)

  Y.km.v1 <- GetPredNuggetAsNoise(km.v1$model)
  Yres.v1 <- Yres2.m1 - Y.km.v1
  noise.var <- ifelse(Y.km.v1 < 0, 0, Y.km.v1)

  km.m2 <- GetKm(X=X, Y=Y, formula=~., noise.var=noise.var)
  Y.km.m2 <- DiceEval::modelPredict(km.m2, newdata=X)
  Yres2.m2 <- (Y - Y.km.m2)^2
  km.v2 <- GetKm(X=X, Y=Yres2.m2, formula=~1, nugget.estim=TRUE)

  center <- apply(X, 2, median)
  if(make.plots){
    file <- paste0('section-views-km-', statistic, '-', var2, '.v1.pdf')
    pdf(file, width=5, height=30)
    DiceView::sectionview.km(model=km.v1$model, center=center, mfrow=c(length(center), 1))
    dev.off()

    file <- paste0('km-training-model-diagnostics-', statistic, '-', var2, '.pdf')
    pdf(file)
    DiceKriging::plot(km.train$model)
    dev.off()

    file <- paste0('noise-var-distribution-', statistic, '-', var2, '.pdf')
    pdf(file)
    par(mfrow=c(3, 1))
    hist(Y.km.v1)
    hist(noise.var)
    qqnorm(noise.var)
    dev.off()

    file <- paste0('km-v1-model-residuals-', statistic, '-', var2, '.pdf')
    pdf(file)
    par(mfrow=c(3,1))
    plot(Y.km.v1~Yres2.m1)
    plot(Yres.v1/sd(Yres.v1))
    qqnorm(Yres.v1)
    dev.off()
  }

  list(m1=km.m1, m2=km.m2, v1=km.v1, v2=km.v2, comparisons=mc, center=center)
}
