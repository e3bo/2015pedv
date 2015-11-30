#+setup, include=FALSE, cache=FALSE
library(knitr)
opts_chunk$set(fig.path='flows-figures/', fig.align='center', fig.show='hold')
options(replace.assign=TRUE,width=80)

#' ## Load Packages

set.seed(4253)
library(lme4)
library(car)
library(Hmisc)
library(pscl)
library(glmmADMB)
library(ggplot2)
print(sessionInfo())

#' ## Data prep
iqr <- function(x) diff(quantile(x, probs=c(0.25, 0.75), na.rm=TRUE))
scaleiqr <- function(x) scale(x, scale=iqr(x))
clog <- function(x) log(x + 0.5)

ObservedMatchingInternalFlows <- function(flows, case.data) {
    unwanted <- c('week', 'totalNumberSwineAccessions', 'Unk')
    ind <- which(!colnames(case.data) %in% unwanted)
    observed <- case.data[, ind]
    avail <- rownames(flows)[!is.na(flows[, 'impInternalFlow'])]
    zeros <- setdiff(avail, colnames(observed))
    observed[, zeros] <- 0
    observed
}

ReshapeAugmentObserved <- function(observed, fiAll, flows, stateCty) {
  om <- reshape(observed, v.names='cases', varying=colnames(observed),
                direction='long', timevar='state', times=colnames(observed),
                idvar='week')
  splt <- split(om, om$state)
  tmpf <- function(x) {
      ord <- order(x$week)
      x <- x[ord, ]
      x$cumCases <- cumsum(x$cases)
      x
  }
  proc <- lapply(splt, tmpf)
  om <- do.call(rbind, args=proc)
  om$nFarms <- fiAll[om$state, ]$X2002.Farms
  om$nSusc <- om$nFarms - om$cumCases
  omprev <- om
  omprev$week <- omprev$week + 1
  omprev$nFarms <- NULL
  omprev$cumCases <- NULL
  c1 <- which(colnames(omprev) == 'cases')
  c2 <- which(colnames(omprev) == 'nSusc')
  colnames(omprev)[c1] <- 'cases1WksAgo'
  colnames(omprev)[c2] <- 'nSusc1WksAgo'
  om <- merge(om, omprev, by=c('week', 'state'), all.x=TRUE)
  omprev$week <- omprev$week + 1
  colnames(omprev)[c1] <- 'cases2WksAgo'
  colnames(omprev)[c2] <- 'nSusc2WksAgo'
  om <- merge(om, omprev, by=c('week', 'state'), all.x=TRUE)
  omprev$week <- omprev$week + 1
  colnames(omprev)[c1] <- 'cases3WksAgo'
  colnames(omprev)[c2] <- 'nSusc3WksAgo'
  om <- merge(om, omprev, by=c('week', 'state'), all.x=TRUE)
  omprev$week <- omprev$week + 1
  colnames(omprev)[c1] <- 'cases4WksAgo'
  colnames(omprev)[c2] <- 'nSusc4WksAgo'
  om <- merge(om, omprev, by=c('week', 'state'), all.x=TRUE)
  om$internalFlow <- flows[om$state, "impInternalFlow"]
  rownames(om) <- paste(om$week, om$state, sep='.')
  test <- !is.na(om$internalFlow)
  om <- om[test, ]
  om <- merge(om, stateCty, by='state')
  rownames(om) <- paste(om$week, om$state, sep='.')
  om$logInternalFlowScaled <- scaleiqr(log(2 * om$internalFlow))
  om$weekCent <- scale(om$week, center=TRUE, scale=FALSE)
  om$logCmedDenseScaled <- scaleiqr(log(om$cmedDense * om$nFarms * om$nFarms))
  om$logMDenseScaled <- scaleiqr(log(om$mDense * om$nFarms * om$nFarms))
  om$clogCases1wa <- clog(om$cases1WksAgo)
  om$clogCases2wa <- clog(om$cases2WksAgo)
  om$clogCases3wa <- clog(om$cases3WksAgo)
  om$clogCases4wa <- clog(om$cases4WksAgo)
  om$clogCases1wa05 <- log(om$cases1WksAgo + 0.05)
  om$clogCases1wa005 <- log(om$cases1WksAgo + 0.005)
  om
}

data('real.case.data', package='sds')
data('internal.flows', package='sds')
data('farms.and.inventory', package='sds')
data('state.summaries', package='sds')

observed <- ObservedMatchingInternalFlows(flows=internal.flows, case.data=real.case.data)
om <- ReshapeAugmentObserved(observed=observed, fiAll=farms.and.inventory,
                             flows=internal.flows, stateCty=state.summaries)

#' ## Model fitting
#'
#' Let's start with a least square fit on the logged response. We'll
#' center the predictors around means to reduce the correlation
#' between slopes and intercepts.

m <- list()
f <- list()
f$lm <- as.formula(clog(cases) ~ clogCases1wa + logInternalFlowScaled + logCmedDenseScaled + weekCent + offset(log(nSusc1WksAgo) - 2*log(nFarms)))
m$lm <- lm(f$lm, data=om)
summary(m$lm)

#' We've included week to account for the general trend of increasing
#' cases over the period studied, possibly due to seasonal changes in
#' manure spreading, virus survivablility, or pig production.

tmpf <- function() {
  allCases <- with(om, tapply(cases, week, sum))
  week <- as.integer(names(allCases))
  plot(week, allCases, log='y')
}
tmpf()

#' Next we use a random effect to account for lack of dependence within states.
f$lme <- update(f$lm, . ~ . + (1|state))
f$lmeN <- update(f$lme,  . ~ . - logInternalFlowScaled)
m$lme <- lmer(f$lme, data=om, REML=FALSE)
m$lmeN <- lmer(f$lmeN, data=om, REML=FALSE)
m$lme

likRatio <- function(null, alt){
    null <- logLik(null)
    alt <- logLik(alt)
    dfn <- attr(null, "df")
    dfa <- attr(alt, "df")
    if(dfn >= dfa) stop("null model has greater or same df than alternative")
    df <- dfa - dfn
    likRatio <- -2 * (as.numeric(null) - as.numeric(alt))
    p <- pchisq(likRatio, df=df, lower.tail = FALSE)/2
    res <- c(lr=likRatio, df=df, p=p)
    signif(res, 3)
}
likRatio(m$lmeN, m$lme)
#' The effect of flows is still statistically significant.

tmpf <- function() {
    df <- model.frame(m$lme)
    iqry <- quantile(df$"clog(cases)", probs=c(.1, .9))
    effectSize <- fixef(m$lme)['logInternalFlowScaled']
    effectSize/diff(iqry)
}
tmpf()
#' The effect size is non-negligible.
#'
#' A count model would be more appropriate for the discrete outcome.

f$ct <- update(f$lm, cases ~ .)
m$pois <- glm(f$ct, family=poisson, data=om)
getDispStat <- function(mod) {
    df <- mod$df.residual
    E <- residuals(mod, type = "pearson")
    sum(E^2)/df
}
getDispStat(m$pois)

#' The dispersion parameter is fairly large. Let's check for outliers.

(otest <- outlierTest(m$pois))

#' There's a fairly large number, given the size of the data set,
#' which suggests that a negative binomial model may be more
#' appropriate. Such a model would also have a theoretical
#' justification, as in the TSIR model.

m$nb <- glm.nb(f$ct, data=om)
summary(m$nb)
likRatio(m$pois, m$nb)

#' The likelihod ratio test supports the negative binomial model.
#' Is this support sensitive to outliers?

tmpf <- function() {
    orows <- names(otest$rstudent)
    test <- !rownames(om) %in% orows
    omno <- om[test, ]
    mpno <- glm(f$ct, family='poisson', data=omno)
    print(getDispStat(mpno))
    print(outlierTest(mpno))
    mnbno <- glm.nb(f$ct, data=omno)
    print(summary(mnbno))
    print(likRatio(mpno, mnbno))
}
tmpf()
rm(otest)
#' No.

#' Does the negative binomial model support the includion of the flow term?

m$nbN <- update(m$nb, . ~ . - logInternalFlowScaled)
likRatio(m$nbN, m$nb)

#' Yes. What if we use a different predictor for geographic density?
#' The current choice is highly correlated with flows on the log scale
#' as shown next, but we'll try another one to be safe.

tmpf <- function(stateCty, fiAll, flows) {
    sel <- unique(om$state)
    rownames(stateCty) <- stateCty$state
    den <- stateCty[sel, c('mDense', 'cmDense', 'medDense', 'cmedDense')]
    nf <- fiAll[sel,]$X2002.Farms
    fl <- flows[sel, 'impInternalFlow']
    tmpff <- function(x) {
        list(flowVdense=cor(log(cbind(flow=fl, den)), method=x),
             flowVdenseNf=cor(log(cbind(flow=fl, den*nf)), method=x),
             flowVdesneNfNf=cor(log(cbind(flow=fl, den*nf*nf)), method=x))
    }
    lapply(c(pearson='pearson', spearman='spearman'), tmpff)
}
tmpf(stateCty=state.summaries, fiAll=farms.and.inventory, flows=internal.flows)

m$nbmDense <- update(m$nb, . ~ . - logMedDenseScaled + logMDenseScaled)
m$nbmDenseN <- update(m$nbmDense, . ~ . - logInternalFlowScaled, control=glm.control(maxit=60))
summary(m$nbmDense)
likRatio(m$nbmDenseN, m$nbmDense)

#' There coefficients and significance level are similar.

#' Let's test different lags.

f$ct2 <- as.formula(cases ~ clogCases2wa + logCmedDenseScaled + logInternalFlowScaled + weekCent + offset(log(nSusc2WksAgo) - 2 * log(nFarms)))
m$nb2 <- glm.nb(f$ct2, data=om)
logLik(m$nb2)
logLik(glm.nb(f$ct, data=om[om$week > 2, ]))

f$ct3 <- as.formula(cases ~ clogCases3wa + logCmedDenseScaled + logInternalFlowScaled + weekCent + offset(log(nSusc3WksAgo) - 2 * log(nFarms)))
m$nb3 <- glm.nb(f$ct3, data=om, control=glm.control(maxit=50))
logLik(m$nb3)
logLik(glm.nb(f$ct, data=om[om$week > 3, ]))

f$ct4 <- as.formula(cases ~ clogCases4wa + logCmedDenseScaled + logInternalFlowScaled + weekCent + offset(log(nSusc4WksAgo) - 2 * log(nFarms)))
m$nb4 <- glm.nb(f$ct4, data=om, control=glm.control(maxit=50))
logLik(m$nb4)
logLik(glm.nb(f$ct, data=om[om$week > 4, ]))

#' A lag of 1 is by far the best by likelihood.

#' Let's try different baseline hazards.

f$nb05 <- update(f$ct, . ~ . - clogCases1wa + clogCases1wa05)
m$nb05 <- glm.nb(f$nb05, data=om, control=glm.control(maxit=100))

f$nb005 <- update(f$ct, . ~ . - clogCases1wa + clogCases1wa005)
m$nb005 <- glm.nb(f$nb005, data=om, control=glm.control(maxit=100))

lapply(m[c('nb', 'nb05', 'nb005')], logLik)

# A baseline hazard of 0.05 is best, although the model seems fairly
# insensitive to this parameter.

#' Let's compare the observed and expected distribution of cases after
#' summing over the variation in the mean.

makeCompPlot <- function(xvals = 0:10, mod, exclude = integer(0), ...) {
    ypred <- predict(mod, type = "response")
    tmpf <- function(x) {
        ret <- dnbinom(xvals, mu = x, size=mod$theta)
        ret/sum(ret)
    }
    ds <- sapply(ypred, tmpf)
    expectedProbs <- rowSums(ds)/ncol(ds)
    resp <- mod$model$cases
    if (length(exclude) > 0) {
        resp <- resp[-exclude]
    }
    obsResp <- table(resp)
    expResp <- expectedProbs * sum(obsResp)
    ye <- expResp
    yo <- obsResp[as.character(xvals)]
    plot(xvals, yo, pch = "o", xlab = "response value", ylab = "number", ...,
         type = "b")
    points(xvals, ye, pch = "e", type = "b")
}
makeCompPlot(mod=m$nb, xvals=0:30, log='y')

#' There's no strong signs of zero inflation, or any other major
#' defects.

m$nbzi <- zeroinfl(f$ct, data=om, dist="negbin")
AIC(m$nb, m$nbzi)

#' One could justify a zero-inflated model based on AIC, but it's not
#' clear to me such a model is better here. What mechanism would cause
#' the zero-inflation, and do we have a good reason to spend degrees
#' of freedom on it?
#'
#' Let's check for evidence of non-linearity in partial residual plots.

crPlots(m$nb)

#' Toward the extreme ends of the data range, there is some sign of
#' non-linearity, but generally linearity assumptions seem reasonable.

#' Next let's check for influential observations.

influencePlot(m$nb)

tmpf <- function() {
    filterInf <- function(mod, df, cutoff=0.2) {
        dfb <- dfbetas(mod)
        test <- rowSums(dfb > cutoff) > 0
        if(any(test)) {
            inf <- dfb[test, ]
            print(inf)
            infNames <- names(which(test))
             df[!rownames(df) %in% infNames, ]
        } else {
            NULL
        }
    }
    i <- 1
    mods <- list(m$nb)
    omni <- filterInf(m$nb, om)
    while(!is.null(omni)) {
        omnip <- omni
        mnbni <- glm.nb(f$ct, data=omnip)
        omni <- filterInf(mnbni, omni)
        i <- i + 1
        mods[[i]] <- mnbni
    }
    print(sapply(mods, coef))
    mnbniNull <- update(mnbni, . ~ . - logInternalFlowScaled, control=glm.control(maxit=100))
    print(likRatio(mnbniNull, mnbni))
}
tmpf()

#' Removing the influential observations changes the coefficients a
#' bit, but the significance of the flows term is not affected.
#'
#' Next, let's account for correlations of data from the same state by
#' using a clustered bootstrap variance estimate.

clusteredBoot <- function(data, statistic, B, cluster){
    ## based on rms::bootcov() and boot::boot()
    n <- nrow(data)
    if (length(cluster) != n) stop("length of cluster does not match # rows used in fit")
    if (any(is.na(cluster))) stop("cluster contains NAs")
    cluster <- as.character(cluster)
    clusters <- unique(cluster)
    nc <- length(clusters)
    Obsno <- split(1:n, cluster)
    pb <- txtProgressBar(min=1,max=B)
    drawStat <- function(x) {
        setTxtProgressBar(pb, x)
        j <- sample(clusters, nc, replace=TRUE)
        obs <- unlist(Obsno[j])
        statistic(data, obs)
    }
    res <- sapply(1:B, drawStat, simplify='array')
    close(pb)
    res
}

tmpf <- function(B) {
    fun <- function(data, i) {
        m <- glm.nb(f$ct, data=data[i, ], control=glm.control(maxit=200))
        coef(summary(m))[, 1:2]
    }
    res <- clusteredBoot(om, statistic=fun, B=B, cluster=om$state)
    ests <- res['logInternalFlowScaled', 'Estimate', ]
    print(quantile(ests, probs=c(0.025, 0.975)))
    print(quantile(ests, prob=1:5/100))
    qqPlot(ests)
    print(splom(t(res[, 'Estimate', ])))
    res
}
system.time(clusterBootDf <- tmpf(B=300))
save.image(file='flows-checkpoint1.RData')
date()

#' It seems like a good idea to also estimate a parameter for a
#' non-linear effect of the number of susceptibles, since we do that
#' for the number of infectives.

tmpf <- function() {
    form <- as.formula(cases ~ clogCases1wa + logInternalFlowScaled + weekCent + log(nSusc1WksAgo) + offset(-2 * log(nFarms)))
    func <- function(data, i) {
        m <- glm.nb(form, data=data[i, ], control=glm.control(maxit=100))
        coef(summary(m))[, 1:2]
    }
    res2 <- clusteredBoot(om, statistic=func, B=20, cluster=om$state)
    print(splom(t(res2[, 'Estimate', ])))
    print(cor(t(res2[c('log(nSusc1WksAgo)', '(Intercept)'), 'Estimate', ])))
}
tmpf()

#' However, the estimate turns out to be highly correlated with the
#' intercept, probably because the number of suceptibles varies little
#' in relation to the number of farms. We leave this parameter out to
#' avoid inflation of the variance of our other estimates.

#' Next we a random effects model using the lme4 package to see if our
#' results are robust.

#' Check that we understand the parameterization of glmer.nb
tmpf <- function() {
    x <- runif(400) - 0.5
    z <- gl(n=40, k=10)
    m <- model.matrix(~x + z)
    u <- rnorm(40, sd=0.5)
    eta <- m %*% c(0, 3, u[-1] - u[1])
    y <- rnbinom(n=length(eta), size=3, mu=exp(eta))
    data.frame(y, x, z)
}
simdf <- tmpf()
glmer.nb(y~x + (1|z), nAGQ=1L, data=simdf)

#' It would seem that nAGQ>1 should increase accuracy, but rather we
#' find it shrinks the estimate of the size parameter for the negative
#' binomial:

glmer.nb(y~x + (1|z), nAGQ=10L, data=simdf)
rm(simdf)

#' Now with the data:

om$off <- with(om, log(nSusc1WksAgo) - 2 * log(nFarms))
om$offCen <- om$off - mean(om$off, na.rm=TRUE)
# centering the offset prevents a scaling error: https://github.com/lme4/lme4/issues/173
om$stateF <- factor(om$state)

f$glme <- as.formula(cases ~ clogCases1wa + logCmedDenseScaled + logInternalFlowScaled + weekCent + (1|stateF) + offset(offCen))
m$nbme <- glmer.nb(f$glme, data=om, nAGQ=1L)
f$glmeN <- update(f$glme, . ~ . -logInternalFlowScaled)
m$nbmeN <- glmer.nb(f$glmeN, data=om, nAGQ=1L)
likRatio(m$nbmeN, m$nbme)
summary(m$nbme)
summary(m$nbmeN)

#' Now we use glmmADMB for comparison

f$admb <- f$glme
omcc <- om[rownames(model.frame(m$nbme)), ]
m$admb <- glmmadmb(f$admb, data=omcc, family='nbinom',
                   admb.opts=admbControl(shess=FALSE, impSamp=100),
                   extra.args='-gh 10')
m$admb200 <- glmmadmb(f$admb, data=omcc, family='nbinom',
                      admb.opts=admbControl(shess=FALSE, impSamp=200))
m$admb0 <- glmmadmb(f$admb, data=omcc, family='nbinom')
summary(m$admb0)
signif(sapply(list(m$nbme, m$admb, m$admb0, m$admb200), fixef), 3)

#' The results look nearly identical, so importance sampling is
#' probably not needed to ensure accuracy of the integration.

f$admbN <- update(f$admb, . ~ . - logInternalFlowScaled)
m$admbN <- glmmadmb(f$admbN, data=omcc, family='nbinom')
likRatio(m$admbN, m$admb0)

#' glmmADMB also supports the inclusion of the flow terms.
#'
#' Output the descriptive stats for this test

omcc$logInternalFlowScaled <- as.numeric(omcc$logInternalFlowScaled)
omcc$clogCases1wa <- as.numeric(omcc$clogCases1wa)
omcc$logCmedDenseScaled <- as.numeric(omcc$logCmedDenseScaled)
omcc$weekCent <- as.numeric(omcc$weekCent)
desc <- describe(omcc[, c('cases', 'clogCases1wa', 'logCmedDenseScaled',
                          'logInternalFlowScaled', 'weekCent', 'stateF',
                          'offCen')])
invisible(latex(desc, file='tsir-describe.tex'))

save.image('flows-checkpoint2.RData')
date()
rm(omcc)

#' Next we plot the response distribution simplify to get some feeling
#' for our fitted model's behavior.

plotResponseDist <- function(mod, size){
    tmpf <- function(y) {
        x <- 1:40
        dnbinom(x, mu=y, size=size)
    }
    q <- quantile(predict(mod, type='response'), probs=c(.05, .25, .5, .75, .95))
    m <- sapply(q, tmpf)
    main <- paste('Negative Binomial, Theta == ', signif(size, 3))
    matplot(m, type='l', xlab='Cases', ylab='Probability', main=main)
    leg <- signif(q, 2)
    col <- seq_len(length(q))
    lty <- col
    legend(x='topright', legend=leg, lty=lty, title='Mean', col=col)
}
plotResponseDist(m$admb0, summary(m$admb0)$alpha)
plotResponseDist(m$nb, m$nb$theta)

#' Next we check for autocorrelation of the residuals in time.

tmpf <- function(mnb) {
    E <-resid(mnb, type='deviance')
    df <- data.frame(wk=mnb$model$week, st=om[names(E), ]$state  , E=E)
    df <- reshape(df, timevar='st', idvar='wk', direction='wide')
    rownames(df) <- df$wk
    df <- df[order(df$wk), ]
    dft <- as.ts(df[, -1])
    maxInd <- ceiling(ncol(dft) / 4)
    makeInd <- function(x) {
        i <- (x - 1)*4 + 1;
        i:(i + 3)
    }
    ind <- lapply(1:maxInd, makeInd)
    plot0 <- function(x) {
        x <- x[x <= ncol(dft)]
        plot.ts(dft[, x], main='', nc=2, plot.type='multiple', ylab=colnames(dft)[x])
    }
    lapply(ind, plot0)
    ac <- lapply(df[, -1], acf, lag.max=4, plot=FALSE)
    par(mfrow=c(2,2))
    plot1 <- function(x, y) if(!is.na(y)) plot(x, main=y)
    makePanel <- function(x) mapply(plot1, ac[x], names(ac)[x])
    lapply(ind, makePanel)
    par(mfrow=c(1,1))
}
tmpf(m$nb)

#' No major autocorrelations are present for those time series having cases.
#'
#' Now let's see if the data support models with of spread between
#' states. First we try incorporating this by weighting cases in the
#' previous week by the amount of flow between that state and other
#' states.

tmpf <- function(rel='directed', flowMat, flows, fiAll, observed) {
    fmo <- flowMat[state.abb, state.abb]
    diag(fmo) <- flows[state.abb, 'impInternalFlow']
    Fund <- fmo + t(fmo)
    Fdir <- fmo
    ## Fdir[i, j] == head sent to i from j
    if(rel == 'directed') {
        F <- fmo
    } else {
        F <- fmo + t(fmo)
    }
    Fsum <- rowSums(F)
    Fnorm <- F/Fsum
    n <- fiAll[state.abb, ]$X2002.Farms
    W <- Fnorm * n
    W <- t(t(W) / n)
    sel <- colnames(observed)
    W <- W[sel, sel]
    Wobs <- W %*% t(observed)
    list(ts=t(Wobs), Fsum=Fsum)
}
data('flows.matrix', package='sds')
wcDir <- tmpf('directed', flowMat=t(flows.matrix), flows=internal.flows,
              fiAll=farms.and.inventory, observed=observed)
wcUnd <- tmpf('undirected', flowMat=t(flows.matrix), flows=internal.flows,
              fiAll=farms.and.inventory, observed=observed)

par(mfrow=c(3,1))
plot.ts(wcDir$ts, plot.type='single')
plot.ts(wcUnd$ts, plot.type='single')
plot.ts(observed, plot.type='single')
par(mfrow=c(1,1))

tmpf <- function() {
    wc <- wcDir$ts
    wc <- reshape(data.frame(wc), v.names='casesDirectedFlowWeighted', varying=colnames(wc),
                  direction='long', timevar='state', times=colnames(wc), idvar='week')
    wc$week <- wc$week + 1
    colnames(wc)[2] <- 'casesLastWeekDirectedFlowWeighted'
    om <- merge(om, wc, by=c('week', 'state'), all.x=TRUE)

    om$directedFlow <- wcDir$Fsum[om$state]
    om$logDirectedFlowScaled <- scaleiqr(log(om$directedFlow))
    om$clogCases1waDir <- clog(om$casesLastWeekDirectedFlowWeighted)

    wc <- wcUnd$ts
    wc <- reshape(data.frame(wc), v.names='casesDirectedFlowWeighted', varying=colnames(wc),
                  direction='long', timevar='state', times=colnames(wc), idvar='week')
    wc$week <- wc$week + 1
    colnames(wc)[2] <- 'casesLastWeekUndirectedFlowWeighted'
    om <- merge(om, wc, by=c('week', 'state'), all.x=TRUE)

    om$undirectedFlow <- wcUnd$Fsum[om$state]
    om$logUndirectedFlowScaled <- scaleiqr(log(om$undirectedFlow))
    om$clogCases1waUnd <- clog(om$casesLastWeekUndirectedFlowWeighted)
    rownames(om) <- paste(om$week, om$state, sep='.')
    om
}
om <- tmpf()

f$ctdir <- update(f$ct, . ~ . + clogCases1waDir + logDirectedFlowScaled - clogCases1wa - logInternalFlowScaled)
m$nbdir <- glm.nb(f$ctdir, data=om)

f$ctundir <- update(f$ct, . ~ . + clogCases1waUnd + logUndirectedFlowScaled - clogCases1wa - logInternalFlowScaled)
m$nbundir <- glm.nb(f$ctundir, data=om)

m[c('nb', 'nbdir', 'nbundir')]

#' The internal flow model has the lowest AIC among the fixed effects
#' models. Both of the other models have negligible time components,
#' but otherwise the coefficients are similar. The undirected model
#' has the highest AIC.

f$glmedir <- update(f$glme, . ~ . + clogCases1waDir + logDirectedFlowScaled - clogCases1wa - logInternalFlowScaled)
m$nbmedir <- glmer.nb(f$glmedir, data=om, nAGQ=1L)
m$admbdir <- glmmadmb(f$glmedir, data=om[!is.na(om$clogCases1waDir), ], family='nbinom')

sapply(m[c('nbmedir', 'admbdir')], fixef)

#' ADMB and lme4 again give similar fits.

f$glmeundir <- update(f$glme, . ~ . + clogCases1waUnd + logUndirectedFlowScaled - clogCases1wa - logInternalFlowScaled)
m$nbmeundir <- glmer.nb(f$glmeundir, data=om, nAGQ=1L)
m$admbundir <- glmmadmb(f$glmeundir, data=om[!is.na(om$clogCases1waUnd), ], family='nbinom')

sapply(m[c('nbmeundir', 'admbundir')], fixef)

#' Another good match.

sort(sapply(m[c('nbme', 'nbmedir', 'nbmeundir')], AIC))

#' The undirected flow matrix has the lowest AIC.
#'
#' One shortcoming of the above comparison is that we are not
#' optimizing the baseline hazard, and so the result may be sensitive
#' to our initial and somewhat arbitrary setting. So next we optimize
#' this parameter.

fint <- function(x,data) {
    data$logInf <- log2(data$cases1WksAgo + x)
    fm <- as.formula(cases ~ logInf + logInternalFlowScaled + weekCent + (1 | stateF) + offset(offCen))
    fit <- glmmadmb(fm, data=data[!is.na(data$logInf), ], family='nbinom')
    ll <- logLik(glmer.nb(fit, data=data))
    print(paste(signif(c(x, ll), 5), collapse=","))
    ll
}
(hint <- optimize(f=fint, interval=c(0.01, 5), data=om, tol=0.1, maximum=TRUE))

fund <- function(x,data) {
    data$logInf <- log(data$casesLastWeekUndirectedFlowWeighted + x)
    fm <- as.formula(cases ~ logInf + logUndirectedFlowScaled + weekCent + (1 | stateF) + offset(offCen))
    fit <- glmmadmb(fm, data=data[!is.na(data$logInf), ], family='nbinom')
    ll <- as.numeric(logLik(fit))
    print(paste(signif(c(x, ll), 5), collapse=","))
    ll
}
(hund <- optimize(f=fund, interval=c(0.01, 5), data=om, tol=0.1, maximum=TRUE))

fdir <- function(x,data) {
    data$logInf <- log(data$casesLastWeekDirectedFlowWeighted + x)
    fm <- as.formula(cases ~ logInf + logDirectedFlowScaled + weekCent + (1 | stateF) + offset(offCen))
    fit <- glmmadmb(fm, data=data[!is.na(data$logInf), ], family='nbinom')
    ll <- as.numeric(logLik(fit))
    print(paste(signif(c(x, ll), 5), collapse=","))
    ll
}
(hdir <- optimize(f=fdir, interval=c(0.01, 5), data=om, tol=0.1, maximum=TRUE))

om$logInfInt <- log(om$cases1WksAgo + hint$maximum)
om$logInfUnd <- log(om$casesLastWeekUndirectedFlowWeighted + hund$maximum)
om$logInfDir <- log(om$casesLastWeekDirectedFlowWeighted + hdir$maximum)

f$glmeundirhmax <- update(f$glmeundir, . ~ . - clogCases1waUnd + logInfUnd)
m$glmeundirhmax <- glmer.nb(f$glmeundirhmax, data=om)

f$glmehmax <- update(f$glme, . ~ . - clogCases1wa + logInfInt)
m$glmehmax <- glmer.nb(f$glmehmax, data=om)

f$glmedirhmax <- update(f$glmedir, . ~ . - clogCases1waDir + logInfDir)
m$glmedirhmax <- glmer.nb(f$glmedirhmax, data=om)

sapply(m[c('glmehmax', 'glmeundirhmax', 'glmedirhmax')], logLik)
abs(unname(diff(sapply(m[c('glmehmax', 'glmeundirhmax')], AIC))))

#' We have a difference of about 11 in the AIC now, with the undirected
#' flow model preferred.

save.image('flows-checkpoint3.RData')
date()

#' Next we make diagnostic plots to verify that no assumptions are
#' badly violated and to see what parts of the data seem to lead to
#' differences in the fits.

fmint <- fortify(m$glmehmax, data=model.frame(m$glmehmax))
fund <- fortify(m$glmeundirhmax, data=model.frame(m$glmeundirhmax))
fdir <- fortify(m$glmedirhmax, data=model.frame(m$glmedirhmax))

cols <- c('stateF', '.fitted', '.resid', '.scresid', 'weekCent', 'logCmedDenseScaled')

D <- rbind(cbind(fmint[, cols], foi=fmint$logInfInt,
                 flow.measure=fmint$logInternalFlowScaled, flow.type='internal'),
           cbind(fund[, cols], foi=fund$logInfUnd,
                 flow.measure=fund$logUndirectedFlowScaled, flow.type='undirected'),
           cbind(fdir[, cols], foi=fdir$logInfDir,
                 flow.measure=fdir$logDirectedFlowScaled, flow.type='directed'))

#'
#+ fig.width=14, fig.height=14, dev.args=list(pointsize=18)
g <- ggplot(D, aes(.fitted, .scresid, colour=flow.type))
g <- g + geom_point() + facet_wrap(~stateF, scales='free') + geom_hline(yintercept=0)
g

g <- ggplot(D, aes(foi, .scresid))
g <- g + geom_point(colour='blue') + facet_grid(.~flow.type)
g <- g + geom_line(aes(group=stateF), alpha=0.4) + geom_smooth(method="loess", span=.9)
g

g <- ggplot(D, aes(flow.measure, .scresid))
g <- g + geom_point(colour='blue') + facet_grid(.~flow.type)
g <- g + geom_line(aes(group=stateF), alpha=0.4) + geom_smooth(method="loess", span=.9)
g

g <- ggplot(D, aes(logCmedDenseScaled, .scresid))
g <- g + geom_point(colour='blue') + facet_grid(.~flow.type)
g <- g + geom_line(aes(group=stateF), alpha=0.4) + geom_smooth(method="loess", span=.9)
g

g <- ggplot(D, aes(weekCent, .scresid))
g <- g + geom_point(colour='blue') + facet_grid(.~flow.type)
g <- g + geom_line(aes(group=stateF), alpha=0.4) + geom_smooth(method="loess", span=.9)
g

#' None of the models have large trends in the Pearson residuals
#' againt the predictors.
#'
#+ fig.width=14, fig.height=14, dev.args=list(pointsize=18)
ggplot(D, aes(stateF, .scresid, colour=flow.type)) + geom_boxplot() + geom_hline(yintercept=0)

ggplot(D, aes(flow.type, .scresid, colour=flow.type)) + geom_violin()

#' The distributions of Pearson residuals also look similar.

tmpf <- function(x) ranef(x)$state[, '(Intercept)']
re <- sapply(m[c('glmehmax', 'glmeundirhmax', 'glmedirhmax')], tmpf)
rownames(re) <- rownames(ranef(m$glmehmax)$stateF)
dotplot(re, pch=c('i', 'u', 'd'))

#' The conditional modes for the random effects are generally similar.

lapply(m[c('glmehmax', 'glmeundirhmax', 'glmedirhmax')], summary)

#' The estimated dispersion parameters and state group variance
#' provides one way of summarizing the difference in the models.

