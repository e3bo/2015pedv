load('common-data.RData')

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

GetRegressionData <- function(case.data, flows, fiAll, stateCty){
  observed <- ObservedMatchingInternalFlows(flows=flows, case.data=case.data)
  ReshapeAugmentObserved(observed=observed, fiAll=fiAll, flows=flows, stateCty=stateCty)
}

