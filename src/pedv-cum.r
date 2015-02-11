#+setup, include=FALSE, cache=FALSE
library(knitr)
opts_chunk$set(fig.path='figure-pedv-cum/', fig.align='center', fig.show='hold')
options(replace.assign=TRUE,width=80)
Sys.setlocale("LC_TIME", "C") #Needed for identical()
Sys.setlocale("LC_COLLATE", "C")

#' ## Packages
#+
library(car)
library(c060)
library(geoR)
library(ggplot2)
library(glmnet)
library(Hmisc)
library(igraph)
library(maps)
library(plyr)
print(sessionInfo())
set.seed(4253)

#' ## Data loading
#'
#' Farm counts by operation types
#+
tmpf <- function(){
  censPath <- file.path('table-25-hogs-and-pigs-by-operation-type.csv')
  cens <- read.csv(censPath,strip.white=TRUE, na.strings='(D)',
                   stringsAsFactors=FALSE)
  # The introduction to the reports says '-' represents 0 and (D) means deleted for privacy
  tmpf <- function(x) {  
      if(is.character(x)){
          ret <- ifelse(x=='-', 0, x)
          type.convert(ret)
      }else{
          x
      }
  }
  cens <- lapply(as.list(cens), tmpf)
  cens <- data.frame(cens)
  test <- cens$GEO %in% state.name
  ind <- grep('.*Farms|GEO|ITEM', colnames(cens))
  stateCounts <- cens[test, ind]

  key <- match(stateCounts$GEO, state.name)
  stateCounts$abb <- state.abb[key]

  test <- stateCounts$ITEM == "Total inventory \\ Farms with 1 to 24"
  smallFarms <- stateCounts[test,]
  test <- stateCounts$ITEM == "Total inventory"
  allFarms <- stateCounts[test,]
  
  stopifnot(allFarms$GEO == smallFarms$GEO)
  ind <- grep('.*Farms', colnames(allFarms))
  nonSmallFarms <- allFarms
  nonSmallFarms[,ind] <- allFarms[,ind] - smallFarms[,ind]
  nonSmallFarms$ITEM <- NULL
  
  nonSmallFarms
}
stateData <- tmpf()

#+
tmpf <- function(sd){
    ind <- grep('.*Farms', colnames(sd))
    dists <- sd[,ind]
    rownames(dists) <- sd$abb
    dists
}
dists <- tmpf(stateData)

#' County land areas
#+
data(county.fips, package='maps')
fp <- file.path('2013_Gaz_counties_national.txt')
rd <- read.delim(fp)
rd <- rd[, c('GEOID', 'ALAND')]

#' Hog counts by county
#+
tmpf <- function(){
  censPath <- file.path('table-12-hogs-and-pigs-by-county.csv')
  cens <- read.csv(censPath,strip.white=TRUE, na.strings='(D)',
                   stringsAsFactors=FALSE)
# The introduction to the reports says '-' represents 0 and (D) means deleted for privacy
  tmpf <- function(x) {  
      if(is.character(x)){
          ret <- ifelse(x=='-', 0, x)
          type.convert(ret)
      }else{
          x
      }
  }
  cens <- lapply(as.list(cens), tmpf)
  cens <- data.frame(cens)
  test <- cens$STCOFIPS %in% county.fips$fips
  cens <- cens[test,]
  ind <- grep('farms, 2007)$', cens$ITEM)
  cens <- cens[ind,]  
  cens
}
countyData <- tmpf()

#' ERS farm resource regions
#+
regs <- read.csv('reglink.csv', skip=2, colClasses=c(NA, NA, 'NULL', 'NULL'))


#' Derive state-level summaries of county-level data
#+
ctyTots <- ddply(countyData, 'STCOFIPS', summarize, totalFarms=sum(DATA), STFIPS=STFIPS[1],
                 smallFarms = sum(DATA["Inventory \\ Total hogs and pigs \\ Farms by inventory \\ 1 to 24 (farms, 2007)" == ITEM]))
ctyTots$nonSmallFarms <- with(ctyTots, totalFarms - smallFarms)
mg <- merge(ctyTots, rd, by.x='STCOFIPS', by.y='GEOID')
mg <- merge(regs, mg, by.x='Fips', by.y='STCOFIPS')
mg$kmsq <- mg$ALAND/1e6
mg$nonSmallDense <- mg$nonSmallFarms / mg$kmsq
stateCty <- ddply(mg, 'STFIPS', summarize, totNonSmall=sum(nonSmallFarms),
                  mDense=mean(nonSmallDense),
                  cmDense=mean(nonSmallDense[nonSmallFarms > 0]),
                  medDense=median(nonSmallDense), maxDense=max(nonSmallDense),
                  reg1=weighted.mean(ERS.resource.region==1, w=nonSmallFarms),
                  reg2=weighted.mean(ERS.resource.region==2, w=nonSmallFarms),
                  reg3=weighted.mean(ERS.resource.region==3, w=nonSmallFarms),
                  reg4=weighted.mean(ERS.resource.region==4, w=nonSmallFarms),
                  reg5=weighted.mean(ERS.resource.region==5, w=nonSmallFarms),
                  reg6=weighted.mean(ERS.resource.region==6, w=nonSmallFarms),
                  reg7=weighted.mean(ERS.resource.region==7, w=nonSmallFarms),
                  reg8=weighted.mean(ERS.resource.region==8, w=nonSmallFarms),
                  reg9=weighted.mean(ERS.resource.region==9, w=nonSmallFarms))
data(state.fips, package='maps')
key <- match(stateCty$STFIPS, state.fips$fips)
stateCty$abb <- state.fips$abb[key]
                                   
#' Case data
#+
caseData <- list()
dataDir <- file.path('.')
tmpf <- function(){
  fn <- file.path(dataDir, 'PEDvweeklyreport-state-age-cummulative-01-08-14.csv')
  ret <- read.csv(fn)
  manuallyChecked <- structure(
      list(State = structure(c(1L, 2L, 3L, 4L, 5L, 6L, 7L,
               8L, 9L, 10L, 11L, 12L, 13L, 14L, 15L, 16L, 17L, 18L, 19L, 20L,
               22L, 23L, 21L), .Label = c("CA", "CO", "IA", "IL", "IN", "KS",
                                   "KY", "MD", "MI", "MN", "MO", "NC", "NE", "NY", "OH", "OK", "PA",
                                   "SD", "TN", "TX", "UNKNOWN", "WI", "WY"), class = "factor"),
           Suckling = c(0L, 11L, 132L, 22L, 10L, 52L, 3L, 0L, 2L, 20L,
               2L, 100L, 0L, 2L, 11L, 55L, 4L, 1L, 0L, 6L, 2L, 1L, 1L),
           Nursery = c(1L, 6L, 149L, 10L, 12L, 14L, 0L, 0L, 3L, 42L,
               7L, 50L, 0L, 0L, 20L, 72L, 8L, 0L, 0L, 4L, 1L, 0L, 2L),
           Grower.Finisher = c(0L, 4L, 251L, 13L, 11L, 41L, 0L, 0L, 2L, 27L, 1L, 67L, 0L, 0L,
               10L, 45L, 11L, 0L, 6L, 12L, 1L, 0L, 1L),
           Sow.Boar = c(0L, 8L, 39L, 9L, 5L, 30L, 0L, 0L, 1L, 5L, 3L, 51L, 0L, 0L, 4L,
               37L, 2L, 0L, 0L, 0L, 0L, 0L, 0L),
           Unknown = c(0L, 1L, 108L, 17L, 20L, 8L, 1L, 1L, 1L, 112L, 5L, 38L, 5L, 0L, 14L, 31L,
               0L, 3L, 0L, 6L, 1L, 0L, 10L)),
      .Names = c("State", "Suckling", "Nursery", "Grower.Finisher", "Sow.Boar", "Unknown"),
      class = "data.frame", row.names = c(NA,-23L))
  stopifnot(identical(manuallyChecked, ret))
  nonreporting <- setdiff(state.abb, c(as.character(ret$State)))
  df <- data.frame(State=nonreporting, Suckling=0, Nursery=0, Grower.Finisher=0, Sow.Boar=0, Unknown=0)
  rbind(ret, df)
}
counts <- tmpf()

#' AK and HI are special cases. Let's just look at the coniguous 48
#' states.
#+
exclude <- c('AK', 'HI')
test <- !(counts$State %in% exclude)
counts <- counts[test,]
cts <- reshape(counts, v.names='cases', idvar='State', varying=2:6,
               direction='long', timevar='age', times=names(counts)[2:6])
ctsTot <- ddply(cts, 'State', summarize, cases=sum(cases))

#' Not a very strong relationship, but there is agreement in that
#' states not reporting any cases had much smaller changes in pig
#' litter

#' Shipment flows data
#+
flowMatGet <- function(){
    fp <- file.path('shipment-flows-origins-on-rows-dests-on-columns.csv')                    
    ep <- read.csv(fp, row.names=1)
    ep <- t(ep)
    key <- na.omit(match(unique(counts$State), colnames(ep)))
    ep <- ep[key, key]
    data.matrix(ep)
}
flowMat <- flowMatGet()



#' Calculate eigenvector centrality. All states except VT (no
#' outshipments) and CT (no inshipments) are in one strongly connected
#' compoent, so it doesn't seem necessary to bother with Katz
#' centrality or the like.
#+
vecGet <- function(symmetric=FALSE, normOutput=FALSE){
    getStateDist <- function(oddsInstate=9, pow=1, ep=flowMat){
        diag(ep) <- rowSums(ep) * 9
        ept <- ep^pow
        if(symmetric & normOutput){
            ret <- rowSums(ept)
            ret <- ret/sum(ret)
        } else {
            if(symmetric){
                ept <- (ept + t(ept))/2
            }
            if(normOutput){
                M <- t(ept)
                kout <- rowSums(M)
                kout <- ifelse(kout==0,1,kout)
                M <- M / kout
                ept <- t(M)
            }
            eig <- eigen(ept, symmetric=symmetric)
            ret <- eig$vector[,1]
            stopifnot(Im(ret)==0)
            ret <- Re(ret)
            names(ret) <- rownames(ep)
            ret <- ret/sum(ret)
        }
        ret
    }
    powers <- as.list(2^(-1*(0:4)))
    tmpf <- function(x) getStateDist(pow=x)
    vecs <- lapply(powers, tmpf)
    vecs <- do.call(rbind, vecs)
    vecs <- t(vecs)
    vecs <- data.frame(State=rownames(vecs), vecs)
    lab <- paste('isSym', symmetric, 'andIsNormed', normOutput, 'andPow', sep='')
    tmpf <- function(x) paste(lab, x, sep='')
    names(vecs)[-1] <- sapply(powers, tmpf)    
    vecs
}
args <- expand.grid(symmetric=c(TRUE, FALSE), normOutput=c(TRUE, FALSE))
vecs <- list()
for(i in 1:nrow(args)){
    vecs[[i]] <- do.call(vecGet, args[i,])
}
vecs <- Reduce(merge, vecs)

#' Balance sheet data
#'
#' Hawaii has a dash for inshipments, which we will assumes means a
#' negligible number.
#'
#' We distribute the IdahoWashington row among ID and WA in proportion
#' to the number of nonsmall farms in those states.
#+
bal <- read.csv('state-hogBalanceSheetDec2011Dec2012.csv',
                na.strings='-', row.names=1, colClasses=c(state2='NULL'))
recWA <- recID <- bal['IdahoWashington',]
iWA <- which(stateCty$abb == 'WA')
iID <- which(stateCty$abb == 'ID')
pWA <- with(stateCty, totNonSmall[iWA]/(totNonSmall[iWA] + totNonSmall[iID]))
recWA <- recWA*pWA
recID <- recID*(1-pWA)
bal <- rbind(bal, "Washington"=recWA, "Idaho"=recID)
key <- match(rownames(bal), state.name)
bal$State <- state.abb[key]
bal <- bal[!is.na(bal$State),]
bal['Hawaii','inshipments'] <- 0

#' Derive sampling weights
#'
#' We are calculating sampling weights for each age class in each
#' state based on random sampling of the farms in each state and then
#' random sampling of the pigs on each sampled farm. Census data for
#' each state provides the number of operations of each type in each
#' state. The age distribution of pigs on each farm are calculated as
#' follows.
#' 
#' Parameters taken from Stalder (2013) "Pork industry productivity analysis".
#+
agtGet <- function(){
    littersPerSowWeek <- 2.31/52
    meanLitter <- 10.3
    meanWeaningAge <- 21.5/7
    meanGoToFeederAge <- (21.5 + 46.0)/7
    meanSlaughterAge <- (21.5 + 46.0 + 121.5)/7
    sowBoarRatio <- 0
    nSow <- 100
    nBoar <- ifelse(sowBoarRatio > 0, nSow / sowBoarRatio, 0)
    nSowBoar <- nSow + nBoar
    flowRate <- nSow * littersPerSowWeek * meanLitter
    nSuckling <- flowRate * meanWeaningAge
    nNursery <- flowRate * (meanGoToFeederAge - meanWeaningAge)
    nFeeder <- flowRate * (meanSlaughterAge - meanGoToFeederAge)
    agt <- diag(c(nSuckling, nNursery, nFeeder, nSowBoar))
    agt
}
agt <- agtGet()

#' We've computed an age distribution based on the birth rate of sows
#' and the residence time in each production stage. We now create an
#' age distribution for each type of farm based on the age classes
#' present in that farm.
#+ 
M0 <- cbind(c(1, 0, 0, 1),
            c(1, 1, 1, 1),
            c(0, 0, 1, 0),
            c(1, 1, 0, 1),
            c(0, 1, 0, 0))
W <- agt %*% M0
N <- colSums(W)
W <- scale(W, center=FALSE, scale=N)
dists$Other.Farms <- NULL

#' Now sample weights are computed for each state using distribution
#' of farm types in each state.
#+
wts <- W %*% t(dists)
wts <- data.frame(age=c("Suckling", "Nursery", "Grower.Finisher", "Sow.Boar"),
             wts)
wts <- reshape(wts, v.names='samplingWeight', idvar='age', varying=2:51,
               direction='long', timevar='State', times=names(wts)[2:51])
rm(M0, W, N)

#' Merge predictors
#+
mg <- merge(wts, cts)
mg <- merge(mg, vecs)
mg <- merge(mg, bal)
mgAge <- mg
save(mgAge, file='mgAge.RData')
q('no')


#' Merge in data on spatial density of farms
#+
mg <- merge(ctsTot, vecs)
mg <- merge(mg, bal)
mg <- merge(mg, stateCty, by.x='State', by.y='abb')
mg$l2CasesOffset <- log2(mg$cases + 0.5)

#' Compare total case counts with change in pig litter
#+
df2 <- read.csv('pigsPerLitterDecemberThruFebrurary2013and2014.csv')
df2$change <- with(df2, pigLitter2014 - pigLitter2013)
key <- match(df2$state, tolower(state.name))
df2$abb <- state.abb[key]
df2$percentDecrease <- df2$change/df2$pigLitter2013 * -100
foo <- df2[!is.na(df2$abb),]
ord <- order(foo$percentDecrease)
dotchart(foo$percentDecrease[ord], labels=foo$abb[ord], xlab='Percent decrease in litter size',
         ylab='State')
v <- df2$percentDecrease[df2$state=='united states']
abline(v=v)
mtext('U.S.', at=v)
v <- df2$percentDecrease[df2$state=='other states']
abline(v=v)
mtext('Other states', at=v)
#
foo <- mg[, c('State', 'cases', 'totNonSmall')]
foo <- merge(foo, df2, by.x='State', by.y='abb', all.x=TRUE)
mg <- merge(mg, df2[, c('abb', 'percentDecrease')], by.x='State', by.y='abb', all.x=TRUE)
#
foo$cumInc <- foo$cases/foo$totNonSmall
with(foo, plot(cumInc ~ percentDecrease, type='n',
               xlab='Percent decrease in litter rate',
               ylab='cases / farms', frame=FALSE))
with(foo, text(percentDecrease, cumInc, labels=State))
summary(lm(cumInc~percentDecrease, data=foo))
ct <- cor.test(foo$percentDecrease, foo$cumInc, method='spearman')
text(2, .6, paste('Spearman rho =', round(ct$estimate, 2)))
text(2, .55, paste('p =', round(ct$p.value, 3)))
dev.off()
#
foo$pdImp <- foo$percentDecrease
foo$pdImp[is.na(foo$pdImp)] <- df2$percentDecrease[df2$state == 'other states']
foo$anyCases <- foo$cases > 0
foo$manyCases <- foo$cases > 8
boxplot(pdImp~ manyCases, data=foo, xlab='More than 8 cases',
            ylab='Percent decrease in pig litter')

wilcox.test(pdImp ~ anyCases, data=foo)
wilcox.test(pdImp ~ manyCases, data=foo)

#' #' Add lat longs
#+
key <- match(mg$State, state.abb)
mg$stateLong <- state.center$x[key]
mg$stateLat <- state.center$y[key]

#' Add predictor based on mean number of cases in states with shared
#' borderline
#+
borderMatGet <- function(){
    nbEdgelist <- read.csv('state_neighbors_fips.txt', header=FALSE)
    g <- graph.data.frame(nbEdgelist, directed=FALSE)
    g <- simplify(g)
    key <- match(V(g)$name, state.fips$fips)
    abb <- state.fips$abb[key]
    V(g)$name <- as.character(abb)
    states <- levels(mg$State)
    vids <- which(V(g)$name %in% states)
    g2 <- induced.subgraph(g, vids)
    nhood <- as.matrix(get.adjacency(g2))
}
borderMat <- borderMatGet()
nbWeightGet <- function(M, flowTrans=identity, normFlows=TRUE){
    M <- flowTrans(M)
    if(!normFlows){
        ret <- (M)
    } else {
        N <- rowSums(M)
        N <- ifelse(N > 0, N, 1)
        ret <- M/N
    }
    ret
}
nbCasesGet <- function(..., caseTrans=function(x) log2(x + 0.5), symmetrize=FALSE){
    W <- nbWeightGet(...)
    avgc <- tapply(mg$cases, mg$State, mean)
    key <- match(colnames(W), names(avgc))
    avgc <- avgc[key]
    stopifnot(names(avgc) == colnames(W))
    if(symmetrize){
        W <- (W + t(W))/2
    }else{
        W <- t(W)
    }
    avgNbC <- W %*% avgc
    avgNbC <- caseTrans(avgNbC)
    data.frame(State=rownames(avgNbC), avgNbC=avgNbC)
}
foo <- nbCasesGet(M=borderMat, normFlows=TRUE)
names(foo)[2] <- 'carNbC'
mg <- merge(mg, foo)
tmpf <- function(y){
   ret <- nbCasesGet(M=flowMat, flowTrans=function(x) x^y, normFlows=FALSE)
   names(ret)[2] <- paste('casesFlowNb', y, sep='')
   ret
}
powers <- as.list(2^(-1*(0:4) ))
foo <- lapply(powers, tmpf)
for(i in seq_along(foo)){
    mg <- merge(mg, foo[[i]])
}
rm(powers, foo)

#' ## Analysis
#' 
#' Examine colinearity
#' 
#+ , fig.width=11.11, fig.height=11.11, out.width=800, out.height=800
getPreds <- function(df, transform=c('identity', 'log'),
                     predNames=c('totNonSmall', 'mDense',
'cmDense', 'medDense', 'maxDense', 'inventory2012', 'pigCrop',
'inshipments', 'marketings')){
    transform <- match.arg(transform)
    df <- df[, predNames]
    x <- data.matrix(df)
    bot <- min(x[x[,'medDense'] > 0, 'medDense']) / 2
    if(transform=='log'){
        x[, 'medDense'] <- x[,'medDense'] + bot
        log(x)
    } else {
        x
    }
}
x <- getPreds(df=mg, transform='id')
xlog <- getPreds(df=mg, transform='log')
scatterplotMatrix(xlog, transform=FALSE, cex.labels=1.5, pch=16,
                  col=c(palette()[3:2], rgb(0,0,0,0.25)), cex.axis=2)

#' Add in the eigenvector centrality predictors, exept for the one set
#' that predicts all of the cases in Vermont (i.e.,
#' isSymFALSEisNormedTRUE* predictors). Vermont only ships to itself,
#' so a random walk gets trapped there in those cases. Clearly that is
#' a bad model so we'll remove them.
#' 
#+
getOtherPreds <- function(D, X){
    xEig <-  D[, grep('^isSym', colnames(D))]
    xEig <- xEig[, -grep("^isSymFALSEandIsNormedTRUE", colnames(xEig))]
    sel <- D[, grep('^reg|^casesFlowNb|^carNb', colnames(D))]
    ret <- data.matrix(cbind(X, sel, xEig))
    rownames(ret) <- as.character(D$State)
    ret
}
xx <- getOtherPreds(D=mg, X=x)
xxlog <- getOtherPreds(D=mg, X=xlog)

#' We have too many predictors for a scatterplot matrix to be
#' useful. Let's do PCA.
#+
pc <- prcomp(xx, scale=TRUE)
pve <- summary(pc)$importance[2,]
cve <- summary(pc)$importance[3,]
par(mfrow=c(1,2))
plot(pve, xlab='Principle component', ylab='PVE', type='o', col='blue')
plot(cve, xlab='Principle component', ylab='Cumulative PVE', type='o', col='blue')
par(mfrow=c(1,1))
biplot(pc)
biplot(pc, c(1,3))

#' The plot is crowded, but you can see that IA is extreme on PC1 and
#' an outlier on PC2. NC on the other hand, is extreme on PC1 and an
#' outlier on PC3.
#'
#' Next let's check the transformed-only data.
#+
pc <- prcomp(xxlog, scale=TRUE)
pve <- summary(pc)$importance[2,]
cve <- summary(pc)$importance[3,]
par(mfrow=c(1,2))
plot(pve, xlab='Principle component', ylab='PVE', type='o', col='blue')
plot(cve, xlab='Principle component', ylab='Cumulative PVE', type='o', col='blue')
par(mfrow=c(1,1))
biplot(pc)
biplot(pc, c(1,3))

#' The data is roughly three dimensional with most of the
#' variation in one dimension but the data seems more evenly spread
#' over the space and perhaps more likely to provide a useful model.
#'
#' Model of all cases
#+
alphas <- 0:5/5
tmpf <- function(alpha){
    cv.glmnet(x=xxlog, y=mg$l2CasesOffset, alpha=alpha)
}
cvsAll <- lapply(alphas, tmpf)

#' PCA of selected variables
#+
b <- drop(coef(cvsAll[[6]]))[-1]
pcAll1 <- prcomp(xxlog[, b>0], scale=TRUE)
biplot(pcAll1)

#' Check prediction error by alpha
#+ 
plot3Alphas <- function(fitList, alphas, which=4:6, ...){
 ## Adpated from http://www.stanford.edu/~hastie/glmnet/glmnet_alpha.html
 opar <- par(mfrow = c(2, 2))
 cv <- fitList[which]
 a <- alphas[which]
 tmpf <- function(x,y) {
     main <- substitute(expression(paste(alpha, '=', y)), list(y=y))
     plot(x, sub=main, ...)
 }
 invisible(mapply(tmpf, x=cv, y=a))
 ls <- lapply(cv, '[[', 'lambda.1se')
 plot(log(cv[[1]]$lambda), cv[[1]]$cvm, pch = 19, col = "red", xlab = "log(Lambda)",
         ylab = cv[[1]]$name, ...)
 abline(v=log(ls[[1]]), col="red")
 points(log(cv[[2]]$lambda), cv[[2]]$cvm, pch = 19, col = "grey")
 abline(v=log(ls[[2]]), col="grey")
 points(log(cv[[3]]$lambda), cv[[3]]$cvm, pch = 19, col = "blue")
 abline(v=log(ls[[3]]), col="blue")
 legend("topleft", legend = paste('alpha =', a), pch = 19,
           col = c("red", "grey", "blue"))
 par(opar)
}
plot3Alphas(cvsAll, alphas, ylim=c(0,10))

#' There's no alpha with a decisive predicition error advantage here.
#'
#' Some Diagnostics
#+
m <- cvsAll[[6]]
yhat <- predict(m, type='response', new=xxlog)
lmfit <- lm(mg$l2CasesOffset~0, offset=yhat)
layout(matrix(1:4, 2))
plot(lmfit, which=1:3)

#' The trend in the mean of the residuals and the small variance for
#' low predictions look a little troublesome and suggest we try a
#' hurdle model.
#'
#' Logistic regression for prob non-zero
#+
mg$anyCases <- factor(mg$cases > 0)
tmpf <- function(alpha){
    cv.glmnet(x=xxlog, y=mg$anyCases, alpha=alpha, family='binomial')
}
cvsBin <- lapply(alphas, tmpf)
lsBin <- lapply(cvsBin, '[[', 'lambda.1se')
mapply(function(x, y) x$cvm[which(x$lambda == y)], x=cvsBin, y=lsBin) 
layout(matrix(1:6, ncol=3))
tmpf <- function(x,y) {
    main <- substitute(expression(paste(alpha, '=', y)), list(y=y))
    plot(x, sub=main)
}
invisible(mapply(tmpf, x=cvsBin, y=alphas))

#' Visualization of coefficients
#' 
#+ 
coefplot <- function(fitList, preds, which=3:6){
  coefMat <- sapply(fitList, function(x) as.numeric(coef(x)))[-1,]
  rownames(coefMat) <- rownames(coef(fitList[[1]]))[-1]
  colnames(coefMat) <- paste('alpha', alphas, sep='')
  coefMat <- t(coefMat)
  par(mfrow=c(1,2),
      mar=c(20,4,1,1))
  barplot(coefMat[which,], beside=T, las=2, ylab='Coefficients')
  scales <- apply(preds, 2, sd)
  coefMatStd <- coefMat
  stopifnot(names(scales)==colnames(coefMat))
  for(i in 1:nrow(coefMatStd)){
      coefMatStd[i,] <- coefMatStd[i,] * scales
  }
  test <- abs(colSums(coefMatStd[which,])) > 0.01
  barplot(coefMatStd[which,test], beside=T, las=2, ylab='Standardized coefficients')
  par(mfrow=c(1,1),
      mar=c(5,4,4,1) + 0.1)

}
coefplot(cvsBin, xxlog, which=4:6)

#' Diagnostics for logistic fit
#+
mbin <- cvsBin[[6]]
yhatb <- predict(mbin, type='link', new=xxlog)
anyCasesInt <- as.integer(as.logical(mg$anyCases))
glmfit <- glm(anyCasesInt~0, offset=yhatb, family=binomial)
layout(matrix(1:4, 2))
plot(glmfit, which=1:3)
Ep <- residuals(glmfit, type='pearson')

#' Next let's check for outliers. Because of our tricky way of
#' producing glmfit, this will be a t-test on the deviance residuals
#' instead of the usual studentized residuals, but that seems unlikely
#' to be an issue.
#+
outlierTest(glmfit)

#' Bubble plot of residuals on U.S. map
#+
df <- data.frame(res=as.vector(Ep), x=mg$stateLong, y=mg$stateLat, state=mg$State)
all_states <- map_data("state")
p <- ggplot()
p <- p + geom_polygon(data=all_states, aes(x=long, y=lat, group = group),colour="white", fill="grey10" )
p <- p + geom_point(data=df, aes(x=x, y=y, size=abs(res)*2, color=as.factor(sign(res))))
p

#' The colors seem to go together somewhat, which means there may be
#' some spatial correlations.
#+
#' Get great circle distance
#+
cx <- df$x
cy <- df$y
E <- df$res
n <- length(cx)
centerDists <- matrix(nrow=n, ncol=n)

#' Calculates the geodesic distance between two points specified by radian latitude/longitude using the
#' Haversine formula (hf)
#' source: http://www.r-bloggers.com/great-circle-distance-calculations-in-r/
#+
gcd.hf <- function(long1, lat1, long2, lat2) {
      R <- 6371 # Earth mean radius [km]
        delta.long <- (long2 - long1)
        delta.lat <- (lat2 - lat1)
        a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
        c <- 2 * asin(min(1,sqrt(a)))
        d = R * c
        return(d) # Distance in km
  }
deg2rad <- function(deg) return(deg*pi/180)
getCenterDist <- function(state1,state2){
    long1 <- deg2rad(cx[state1])
    long2 <- deg2rad(cx[state2])
    lat1 <- deg2rad(cy[state1])
    lat2 <- deg2rad(cy[state2])
    gcd.hf(long1, lat1, long2, lat2)
}
for(i in seq_len(n)){
    for(j in seq_len(n)){
        centerDists[i,j] <- getCenterDist(i, j)
    }
}
colnames(centerDists) <- rownames(centerDists) <- df$state

#' Now lets find points on 2-d plane that satisfy distances
#+
loc <- cmdscale(centerDists)
x <- -loc[, 1]
y <- loc[,2]
plot(x, y, type='n', asp=1)
text(x, y, rownames(loc), cex=.6)

#' Now lets check for correlations in space
#+
sE <- E/sd(E)
vg <- variog(coords=loc, data=sE, max.dist=2000)
olsFit <- variofit(vg, cov.model='exponential')
likFit <- likfit(coords=loc, data=sE, cov.model='exponential', ini.cov.pars=c(.1, 194))
par(mfrow=c(1,1))
plot(vg)
lines(olsFit)
lines(likFit, lwd=2)

#' The densley packed, spatially autocorrelated states in the North
#' East might be exerting too much influence. Let's try to identify
#' groups of close states and combine their data into a single
#' observation.
#'
#+
clus <- hclust(as.dist(centerDists), method='single')
plot(clus)
abline(h=175)

#' If we cut the single-linkage tree at 175, it seems lead to sensible
#' groupings. So let's create a new data set with rows from those
#' states averaged together.
#+
groups <- cutree(clus, h=175)
tmpf <- function(x) names(groups)[groups==x]
groupL <- lapply(unique(groups), tmpf)
tmpf <- function(x){
  rows <- mg$State %in% x
  M <- data.matrix(mg[rows,-1])
  colMeans(M)
}
rowsAvg <- lapply(groupL, tmpf)
mgAvg <- do.call(rbind, rowsAvg)
rownames(mgAvg) <- sapply(groupL, paste, collapse='-')
mgAvg <- data.frame(mgAvg)
mgAvg[, 'anyCases'] <- factor(mgAvg[, 'cases'] > 0)

xlogAvg <- getPreds(df=mgAvg, transform='log')
xxlogAvg <- getOtherPreds(D=mgAvg, X=xlogAvg)
rownames(xxlogAvg) <- rownames(mgAvg)

mbinavg <- cv.glmnet(x=xxlogAvg, y=mgAvg$anyCases, alpha=1, family='binomial')
mbinavg0.8 <- cv.glmnet(x=xxlogAvg, y=mgAvg$anyCases, alpha=0.8, family='binomial')

cbind(drop(coef(mbin)), drop(coef(mbinavg)))
cbind(drop(coef(cvsBin[[5]])), drop(coef(mbinavg0.8)))

#' That weighting has resulted in minor changes to the estimated
#' coefficients and some changes to the selected variables. Notably,
#' median farm density has been selected when alpha=0.8 for the
#' averaged data, and inshipments has been selected instead of
#' eigenvector centrality when alpha=1.
#+
yhatavg <- predict(mbinavg, type='link', new=xxlogAvg)
anyCasesIntAvg <- as.integer(as.logical(mgAvg$anyCases))
glmfitavg <- glm(anyCasesIntAvg~0, offset=yhatavg, family=binomial)
layout(matrix(1:4, 2))
plot(glmfitavg, which=1:3)
Epavg <- residuals(glmfitavg, type='pearson')
sEpavg <- Epavg/sd(Epavg)
locAvg <- with(mgAvg, cbind(x=stateLong, y=stateLat))
vgavg <- variog(coords=locAvg, data=sEpavg, max.dist=2000)
plot(vgavg)

#' The variogram looks much better.
#'
#' Next we model the positive counts
#+
test <- as.logical(mgAvg$anyCases)
sub <- mgAvg[test, ]
subx <- xxlogAvg[test, ]
subl <- locAvg[test,]
alphas <- 0:5/5
tmpf <- function(alpha){
    cv.glmnet(x=subx, y=log(sub$cases), alpha=alpha, nfolds=5)
}
cvsPos <- lapply(alphas, tmpf)

#' The diagnostics
#+
m <- cvsPos[[6]]
yhat <- predict(m, type='response', new=subx)
lmfit <- lm(log(sub$cases)~0, offset=yhat)
layout(matrix(1:4, 2))
plot(lmfit, which=1:3)
Ep <- residuals(lmfit, type='pearson')
vg <- variog(coords=subl, data=Ep/sd(Ep), max.dist=2000)
olsFit <- variofit(vg, cov.model='exponential')
likFit <- likfit(coords=subl, data=Ep/sd(Ep), cov.model='exponential', ini.cov.pars=c(.01, 400))
plot(vg)
lines(olsFit)
lines(likFit, lwd=2)

#' There's a trend in the residuals vs fitted, but that's probably
#' just showing us that the shrinkage is protecting us from
#' overestimating the coefficients. The semivariogram looks good.
#'
#' Let's see if there are any concerning patterns in plots of the
#' residuals versus predictors.
#+
for(i in colnames(subx)){
    plot(Ep~subx[, i], main=i, pch=16, col=rgb(0,0,0,0.7), xpd=TRUE)
    text(subx[, i], jitter(Ep, amount=0.1), labels=sub$State)
}


#' It's a little surprising that the region 1 doesn't get selected,
#' but that's probably because it is correlated with population size.
#+
plot(reg1~inventory2012,data=mgAvg)
cor.test(mgAvg$reg1, mgAvg$inventory2012)

#' The eigenvector predictors often seem to be grouped together with
#' the balance sheet ones. That is to be expected, given that the
#' flows are derived in part from the balance sheet data. Perhaps it
#' is not worth presenting the eigenvector results, which will take a
#' while to properly describe and not really change the
#' conclusions. So, let's do a final round of fitting without them.
#' 
#' ## Fits without eigenvector predictors
#+
xf <- xxlogAvg[, -grep('^isSym', colnames(xx))]
xfpc <- prcomp(xf, scale=TRUE)
summary(xfpc)
biplot(xfpc)
biplot(xfpc, choices=c(1,3))

#' Balance sheet and farm count contibutes primarily to the first PC
#' and to the third PC to a lesser extent. Several region variables
#' make the largest contributions to the second and third PCs. The second PC
#' comprises primarily neighboring cases, density, and region. There
#' should be some power to detect the effect of these variables
#' independent of the population-size variables in the first PC.
#'
#+
test <- as.logical(mgAvg$anyCases)
subxf <- xf[test, ]
alphas <- 0:5/5
tmpf <- function(alpha){
    cv.glmnet(x=subxf, y=log(sub$cases), alpha=alpha, nfolds=5)
}
cvsPosF <- lapply(alphas, tmpf)
plot3Alphas(cvsPosF, alphas)
coefplot(cvsPosF, subxf, which=4:6)

#' Diagnostics for case counts
#+
m <- cvsPosF[[5]]
yhatcf <- predict(m, type='response', new=subxf)
lmfit <- lm(log(sub$cases)~0, offset=yhatcf)
layout(matrix(1:4, 2))
plot(lmfit, which=1:3)
Ep <- residuals(lmfit, type='pearson')
Estd <- Ep/sd(Ep)
vg <- variog(coords=subl, data=Estd, max.dist=2000)
olsFit <- variofit(vg, cov.model='exponential')
likFit <- likfit(coords=subl, data=Estd, cov.model='exponential', ini.cov.pars=c(.01, 300))
plot(vg, xlab='distance')
lines(olsFit)
lines(likFit, lwd=2)
par(mfrow=c(1,1))

#' Logistic regression for prob non-zero
#+
alphas <- 0:5/5
tmpf <- function(alpha){
    cv.glmnet(x=xf, y=mgAvg$anyCases, alpha=alpha, family='binomial')
}
cvsBinF <- lapply(alphas, tmpf)
plot3Alphas(cvsBinF, alphas)
coefplot(cvsBinF, xf, which=4:6)

#' Diagnostics for probability non-zero
#+
m <- cvsBinF[[5]]
yhatb <- predict(m, type='link', new=xf)
glmfit <- glm(anyCasesIntAvg~0, offset=yhatb, family=binomial)
layout(matrix(1:4, 2))
plot(glmfit, which=1:3)
Ep <- residuals(glmfit, type='pearson')
sEp <- Ep/sd(Ep)
vg <- variog(coords=locAvg, data=sEp, max.dist=2000)
plot(vg)
likFit <- likfit(coords=locAvg, data=sEp, cov.model='exponential', ini.cov.pars=c(.01, 300))
lines(likFit)

#' Pig-litter--decrease model
#' 
#+
alphas <- 0:5/5
mgAvg$largeDec <- mgAvg$percentDecrease > 2
mgAvg$largeDec[is.na(mgAvg$largeDec)] <- FALSE
tmpf <- function(alpha){
    cv.glmnet(x=xf, y=mgAvg$largeDec, alpha=alpha, family='binomial')
}
cvsBinFD <- lapply(alphas, tmpf)
plot3Alphas(cvsBinFD, alphas)
coefplot(cvsBinFD, xf, which=4:6)
m <- cvsBinFD[[6]]
yhatb <- predict(m, type='link', new=xf)
ldint <- as.integer(as.logical(mgAvg$largeDec))
glmfit <- glm(ldint~0, offset=yhatb, family=binomial)
layout(matrix(1:4, 2))
plot(glmfit, which=1:3)
Ep <- residuals(glmfit, type='pearson')
sEp <- Ep/sd(Ep)
vg <- variog(coords=locAvg, data=sEp, max.dist=2000)
plot(vg)
likFit <- likfit(coords=locAvg, data=sEp, cov.model='exponential',
                 ini.cov.pars=c(.01, 300))
lines(likFit)

#' Same variables selected, same relative sizes.  Let's check deviance
#' ratio to evaluate fit.
#'
#+
getFitStats <- function(x){
    lam <- x$lambda.1se
    ft <- x$glmnet.fit
    ind <- which(ft$lambda == lam)
    ind2 <- which(x$lambda==lam)
    res <- list(lambda=lam)
    res$a0 <- ft$a0[ind]
    res$dr <- ft$dev.ratio[ind]
    res$nzero <- x$nzero[ind]
    res$cvm <- x$cvm[ind2]
    res$cvsd <- x$cvsd[ind2]
    res$nobs <- ft$nobs
    res    
}
mods <- list(loglog=cvsPosF, litter=cvsBinFD, presence=cvsBinF)
tmpf <- function(x){
    foo <- sapply(x[5:6], getFitStats)
    foo <- t(foo)
    as.data.frame(foo)
}
stats <- lapply(mods, tmpf)
statsdf <- do.call(rbind, stats)
statsdf <- cbind(alpha=c(0.8, 1), statsdf)
statsdf <- round(data.matrix(statsdf), 2)

invisible(latex(statsdf, file='stats-table.tex'))

#' ## Stability selection
#'
#' Stability selection seems more appopriate for identifying important
#' variables than cross-validation. We will check for sensitivity to
#' alpha.

tmpf <- function(alpha,...){
    stabpath(alpha=alpha, weakness=1, steps=1000, ...)
}
tmpff <- function(spec){
    alphas <- c(0.01,.2,.5,.8,1)
    names(alphas)=paste('alpha', format(alphas), sep='=')
    lapply(alphas, tmpf, y=spec$y, family=spec$family, x=spec$x)
}
specs <- list(litterRateDecrease=list(y=mgAvg$largeDec, family='binomial', x=xf),
              cumCases=list(y=log(sub$cases), family='gaussian', x=subxf),
              anyCases=list(y=mgAvg$anyCases, family='binomial', x=xf))
par(mfrow=c(3,5))
res <- lapply(specs, tmpff)
stab <- lapply(res, function(x) lapply(x, plot, type='pcer', error=0.05))
tmpfff <- function(x) lapply(x, '[[', 'stable')
lapply(stab, tmpfff)

save.image('pedv-cum.RData')
