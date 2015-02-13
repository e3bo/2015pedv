#!/usr/bin/Rscript

library(lme4)
library(coefplot2)
library(glmmADMB)
library(Hmisc)
library(ggplot2)

load('flows-checkpoint3.RData')

## Refit with glmmadmb due to warnings about gradients above tolerance with lme4 fits

modNames <- c('glmeundirhmax', 'glmedirhmax', 'glmehmax')
fFin <- f[modNames]

tmpf <- function(x) {
    test <- !is.na(om$cases1WksAgo)
    data <- om[test, ]
    glmmadmb(x, data=data, family='nbinom')
}
mFin <- lapply(fFin, tmpf)
save.image('flows-checkpoint4.RData')

## Plot of predicted cases vs flow

scales <- c(week=iqr(om$weekCent),
            internal=iqr(log(om$internalFlow)),
            cmedDense=iqr(log(om$cmedDense*om$nFarms*om$nFarms)),
            undirected=iqr(log(om$undirectedFlow)),
            directeted=iqr(log(om$directedFlow)))
sink(file='scales.txt')
print(scales)
sink()

centers <- c(week=mean(om$weekCent),
            internal=mean(log(om$internalFlow)),
            cmedDense=mean(log(om$cmedDense*om$nFarms*om$nFarms)),
            undirected=mean(log(om$undirectedFlow)),
            directeted=mean(log(om$directedFlow)))

fund <- fortify(m$glmeundirhmax, data=model.frame(m$glmeundirhmax))
stopifnot(all.equal(attributes(fund$logUndirectedFlowScaled)[["scaled:center"]], centers[['undirected']]))
stopifnot(all.equal(attributes(fund$logUndirectedFlowScaled)[["scaled:scale"]][["75%"]], scales[['undirected.75%']]))

theme_set(new=theme_classic())
fund$log10UndirectedFlow <- log10(exp(scales['undirected.75%'] * as.numeric(fund$logUndirectedFlowScaled) + centers['undirected']))
fund$predCases <- exp(fund$.fitted)
fund$smoothpred <- predict(loess(predCases~log10UndirectedFlow, data=fund, span=.5))
fund$casesJittered <- runif(fund$cases, min=-1, max=1)*0.2 + fund$cases

g <- ggplot(data=fund, aes(x=log10UndirectedFlow, y=casesJittered))
g <- g + geom_point(col="#4477AA", alpha=.5)
g <- g + geom_line(aes(y=smoothpred), col='red')
g <- g + labs(x='Log10(transport flow)', y='Cases')
g <- g + coord_trans(y="log1p")

ggsave('flows-prediction.eps', width=4, height=4, pointsize=12, device=cairo_ps)
ggsave('flows-prediction.pdf', width=4, height=4, pointsize=12, device=cairo_pdf)

## Plots of fixed effects

ct <- lapply(mFin, function(x) coeftab(x)[, 1:2])

extractEsts <- function(pattern, type) {
  tmpf <- function(x) grep(pattern, rownames(x))
  inds <- lapply(ct, tmpf)
  res <- list()
  mapply(function(x, y) x[y, type], ct, inds)
}
  
patterns <- list(flow='Flow', dense='logCmedDenseScaled', week='weekCent', inf='Inf')
ests <- sapply(patterns, extractEsts, type='Estimate')
sds <- sapply(patterns, extractEsts, type='Std. Error')
weekIqr <- iqr(model.frame(m$glmedirhmax)$weekCent)
ests[, 'week'] <- ests[, 'week']*weekIqr
sds[, 'week'] <- sds[, 'week']*weekIqr

stopifnot(rownames(ests) == c('glmeundirhmax', 'glmedirhmax', 'glmehmax'))
ests <- cbind(ests, baseline=c(hund$maximum, hdir$maximum, hint$maximum))
sds <- cbind(sds, baseline=0)

tmpf <- function() {
    longnames <- c(flow='Scaled transport flow', inf='Cases last week',
                   week='Scaled week', dense='Scaled farm density', baseline='Baseline risk') 
    pal <- c("blue"="#4477AA", "red"="#DDCC77", "yellow"="#CC6677")
    pch <- 15:17
    for(i in seq_len(nrow(ests))){
        coefplot2(ests[i, ], sds[i, ], varnames=longnames[colnames(ests)],
                  offset=(i-1)/10, col=pal[i], add=i>1, pch=pch[i],
                  main="Regression estimates")
    }
    legend('topright', legend=c('undirected', 'directed', 'none'), col=pal, pch=pch, title='Interstate flow')
}

pdf('coefplot.pdf', width=5,height=4)
tmpf()
dev.off()
## This works better than making the ps with R for some reason
system('pdftops coefplot.pdf')

res <- 100
png('coefplot.png', width=5*res, height=4*res, res=res)
tmpf()
dev.off()


## Table of likelihoods, dispersion, parameters, random effects

getNBTheta <- function(x) {
    if(inherits(x, 'glmmadmb')) {
        summary(x)$alpha
    } else {
        as.numeric(strsplit(family(x)$family, split='\\(|\\)')[[1]][2])
    }
}
getLL <- function(x) as.numeric(logLik(x))
getSD <- function(x) as.numeric(sqrt(VarCorr(x)$stateF))
getdf <- function(x) attr(logLik(x), 'df')
getintercept <- function(x) fixef(x)['(Intercept)']

modList <- c(mFin, m[c('nbme', 'nbmeN')])
disp <- sapply(modList, getNBTheta)
llik <- sapply(modList, getLL)
resd <- sapply(modList, getSD)
df <- sapply(modList, getdf)
int <- sapply(modList, getintercept)

tmpf <- function() {
  tab <- data.frame(baseline=I('no'), int=unname(int), disp, resd, df, llik)
  rownames(tab) <- names(modList)
  tab[c('glmeundirhmax', 'glmedirhmax', 'glmehmax'), 'baseline'] <- 'yes'
  tab[, 'df'] <- tab[, 'df'] + ifelse(tab[, 'baseline'] == 'yes', 1, 0)
  aic <- 2*-tab[, 'llik'] + 2*tab[, 'df']
  tab <- cbind(tab, deltaAIC = aic - min(aic))
  flowTerm <- c(glmeundirhmax='undirected', glmedirhmax='directed', glmehmax='internal',
                nbme='internal', nbmeN='none')
  tab <- data.frame(flowTerm=flowTerm[rownames(tab)], tab)
  tf <- format.df(tab)
  align <- paste(attr(tf, "col.just"), collapse="|")
  align <- paste('|', align, '|', sep='')
  longnames <- c(flowTerm='\\bf Flow term', baseline='\\bf Fit $\\eta$',
                 disp='\\bf Intercept', df='\\bf d.f.', llik='\\bf Log lik.',
                 deltaAIC='\\bf $\\Delta$ AIC', resd='\\bf $\\hat{\\sigma}$',
                 int='{\\bf Intercept}')
  colnames(tab) <- longnames[colnames(tf)]
  tab <- latexTabular(tab, helvetica=FALSE, align=align, translate=FALSE,
                   cdec=c(0,0,1,2,2,0,1,1))
  tab <- gsub('\\\\multicolumn\\{1\\}\\{c\\}', '', tab)
  tab <- gsub('(\\\\)\\s?(\n)', '\\1\\\\hline\\2', tab)
  tab <- gsub('(\\n)(\\{\\\\bf Flow term\\})', '\\1\\\\hline\\2', tab)
  cat(tab, file='tab.tex')
  list(aic=aic, tab=tab)
}
tables <- tmpf()

## Dotplot of AICs

tmpf <- function() {
    tab <- tables$tab
    aic <- tables$aic
    flowTerm <- c(glmeundirhmax='Within-state + undirected between state',
                  glmedirhmax='Within-state + directed between state',
                  glmehmax='Within-state only',
                  nbme='Within-state only, fixed other risks ',
                  nbmeN='No flow, fixed other risks')
    labels <- flowTerm[rownames(tab)]
    dotchart2(tables$aic, labels=labels, xlab='AIC (i.e., Estimated information loss)', dotsize=2)
}
png('aic.png', width=7*res, height=4*res, res=res)
tmpf()
dev.off()

save.image('flows-checkpoint5.RData')

                                             
