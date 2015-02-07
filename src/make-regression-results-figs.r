library(lme4)
library(coefplot2)
library(glmmADMB)
library(Hmisc)

load('checkpoint3.RData')

## Refit with glmmadmb due to warnings about gradients above tolerance with lme4 fits

modNames <- c('glmeundirhmax', 'glmedirhmax', 'glmehmax')
fFin <- f[modNames]

tmpf <- function(x) {
    test <- !is.na(om$cases1WksAgo)
    data <- om[test, ]
    glmmadmb(x, data=data, family='nbinom')
}
mFin <- lapply(fFin, tmpf)

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
    longnames <- c(flow='Scaled swine flow', inf='Cases last week', week='Scaled week', dense='Scaled farm density', baseline='Baseline risk') 
    pal <- c('black', 'orange', 'blue')
    pch <- 15:17
    for(i in seq_len(nrow(ests))){
        coefplot2(ests[i, ], sds[i, ], varnames=longnames[colnames(ests)],
                  offset=(i-1)/10, col=pal[i], add=i>1, pch=pch[i],
                  main="Regression estimates")
    }
    legend('topright', legend=c('undirected', 'directed', 'none'), col=pal, pch=pch, title='Interstate flow')
}
#tmpf()

pdf('coefplot.pdf', width=5,height=4)
tmpf()
dev.off()
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

tab <- data.frame(baseline=I('no'), int=unname(int), disp, resd, df, llik)
rownames(tab) <- names(modList)
tab[c('glmeundirhmax', 'glmedirhmax', 'glmehmax'), 'baseline'] <- 'yes'
tab[, 'df'] <- tab[, 'df'] + ifelse(tab[, 'baseline'] == 'yes', 1, 0)
aic <- 2*-tab[, 'llik'] + 2*tab[, 'df']
tab <- cbind(tab, deltaAIC = aic - min(aic))
flowTerm <- c(glmeundirhmax='undirected', glmedirhmax='directed', glmehmax='internal',
              nbme='internal', nbmeN='none')
tab <- data.frame(flowTerm=flowTerm[rownames(tab)], tab)

tf <- format.df(tab, cdec=c(0,0,1,2,2,0,1,1))
align <- paste(attr(tf, "col.just"), collapse="|")
align <- paste('|', align, '|', sep='')
longnames <- c(flowTerm='Flow term', baseline='Fit eta', disp='hat theta',
               df='d.f.', llik='log lik.', deltaAIC='Delta AIC', resd='hat sigma',
               int='Intercept')
colnames(tf) <- longnames[colnames(tf)]
cat(latexTabular(tf, helvetica=FALSE, align=align, translate=FALSE), file='tab.tex')

scales <- c(week=iqr(om$weekCent),
            internal=iqr(log(om$internalFlow)),
            cmedDense=iqr(log(om$cmedDense*om$nFarms*om$nFarms)),
            undirected=iqr(log(om$undirectedFlow)),
            directeted=iqr(log(om$directedFlow)))
sink(file='scales.txt')
print(scales)
sink()

flowTerm <- c(glmeundirhmax='Within-state + undirected between state',
              glmedirhmax='Within-state + directed between state',
              glmehmax='Within-state only',
              nbme='Within-state only, fixed other risks ',
              nbmeN='No flow, fixed other risks')

png('aic.png', width=7*res, height=4*res, res=res)
dotchart2(aic, labels=flowTerm[rownames(tab)], xlab='AIC (i.e., Estimated information loss)', dotsize=2)
dev.off()

library(ggplot2)
fund <- fortify(m$glmeundirhmax, data=model.frame(m$glmeundirhmax))
theme_set(new=theme_classic())

g <- ggplot(data=fund, aes(x=exp(logUndirectedFlowScaled), y=exp(.fitted)))
g <- g + geom_smooth(method='loess', alpha=0, size=2)
g <- g + geom_point(aes(x=I(runif(nrow(fund), min=-.1,max=.1) + exp(logUndirectedFlowScaled)),
                        y=cases), col='red', alpha=0.5)
g <- g + labs(x='Flow (swine / pairs / year)', y='Cases')
g <- g + coord_trans(y="log1p")
g
ggsave('flow-prediction.pdf', width=7, height=5, pointsize=18)

save.image('checkpoint4.RData')

                                             
