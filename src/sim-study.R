#!/usr/bin/Rscript
library(methods) #for raster

target.mean.deg.grid <- seq(0.3, 1.2, by=0.1)
raster.cell.side.grid <- c(16000, 32000)

ag <- Map(sds::CreateAgents,
          target.mean.deg=target.mean.deg.grid[1],
          raster.cell.side.meters=raster.cell.side.grid,
          census.dilation=0.1)

net.nbs <- lapply(target.mean.deg.grid[-1], sds::GetNetNbs,
                  block.labels=ag[[1]]$adf$abb)
net.nbs <- c(list(ag[[1]]$net.nbs), net.nbs)
ag.data.inds <- expand.grid(ag.ind=seq_along(ag), net.nbs.ind=seq_along(net.nbs))

nsim <- 100
par.ranges <- list(nstarters=c(0, 10),
                   prep=c(0.01, 1),
                   size=c(0.7, 1),
                   rprob=c(0, 1),
                   seasonal.amplitude=c(0, 1),
                   starting.grid.x=c(0, 1),
                   starting.grid.y=c(0, 1),
                   tprob.outside=c(0, 0.1),
                   tprob.net=c(0, 1),
                   tprob.sp=c(0, 1))

des <- sensitivity::parameterSets(par.ranges=par.ranges,
                                  samples=nsim, method='sobol')
colnames(des) <- names(par.ranges)

df <- do.call(cbind, c(list(des), as.list(ag.data.inds)))
df <- as.data.frame(df)
df$target.mean.deg <- target.mean.deg.grid[df$net.nbs.ind]
df$raster.cell.side <- raster.cell.side.grid[df$ag.ind]

system.time(res <- Map(sds::SimulateAndSummarize,
                       agent.data=ag[df$ag.ind],
                       lags.sel=1,
                       net.nbs=net.nbs[df$net.nbs.ind],
                       nstarters=df$nstarters,
                       permutations=2,
                       prep=df$prep,
                       rprob=df$rprob,
                       size=df$size,
                       seasonal.amplitude=df$seasonal.amplitude,
                       starting.grid.nx=10,
                       starting.grid.ny=2,
                       starting.grid.x=df$starting.grid.x,
                       starting.grid.y=df$starting.grid.y,
                       tprob.net=df$tprob.net,
                       tprob.outside=df$tprob.outside,
                       tprob.sp=df$tprob.sp))

tstats <- lapply(res, '[[', 'mantel.tests')
dfsplt <- split(df, seq(nrow(df)))
resplt <- Map(cbind, dfsplt, tstats)
resplt <- Filter(function(x) ncol(x) > 7, resplt)
resall <- do.call(rbind, resplt)

test <- with(resall,
             symmetrize == TRUE & mat2.name == 'shipment' & method== 'spearman' &
               mat1.name == 'lag1')
sub <- resall[test, ]

input.pars <- c(names(par.ranges), 'target.mean.deg', 'raster.cell.side')
formula <- as.formula(paste('~', paste(input.pars, collapse='+')))
m <- DiceKriging::km(formula=formula, design=sub[, input.pars],
                     response=sub$r, nugget.estim=TRUE, covtype='matern3_2')

nmeta <- 1e4
extra.par.ranges <- list(target.mean.deg=range(target.mean.deg.grid),
                         raster.cell.side=range(raster.cell.side.grid))
all.par.ranges <- c(par.ranges, extra.par.ranges)

GetRandLHSDes <- function(n, ranges){
  X <- lhs::randomLHS(nmeta, length(ranges))
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
X1 <- GetRandLHSDes(nmeta, all.par.ranges)
X2 <- GetRandLHSDes(nmeta, all.par.ranges)

wrapper <- function(X){
  predict(m, newdata=X, type='UK')$m
}
system.time(sob <- sensitivity::sobol(model=wrapper, X1=X1, X2=X2, order=2, nboot=100))

vym <- sob$V['global', 'original']
vyd <- DiceKriging::coef(m)$sd2 + DiceKriging::coef(m)$nugget
vy <- vym + vyd

sob.inds <- sob$S[, 'original'] * vym / vy
names(sob.inds) <- rownames(sob$S)
sob.ind.rand <- vyd / vy

newdata <- head(X1, n=1000)
p <- predict(m, newdata=newdata, type='UK')

x <- newdata[, 'tprob.net']

plot(x, p$mean, ylim=c(0, 0.5))
points(x, p$lower95, col='grey')
points(x, p$upper95, col='grey')

points(sub$tprob.net, sub$r, col=2)

plotdes <- sensitivity::parameterSets(par.ranges=par.ranges[c('tprob.net',
                                          'tprob.sp')],
                                      samples=sqrt(nsim), method='grid')
des2 <- df[, names(all.par.ranges)]
des2[, 'tprob.net'] <- plotdes[1:nrow(des), 'tprob.net']
des2[, 'tprob.sp'] <- plotdes[1:nrow(des), 'tprob.sp']
pmean <- predict(m, newdata=des2, type='UK')$m

#lattice::levelplot(pmean~tprob.net*tprob.sp, data=as.data.frame(des2))
