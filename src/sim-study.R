#!/usr/bin/Rscript

agent.args <- list(target.mean.deg=1, census.dilation=.1)
ag <- do.call(aporkalypse::CreateAgents, agent.args)

nsim <- 100
par.ranges <- list(tprob.net=c(0, 0.1), tprob.sp=c(0, 0.1),
                   rprob=c(0, 1), tprob.outside=c(0, 0.1))
des <- sensitivity::parameterSets(par.ranges=par.ranges,
                                  samples=nsim, method='sobol')
colnames(des) <- names(par.ranges)
sim.args <- data.frame(starting.state='OH', des, stringsAsFactors=FALSE,
                       permutations=2)
df <- cbind(agent.args, sim.args)

system.time(res <- Map(aporkalypse::SimulateAndSummarize,
                       starting.state=df$starting.state,
                       tprob.sp=df$tprob.sp, tprob.net=df$tprob.net,
                       rprob=df$rprob, tprob.outside=df$tprob.outside,
                       dynamic=list(ag), lags.sel=1))

tstats <- lapply(res, '[[', 'mantel.tests')
dfsplt <- split(df, seq(nrow(df)))
resplt <- Map(cbind, dfsplt, tstats)
resplt <- Filter(function(x) ncol(x) > 7, resplt)
resall <- do.call(rbind, resplt)

test <- with(resall,
 symmetrize == TRUE & mat2.name == 'shipment' & method== 'spearman' & mat1.name == 'lag1')
sub <- resall[test, ]

m <- DiceKriging::km(~tprob.net + tprob.sp + rprob, design=sub[, names(par.ranges)],
                     response=sub$r, nugget.estim=TRUE, covtype='matern3_2')

nmeta <- 10000
foo <- function(y) runif(n=nmeta, min=y[1], max=y[2])
X1 <- sapply(par.ranges, foo)
X2 <- sapply(par.ranges, foo)

wrapper <- function(X){
  predict(m, newdata=X, type='UK')$m
}
sob <- sensitivity::sobol(model=wrapper, X1=X1, X2=X2, order=2, nboot=1000)

vym <- sob$V['global', 'original']
vyd <- DiceKriging::coef(m)$sd2 + DiceKriging::coef(m)$nugget
vy <- vym + vyd

sob.inds <- sob$S[, 'original'] * vym / vy
names(sob.inds) <- rownames(sob$S)
sob.ind.rand <- vyd / vy

newdata <- head(X1, n=500)
p <- predict(m, newdata=newdata, type='UK')

x <- newdata[, 1]
plot(x, p$mean, ylim=c(0, 0.5))
points(x, p$lower95, col='grey')
points(x, p$upper95, col='grey')

points(sub$rprob, sub$r, col=2)
