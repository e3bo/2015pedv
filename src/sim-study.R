#!/usr/bin/Rscript

agent.args <- list(target.mean.deg=1, census.dilation=.1)
ag <- do.call(aporkalypse::CreateAgents, agent.args)

nsim <- 10
par.ranges <- list(tprob.net=c(0, 0.1), tprob.sp=c(0, 0.1),
                   rprob=c(0.1, 0.9))
des <- sensitivity::parameterSets(par.ranges=par.ranges,
                                  samples=nsim, method='sobol')
colnames(des) <- names(par.ranges)
sim.args <- data.frame(starting.state='OH', des, stringsAsFactors=FALSE,
                       tprob.outside=0.1)
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

m <- DiceKriging::km(~tprob.net + tprob.sp, design=sub[, names(par.ranges)],
                     response=sub$r, nugget.estim=TRUE, covtype='matern3_2')

nmeta <- 10000
X1 <- sensitivity::parameterSets(par.ranges=par.ranges,
                                 samples=nmeta, method='sobol')
X2 <- sensitivity::parameterSets(par.ranges=par.ranges,
                                 samples=nmeta, method='sobol')
colnames(X1) <- colnames(X2) <- names(par.ranges)

wrapper <- function(X){
  predict(m, newdata=X, type='UK')$m
}
sob <- sensitivity::sobol(model=wrapper, X1=X1, X2=X2, order=1, nboot=100)

vym <- sob$V['global', 'original']
vyd <- coef(m)$sd2 + coef(m)$nugget
vy <- vym + vyd

sob <- sob$S[, 'original'] * vym / vy
names(sob) <- rownames(x2$S)
sob.rand <- vyd / vy

newdata <- head(X1, n=500)
p <- predict(m, newdata=newdata, type='UK')

x <- newdata[, 3]
plot(x, p$mean, ylim=c(0, 0.5))
points(x, p$lower95, col='grey')
points(x, p$upper95, col='grey')
points(sub$rprob, sub$r, col=2)
