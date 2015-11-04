#!/usr/bin/Rscript

agent.args <- list(target.mean.deg=1, census.dilation=0.1)
ag <- do.call(aporkalypse::CreateAgents, agent.args)

nsim <- 100
par.ranges <- list(nstarters=c(0, 10),
                   prep=c(0.01, 1),
                   size=c(0.7, 1),
                   rprob=c(0, 1),
                   seasonal.amplitude=c(0, 1),
                   starting.grid.x=c(0, 1),
                   starting.grid.y=c(0, 1),
                   tprob.outside=c(0, 0.1),
                   tprob.net=c(0, 0.1),
                   tprob.sp=c(0, 1))

des <- sensitivity::parameterSets(par.ranges=par.ranges,
                                  samples=nsim, method='sobol')
colnames(des) <- names(par.ranges)
sim.args <- data.frame(des, permutations=2)
df <- cbind(agent.args, sim.args)

system.time(res <- Map(aporkalypse::SimulateAndSummarize,
                       dynamic=list(ag),
                       lags.sel=1,
                       rprob=df$rprob,
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

formula <- as.formula(paste('~', paste(names(par.ranges), collapse='+')))
m <- DiceKriging::km(formula=formula, design=sub[, names(par.ranges)],
                     response=sub$r, nugget.estim=TRUE, covtype='matern3_2')

nmeta <- 1e4
foo <- function(y) runif(n=nmeta, min=y[1], max=y[2])
X1 <- sapply(par.ranges, foo)
X2 <- sapply(par.ranges, foo)

wrapper <- function(X){
  predict(m, newdata=X, type='UK')$m
}
sob <- sensitivity::sobol(model=wrapper, X1=X1, X2=X2, order=1, nboot=100)

vym <- sob$V['global', 'original']
vyd <- DiceKriging::coef(m)$sd2 + DiceKriging::coef(m)$nugget
vy <- vym + vyd

sob.inds <- sob$S[, 'original'] * vym / vy
names(sob.inds) <- rownames(sob$S)
sob.ind.rand <- vyd / vy

newdata <- head(X1, n=500)
p <- predict(m, newdata=newdata, type='UK')

x <- newdata[, 'tprob.net']
plot(x, p$mean, ylim=c(0, 0.5))
points(x, p$lower95, col='grey')
points(x, p$upper95, col='grey')

points(sub$tprob.net, sub$r, col=2)

test <- with(resall,
 symmetrize == TRUE & mat2.name == 'sharedBord' & method== 'spearman' & mat1.name == 'lag1')
sub2 <- resall[test, ]
plot(r~tprob.sp, data=sub2)
