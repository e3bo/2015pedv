#!/usr/bin/Rscript

agent.args <- list(target.mean.deg=1, census.dilation=.1)
ag <- do.call(aporkalypse::CreateAgents, agent.args)

samp.lhs <- lhs::optimumLHS(100, 3)
sim.args <- data.frame(starting.state='OH',
                       tprob.sp=samp.lhs[, 1] * 0.1,
                       tprob.outside=samp.lhs[, 2] * 0.01,
                       tprob.net=samp.lhs[, 3] * 0.1,
                       rprob=0.5, stringsAsFactors=FALSE)
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

test <- with(resall, symmetrize == TRUE & mat2.name == 'shipment' & method== 'spearman' & mat1.name == 'lag1')
sub <- resall[test, ]

m <- DiceKriging::km(~tprob.net + tprob.sp, design=sub[, c('tprob.net', 'tprob.sp')],
                     response=sub$r, nugget.estim=TRUE, covtype='matern3_2')

n <- 10000
X1 <- data.frame(matrix(runif(2 * n, max=0.1), nrow=n))
X2 <- data.frame(matrix(runif(2 * n, max=0.1), nrow=n))
colnames(X1) <- colnames(X2) <- c('tprob.net', 'tprob.sp')

wrapper <- function(X){
  predict(m, newdata=X, type='UK')$m
}
x2 <- sensitivity::sobol(model= wrapper, X1 = X1, X2=X2, order=2, nboot=100)

vym <- x2$V['global', 'original']
vyd <- coef(m)$sd2 + coef(m)$nugget
vy <- vym + vyd

sob <- x2$S[, 'original'] * vym / vy
names(sob) <- rownames(x2$S)
sob.rand <- vyd / vy

newdata <- head(X1, n=500)
p <- predict(m, newdata=newdata, type='UK')
x <- newdata[, 1]
plot(x, p$mean, ylim=c(0, 0.5))
points(x, p$lower95, col='grey')
points(x, p$upper95, col='grey')
points(sub$tprob.net, sub$r, col=2)
