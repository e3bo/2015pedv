#!/usr/bin/Rscript

agent.args <- list(target.mean.deg=1, census.dilation=.1)
ag <- do.call(aporkalypse::CreateAgents, agent.args)

sim.args <- data.frame(starting.state='OH', tprob.sp=0, tprob.net=seq(0.0, 0.08, len=15),
                       rprob=0.5, stringsAsFactors=FALSE)
df <- cbind(agent.args, sim.args)

res <- Map(aporkalypse::SimulateAndSummarize, starting.state=df$starting.state,
           tprob.sp=df$tprob.sp, tprob.net=df$tprob.net, rprob=df$rprob,
           dynamic=list(ag), lags.sel=1, tprob.outside=.01)

tstats <- lapply(res, '[[', 'mantel.tests')
dfsplt <- split(df, seq(nrow(df)))
resplt <- Map(cbind, dfsplt, tstats)
resplt <- Filter(function(x) ncol(x) > 7, resplt)
resall <- do.call(rbind, resplt)

test <- with(resall, symmetrize == TRUE & mat2.name == 'shipment' & method== 'spearman' & mat1.name == 'lag1')
sub <- resall[test, ]

m <- DiceKriging::km(~tprob.net, design=sub[, 'tprob.net', drop=FALSE], response=sub$r,
                     nugget.estim=TRUE, covtype='matern3_2')

newdata <- data.frame(tprob.net=seq(0, to=.1, len=100))
p <- predict(m, newdata=newdata, type='UK')
x <- newdata[, 1]
plot(x, p$mean, type='l')
lines(x, p$lower95, col='grey')
lines(x, p$upper95, col='grey')
