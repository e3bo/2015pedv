#!/usr/bin/Rscript

agent.args <- list(target.mean.deg=1, census.dilation=.1)
ag <- do.call(aporkalypse::CreateAgents, agent.args)

sim.args <- data.frame(starting.state='NC', tprob.sp=1, tprob.net=seq(0, 0.4, len=10),
                       rprob=0.5, stringsAsFactors=FALSE)
df <- cbind(agent.args, sim.args)

res <- Map(aporkalypse::SimulateAndSummarize, starting.state=df$starting.state,
           tprob.sp=df$tprob.sp, tprob.net=df$tprob.net, rprob=df$rprob,
           dynamic=list(ag))

tstats <- lapply(res, '[[', 'mantel.tests')
dfsplt <- split(df, seq(nrow(df)))
resplt <- Map(cbind, dfsplt, tstats)
resall <- do.call(rbind, resplt)

test <- resall$symmetrize == TRUE & resall$mat2.name == 'shipment' & resall$method== 'spearman'
sub <- resall[test, ]


fit <- kernlab::gausspr(r~tprob.net, data=sub, var=2, kernel='laplacedot')
plot(r~tprob.net, data=sub)
x <- sub[, 'tprob.net', drop=FALSE]
pred <- predict(fit, x)
lines(cbind(x, pred))

plot(pValues~tprob.net, data=sub)
