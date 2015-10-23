#!/usr/bin/Rscript

agent.args <- list(target.mean.deg=1, census.dilation=.1)
ag <- do.call(aporkalypse::CreateAgents, agent.args)

sim.args <- data.frame(rep=1:2, starting.state='IA', tprob.sp=1, tprob.net=1,
                       rprob=0.5, stringsAsFactors=FALSE)
df <- cbind(agent.args, sim.args)

res <- Map(aporkalypse::SimulateAndSummarize, starting.state=df$starting.state,
           tprob.sp=df$tprob.sp, tprob.net=df$tprob.net, rprob=df$rprob,
           dynamic=list(ag))

tstats <- lapply(res, '[[', 'mantel.tests')
dfsplt <- split(df, seq(nrow(df)))
resplt <- Map(cbind, dfsplt, tstats)
do.call(rbind, resplt)

