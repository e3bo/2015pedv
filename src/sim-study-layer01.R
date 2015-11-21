#!/usr/bin/Rscript
library(methods) #for raster

set.seed(122, "L'Ecuyer")

mc.cores <- ifelse(Sys.info()['sysname'] == "Linux",
                   parallel::detectCores() - 1, 1)
mc.cores <- ifelse(mc.cores > 20, 20, mc.cores)
mc.cores <- ifelse(mc.cores == 0, 1, mc.cores)
options('mc.cores'=mc.cores)

target.mean.deg.grid <- 0.1 #seq(0.1, 10.1, by=10)
raster.cell.side.grid <- 16000 # c(1600, 3200, 4800)

ag <- parallel::mcMap(sds::CreateAgents,
                      target.mean.deg=target.mean.deg.grid[1],
                      raster.cell.side.meters=raster.cell.side.grid,
                      census.dilation=0.01)

if(length(target.mean.deg.grid) > 1){
  net.nbs <- lapply(target.mean.deg.grid[-1], sds::GetNetNbs,
                    block.labels=ag[[1]]$adf$abb)
  net.nbs <- c(list(ag[[1]]$net.nbs), net.nbs)
} else {
  net.nbs <- list(ag[[1]]$net.nbs)
}

ag.data.inds <- expand.grid(ag.ind=seq_along(ag), net.nbs.ind=seq_along(net.nbs))
save.image('sim-study-checkpoint1.rda')

nsim <- 1000
par.ranges <- list(prep=c(0, 1),
                   rprob=c(0, 1),
                   seasonal.amplitude=c(0, 1),
                   size=c(0.5, 100),
                   tprob.outside=c(0, 1e-2),
                   tprob.net=c(0, 1),
                   tprob.sp=c(0, 1))

des <- sensitivity::parameterSets(par.ranges=par.ranges,
                                  samples=nsim, method='sobol')
colnames(des) <- names(par.ranges)

desunif <- matrix(runif(nsim * ncol(des)), nrow=nsim)

GetDesignStats <- function(X){
    list(coverage=DiceDesign::coverage(X),
         meshRatio=DiceDesign::meshRatio(X),
         mindist=DiceDesign::mindist(X),
         discrepancyCriteriaL2=DiceDesign::discrepancyCriteria(X, type='L2')$DisL2)
}
(design.stats <- sapply(list(sobol=des, unif=desunif), GetDesignStats))
save.image('sim-study-checkpoint2.rda')

df <- do.call(cbind, c(list(des), as.list(ag.data.inds)))
df <- as.data.frame(df)

df$raster.cell.side <- raster.cell.side.grid[df$ag.ind]
df$target.mean.deg <- target.mean.deg.grid[df$net.nbs.ind]
df$lags.sel <- 1
df$nstarters <- 1
df$permutations <- 2
Wrapper <- function(...) try(sds::SimulateAndSummarize(...))

system.time(res <- parallel::mcMap(Wrapper,
                                   agent.data=ag[df$ag.ind],
                                   lags.sel=df$lags.sel,
                                   net.nbs=net.nbs[df$net.nbs.ind],
                                   nstarters=df$nstarters,
                                   permutations=df$permutations,
                                   prep=df$prep,
                                   rprob=df$rprob,
                                   size=df$size,
                                   seasonal.amplitude=df$seasonal.amplitude,
                                   tprob.net=df$tprob.net,
                                   tprob.outside=df$tprob.outside,
                                   tprob.sp=df$tprob.sp))
recs <- lapply(res, '[[', 'record')
recs <- do.call(rbind, recs)
resall <- cbind(df, recs)
save.image('sim-study-checkpoint3.rda')
