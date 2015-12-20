#!/usr/bin/Rscript
library(methods) #for raster

options('mc.cores'=sds::GetCores())
load('sim-study-checkpoint2.rda')
set.seed(1233, "L'Ecuyer-CMRG", "Inversion")

df <- do.call(cbind, c(list(des), as.list(ag.data.inds)))
df <- as.data.frame(df)

df$prep <- 0.1

df$raster.ncol <- raster.ncol.grid[df$ag.ind]
df$rprob <- 0.5
                                        #df$raster.nrow <- raster.nrow.grid[df$ag.ind] , ## commented out because don't expect to vary nrow independently of ncol
df$size <- 1.75
df$target.mean.deg <- target.mean.deg.grid[df$net.nbs.ind]
df$lags.sel <- 1
df$nstarters <- 1
df$permutations <- 1000
df$starting.grid.nx <- 10
df$starting.grid.ny <- 2
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
				   starting.grid.nx=df$starting.grid.nx,
				   starting.grid.ny=df$starting.grid.ny,
				   starting.grid.x=df$starting.grid.x,
				   starting.grid.y=df$starting.grid.y,
				   tprob.net=df$tprob.net,
                                   tprob.outside=df$tprob.outside,
                                   tprob.sp=df$tprob.sp))
recs <- lapply(res, '[[', 'record')
recs <- do.call(rbind, recs)
resall <- cbind(df, recs)
save.image('sim-study-checkpoint3.rda')
