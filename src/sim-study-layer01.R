#!/usr/bin/Rscript
library(methods) #for raster

set.seed(121, "L'Ecuyer")

options('mc.cores'=GetCores())

target.mean.deg.grid <- 10 #seq(0.1, 10.1, by=10)
raster.cell.side.grid <- 3200 # c(1600, 3200, 4800)

ag <- parallel::mcMap(sds::CreateAgents,
                      target.mean.deg=target.mean.deg.grid[1],
                      raster.cell.side.meters=raster.cell.side.grid,
                      census.dilation=1)

if(length(target.mean.deg.grid) > 1){
  net.nbs <- lapply(target.mean.deg.grid[-1], sds::GetNetNbs,
                    block.labels=ag[[1]]$adf$abb)
  net.nbs <- c(list(ag[[1]]$net.nbs), net.nbs)
} else {
  net.nbs <- list(ag[[1]]$net.nbs)
}

ag.data.inds <- expand.grid(ag.ind=seq_along(ag), net.nbs.ind=seq_along(net.nbs))

save.image('sim-study-checkpoint1.rda')
