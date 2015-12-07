#!/usr/bin/Rscript
library(methods) #for raster

set.seed(121, "L'Ecuyer")

options('mc.cores'=sds::GetCores())

target.mean.deg.grid <- 1:10
raster.ncol.grid <- c(1425, 570) # 285 x 178 corresponds to roughly (16 km)^2 square cells
raster.nrow.grid <- c(890, 356)

ag <- parallel::mcMap(sds::CreateAgents,
                      target.mean.deg=target.mean.deg.grid[1],
                      raster.ncol=raster.ncol.grid,
                      raster.nrow=raster.nrow.grid,
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
