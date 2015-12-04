#!/usr/bin/Rscript
library(methods) #for raster

set.seed(121, "L'Ecuyer")

mc.cores <- ifelse(Sys.info()['sysname'] == "Linux",
                   parallel::detectCores() - 1, 1)
mc.cores <- ifelse(mc.cores > 20, 20, mc.cores)
mc.cores <- ifelse(mc.cores == 0, 1, mc.cores)
options('mc.cores'=mc.cores)

target.mean.deg.grid <- 10 #seq(0.1, 10.1, by=10)
raster.ncol.grid <- 2850 / 2 # 285 x 178 corresponds to roughly 16 km^2 square cells
raster.nrow.grid <- 1780 / 2

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
