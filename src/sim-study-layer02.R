#!/usr/bin/Rscript
library(methods) #for raster
options('mc.cores'=sds::GetCores())
load('sim-study-checkpoint1.rda')
set.seed(122, "L'Ecuyer", "Inversion")

nsim <- 1000
par.ranges <- list(seasonal.amplitude=c(0, 1),
		   starting.grid.x=c(0, 1),
                   starting.grid.y=c(0, 1),
                   tprob.outside=c(0, 1e-4),
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

extra.par.ranges <- list(target.mean.deg=range(target.mean.deg.grid),
                         raster.ncol=range(raster.ncol.grid),
                         raster.nrow=range(raster.nrow.grid))
all.par.ranges <- c(par.ranges, extra.par.ranges)

saveRDS(all.par.ranges, "all.par.ranges.rds")
