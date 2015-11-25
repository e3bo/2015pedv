#!/usr/bin/Rscript
library(methods) #for raster

set.seed(122, "L'Ecuyer")

mc.cores <- ifelse(Sys.info()['sysname'] == "Linux",
                   parallel::detectCores() - 1, 1)
mc.cores <- ifelse(mc.cores > 20, 20, mc.cores)
mc.cores <- ifelse(mc.cores == 0, 1, mc.cores)
options('mc.cores'=mc.cores)

load('sim-study-checkpoint1.rda')

nsim <- 1000
par.ranges <- list(prep=c(0, 1),
                   rprob=c(0, 1),
                   seasonal.amplitude=c(0, 1),
                   size=c(0.5, 10),
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
