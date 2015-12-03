#!/usr/bin/Rscript
library(methods) #for raster

set.seed(125, "L'Ecuyer")
options('mc.cores'=1)

var2 <- Sys.getenv('var2')
if (nchar(var2) == 0) {
  var2 <- 'gcd'
}
kmm2 <- readRDS(paste0("kmm2-", var2, ".rds"))
kmv2 <- readRDS(paste0("kmv2-", var2, ".rds"))
all.par.ranges <- readRDS("all.par.ranges.rds")
nmeta <- 1e3

sob.out <- sds::RunSobol(nmeta, kmm2=kmm2, kmv2=kmv2, all.par.ranges=all.par.ranges,
                         order=2)
save.image(paste0('sim-study-checkpoint5-', var2, '.rda'))

sob.out
