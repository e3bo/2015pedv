#!/usr/bin/Rscript
library(methods) #for raster

set.seed(125, "L'Ecuyer", "Inversion")
options('mc.cores'=1)

vars <- Sys.getenv('vars')
if (nchar(vars) == 0) {
  var2 <- 'gcd'
  stat <- 'r'
} else {
  vars <- strsplit(vars, split='-')[[1]]
  var2 <- vars[2]
  stat <- vars[1]
}

kmm2 <- readRDS(paste0("kmm2-", stat, "-", var2, ".rds"))
kmv2 <- readRDS(paste0("kmv2-", stat, "-", var2, ".rds"))
all.par.ranges <- readRDS("all.par.ranges.rds")
nmeta <- 1e5

sob.out <- sds::RunSobol(nmeta, kmm2=kmm2, kmv2=kmv2, all.par.ranges=all.par.ranges,
                         order=2)
save.image(paste0('sim-study-checkpoint5-', stat, "-", var2, '.rda'))

sob.out
saveRDS(sob.out[[1]]$S, paste0("sobol-indices-", stat, "-", var2, ".rds"))
