#!/usr/bin/Rscript
library(methods) #for raster

set.seed(125, "L'Ecuyer")
options('mc.cores'=sds::GetCores())

kmm2 <- readRDS("kmm2.rds")
kmv2 <- readRDS("kmv2.rds")
center <- readRDS("center.rds")
all.par.ranges <- readRDS("all.par.ranges.rds")
nmeta <- 1e4

sob.out <- sds::RunSobol(nmeta, kmm2=kmm2, kmv2=kmv2, all.par.ranges=all.par.ranges,
                         order=2)
save.image('sim-study-checkpoint5.rda')

sob.out
