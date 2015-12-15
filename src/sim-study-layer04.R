#!/usr/bin/Rscript
library(methods) #for raster
library(rgenoud) #for optim.method='gen' to work in GetMetaModels

options('mc.cores'=1)
load('sim-study-checkpoint3.rda')
set.seed(124, "L'Ecuyer", "Inversion")

vars <- Sys.getenv('vars')
if (nchar(vars) == 0) {
  var2 <- 'gcd'
  stat <- 'r'
} else {
  vars <- strsplit(vars, split='-')[[1]]
  var2 <- vars[2]
  stat <- vars[1]
}

kms <- sds::GetMetaModels(resall, df, var2=var2, statistic=stat, cortype="spearman")
save.image(paste0('sim-study-checkpoint4-', stat, '-', var2, '.rda'))
saveRDS(kms$m2$model, paste0("kmm2-", stat, "-", var2, ".rds"))
saveRDS(kms$var.m2, paste0("var.m2-", stat, "-", var2, ".rds"))

tp.ind <- which(names(kms$center) == 'tprob.net')
tp.sp.ind <- which(names(kms$center) == 'tprob.sp')
sa.ind <- which(names(kms$center) == 'seasonal.amplitude')
file <- paste0('km-m2-views-', stat, "-", var2, '.pdf')
pdf(file)
DiceView::sectionview(kms$m2$model, axis=tp.ind, center=kms$center, mfrow=c(1,1))
DiceView::sectionview(kms$m2$model, axis=tp.sp.ind, center=kms$center, mfrow=c(1,1))
DiceView::sectionview(kms$m2, axis=tp.ind, center=kms$center, mfrow=c(1,1))
DiceView::contourview(kms$m2, axis=matrix(c(sa.ind, tp.sp.ind), nrow=1), center=kms$center, mfrow=c(1, 1))
dev.off()
