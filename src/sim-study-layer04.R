#!/usr/bin/Rscript
library(methods) #for raster

set.seed(124, "L'Ecuyer")
options('mc.cores'=sds::GetCores())

load('sim-study-checkpoint3.rda')

kms <- sds::GetMetaModels(resall, df, var2="gcd")
save.image('sim-study-checkpoint4.rda')
saveRDS(kms$m2$model, "kmm2.rds")
saveRDS(kms$v2$model, "kmv2.rds")
saveRDS(kms$center, "center.rds")

tp.ind <- which(names(kms$center) == 'tprob.net')
tp.sp.ind <- which(names(kms$center) == 'tprob.sp')
sa.ind <- which(names(kms$center) == 'seasonal.amplitude')
pdf('km-m2-views.pdf')
DiceView::sectionview(kms$m2$model, axis=tp.ind, center=kms$center, mfrow=c(1,1))
DiceView::sectionview(kms$m2$model, axis=tp.sp.ind, center=kms$center, mfrow=c(1,1))
DiceView::sectionview(kms$m2, axis=tp.ind, center=kms$center, mfrow=c(1,1))
DiceView::contourview(kms$m2, axis=matrix(c(sa.ind, tp.sp.ind), nrow=1), center=kms$center, mfrow=c(1, 1))
dev.off()
