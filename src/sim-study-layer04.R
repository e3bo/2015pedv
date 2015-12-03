#!/usr/bin/Rscript
library(methods) #for raster

set.seed(124, "L'Ecuyer")
options('mc.cores'=1)

load('sim-study-checkpoint3.rda')

var2 <- Sys.getenv('var2')
if (nchar(var2) == 0) {
  var2 <- 'gcd'
}
kms <- sds::GetMetaModels(resall, df, var2=var2)
save.image(paste0('sim-study-checkpoint4-', var2, '.rda'))
saveRDS(kms$m2$model, paste0("kmm2-", var2, ".rds"))
saveRDS(kms$v2$model, paste0("kmv2-", var2, ".rds"))

tp.ind <- which(names(kms$center) == 'tprob.net')
tp.sp.ind <- which(names(kms$center) == 'tprob.sp')
sa.ind <- which(names(kms$center) == 'seasonal.amplitude')
file <- paste0('km-m2-views-', var2, '.pdf')
pdf(file)
DiceView::sectionview(kms$m2$model, axis=tp.ind, center=kms$center, mfrow=c(1,1))
DiceView::sectionview(kms$m2$model, axis=tp.sp.ind, center=kms$center, mfrow=c(1,1))
DiceView::sectionview(kms$m2, axis=tp.ind, center=kms$center, mfrow=c(1,1))
DiceView::contourview(kms$m2, axis=matrix(c(sa.ind, tp.sp.ind), nrow=1), center=kms$center, mfrow=c(1, 1))
dev.off()
