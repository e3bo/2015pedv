#!/usr/bin/Rscript

data(flows.matrix, package='sds')

M <- round(flows.matrix)
rownames(M)[54] <- colnames(M)[54] <- 'Other states'
for(i in 1:6){
  inds <- seq((i - 1) * 9 + 1, i * 9)
  cap1 <- "USDA ERS estimates of number of swine moved between states in 2001"
  cap1B <- "for feeding and breeding. Origin in row, destination in column. Table"
  cap2 <- "of 6."
  cap <- paste(cap1, cap1B, i, cap2)
  xt <- xtable::xtable(M[ , inds], digits=0, caption=cap)
  file <- paste0('flows-tab-', i, '.tex')
  print(xt, file=file, caption.placement='top')
}
