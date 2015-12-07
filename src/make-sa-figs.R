#!/usr/bin/Rscript

library(ggplot2)

var2levs <- c('shipment', 'gcd', 'sharedBord')

GetSobInds <- function(var2, stat='r'){
  file <- paste0('sim-study-checkpoint5-', var2, '.rda')
  load(file)
  sob.out[[1]]$S
}
S <- lapply(var2levs, GetSobInds)

data <- lapply(S, '[', c('original', 'min. c.i.', 'max. c.i.'))
Joiner <- function(x, y) cbind(Input=rownames(x), x, var2=y, stat='r')
data <- Map(Joiner, data, as.list(var2levs))
data <- do.call(rbind, data)

test <- data$original > 0.05
dt <- data[test, ]
g <- ggplot(data=dt,
            aes(x=Input, y=original, ymin=`min. c.i.`, ymax=`max. c.i.`))
g <- g + geom_pointrange(aes(colour=var2, shape=var2),
                        position=position_dodge(width=0.2))
g <- g + coord_flip()
ggsave('sobol-indices.pdf')
