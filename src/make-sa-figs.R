#!/usr/bin/Rscript

library('ggplot2')

var2levs <- c('shipment', 'sharedBord', 'gcd')

GetSobInds <- function(var2, stat='r'){
  file <- paste0('sobol-indices-', stat, "-", var2, '.rds')
  readRDS(file)
}
S <- lapply(var2levs, GetSobInds)

data <- lapply(S, '[', c('original', 'min. c.i.', 'max. c.i.'))
Joiner <- function(x, y) cbind(Input=rownames(x), x, var2=y, stat='r')
data <- Map(Joiner, data, as.list(var2levs))
data <- do.call(rbind, data)

pal <- c("#4477AA", "#DDCC77", "#CC6677")
test <- data$`min. c.i.` > 0
dt <- data[test, ]

dt$pretty.var <- factor(dt$var2, levels=c('gcd', 'shipment', 'sharedBord'),
                        labels=c('distance (km)', 'transport\nflows', 'shared\nborder'))

levs <- names(sort(tapply(dt$original, dt$Input, mean, na.rm=TRUE)))
stopifnot(levs == c("seasonal.amplitude*tprob.net", "tprob.sp*raster.ncol", "prep",
              "seasonal.amplitude*raster.ncol", "tprob.sp", "seasonal.amplitude",
              "raster.ncol", "tprob.net"))
labs <- levs
labs <- gsub('seasonal.amplitude', 'Seasonal amplitude', labs)
labs <- gsub('tprob.net', 'P(trans. | transport edge)', labs)
labs <- gsub('tprob.sp', 'P(trans. | spatial edge)', labs)
labs <- gsub('raster.ncol', 'No. raster columns', labs)
labs <- gsub('prep', 'P(reporting | infection)', labs)
labs <- gsub('\\*', ' \\*\n  ', labs)
dt$Input <- factor(dt$Input, levels=levs, labels=labs)

g <- ggplot(data=dt,
            aes(x=Input, y=original, ymin=`min. c.i.`, ymax=`max. c.i.`))
g <- g + geom_pointrange(aes(colour=pretty.var, shape=pretty.var),
                         position=position_dodge(width=0.2))
g <- g + scale_shape(name='Matrix')
g <- g + scale_color_manual(name='Matrix', values=pal)
g <- g + coord_flip()
g <- g + theme_classic()
g <- g + theme(legend.position='top', plot.margin = grid::unit(c(-4, 1, 0, 0), "mm"))
g <- g + ylab("\nGlobal sensitivity index")

ggsave('sobol-indices.pdf', width=13.5/2.54, height=4)
ggsave('sobol-indices.eps', width=13.5/2.54, height=4)
