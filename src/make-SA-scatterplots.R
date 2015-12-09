#!/usr/bin/Rscript

library('ggplot2')
load('sim-study-checkpoint3.rda')

pal <- c("#4477AA", "#DDCC77", "#CC6677")
response.cols <- grep('mantel\\.[rp]\\.lag1\\..*spearman.TRUE\\.1000', colnames(resall))
input.cols <- c('tprob.net', 'tprob.sp', 'seasonal.amplitude', 'raster.ncol', 'prep')
sel <- cbind(resall[, input.cols], resall[, response.cols])


m <- reshape2::melt(sel, id=input.cols)
test <- grepl('mantel\\.r', m$variable)
m$Response <- ifelse(test, 'Spearman r', 'p value')
m$Matrix <- 'transport flows'
test <- grepl('gcd', m$variable)
m$Matrix[test] <- 'distance (km)'
test <- grepl('sharedBord', m$variable)
m$Matrix[test] <- 'shared border'

g <- ggplot(data=m, aes(x=tprob.net, y=value))
g <- g + geom_point(alpha=0.5, aes(colour=cut(seasonal.amplitude, breaks=c(0, 0.5, 1))))
g <- g + facet_grid(Response~Matrix, scales='free_y', labeller='label_both')
g <- g + xlab('\nTransmission probability for a transport contact')
g <- g + ylab('Response value\n')
g <- g + scale_x_continuous(labels=as.character)
g <- g + scale_color_manual(name='Seasonal\namplitude', values=pal)
g <- g + theme_classic()
ggsave('scatter-plots-facetted-tprob.net.pdf', plot=g, width=18/2.54, height=6)

g <- ggplot(data=m, aes(x=tprob.sp, y=value))
g <- g + geom_point(alpha=0.5, aes(colour=as.factor(raster.ncol)))
g <- g + facet_grid(Response~Matrix, scales='free_y', labeller='label_both')
g <- g + xlab('\nTransmission probability for a spatial contact')
g <- g + ylab('Response value\n')
g <- g + scale_x_continuous(labels=as.character)
g <- g + scale_color_manual(name='Raster\ncolumns', values=pal)
g <- g + theme_classic()
ggsave('scatter-plots-facetted-tprob.sp.pdf', plot=g, width=18/2.54, height=6)

g <- ggplot(data=m, aes(x=as.factor(raster.ncol), y=value))
g <- g + geom_jitter(alpha=0.5, aes(colour=cut(seasonal.amplitude, breaks=c(0, 0.5, 1))))
g <- g + facet_grid(Response~Matrix, scales='free_y', labeller='label_both')
g <- g + xlab('Number of columns in raster')
g <- g + ylab('Response value')
g <- g + scale_color_manual(name='Seasonal\namplitude', values=pal)
g <- g + theme_classic()
ggsave('scatter-plots-facetted-raster.ncol.sp.pdf', plot=g, width=18/2.54, height=6)

g <- ggplot(data=m, aes(x=as.factor(raster.ncol), y=value))
g <- g + geom_jitter(alpha=0.5, aes(colour=cut(seasonal.amplitude, breaks=c(0, 0.5, 1))))
g <- g + facet_grid(Response~Matrix, scales='free_y', labeller='label_both')
g <- g + xlab('Number of columns in raster')
g <- g + ylab('Response value')
g <- g + scale_color_manual(name='Seasonal\namplitude', values=pal)
g <- g + theme_classic()
ggsave('scatter-plots-facetted-raster.ncol.sp.pdf', plot=g, width=18/2.54, height=6)

g <- ggplot(data=m, aes(x=prep, y=value))
g <- g + geom_point(alpha=0.5)
g <- g + facet_grid(Response~Matrix, scales='free_y', labeller='label_both')
g <- g + xlab('\nProbability of reporting')
g <- g + ylab('Response value\n')
g <- g + scale_x_continuous(labels=as.character)
g <- g + theme_classic()
ggsave('scatter-plots-facetted-prep.pdf', plot=g, width=18/2.54, height=6)
