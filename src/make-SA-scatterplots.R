
test <- resall$season < .5 & resall$raster > 600
cols <- grep('mantel\\.[rp]\\.lag1\\..*spearman.TRUE\\.1000', colnames(resall))
sel <- data.frame(resall[, cols], tprob.net=resall$tprob.net)

m <- reshape2::melt(sel[test, ], id='tprob.net')
test <- grepl('mantel\\.r', m$variable)
m$Response <- ifelse(test, 'Spearman r', 'p value')
m$Matrix <- 'transport flows'
test <- grepl('gcd', m$variable)
m$Matrix[test] <- 'distance (km)'
test <- grepl('sharedBord', m$variable)
m$Matrix[test] <- 'shared border'

g <- ggplot2::ggplot(data=m, aes(x=tprob.net, y=value))
g <- g + geom_point(alpha=0.5)
g <- g + facet_grid(Response~Matrix, scales='free_y', labeller='label_both')
g <- g + xlab('Transmission probability for a transport contact')
g <- g + ylab('Response value')
g <- g + theme_classic()
ggplot2::ggsave('scatter-plots-facetted.pdf', width=18/2.54, height=6)
