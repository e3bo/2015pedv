#!/usr/bin/Rscript

library(fields)
library(GGally)
library(ggplot2)
library(Hmisc)
library(igraph)
library(mapproj)
library(maps)
library(pander)
library(plyr)
library(RColorBrewer)
library(reshape2)
library(scales)
source('custom-plot-functions.r')

load(file='mantel-testing-checkpoint1.RData')

## plotting

### jflowmap
projection <- 'azequalarea'
orientation <- c(30.82,-98.57,0)
#scale of 1 does not allow good edge-bundling with this projection
scale <- 100

adjmatrix <- data.matrix(epO[state.abb, state.abb])
diag(adjmatrix) <- 0
g <- graph.adjacency(adjmatrix, mode='directed', weighted='flow')
stopifnot(colnames(adjmatrix) == V(g)$name)
V(g)$PEDVpositive <- V(g)$name %in% colnames(CC)
key <- match(V(g)$name, state.abb)
proj <- mapproject(x=state.center$x[key], y=state.center$y[key],
                   projection=projection, orientation=orientation)
V(g)$x <- proj$x*scale
V(g)$y <- proj$y*scale
adjmatrix <- getMat('cor')
diag(adjmatrix) <- 0
g2 <- graph.adjacency(adjmatrix, mode='directed', weighted='crossCorrelation')
E(g2)$ccShifted <- E(g2)$crossCorrelation + 1
g3 <- graph.union(g, g2)
adjmatrix <- (adjmatrix + t(adjmatrix))/2
g2 <- graph.adjacency(adjmatrix, mode='undirected', weighted='symCC')
#g3 <- graph.union(g3, g2)
#V(g3)$name <- paste('US', V(g3)$name, sep='-')
g3 <- delete.vertices(g3, c('AK', 'HI'))
write.graph(g3, file='swine-flows-CC.xml', format='graphml')

scl <- 1.8
png(file='us-background.png', width=535*scl, height=283*scl)
par(mai=rep(0,4))
map('state', mar=rep(0,4), projection=projection, orientation=orientation,
    resolution=0, myborder=0, col='dark grey', lwd=3)
bbox <- par('usr')*scale
# uncomment to verify bbox
#do.call(rect, c(as.list(bbox[c(1,3,2,4)]/scale), lwd=30))
dev.off()

# bbox as read by jflowmap
paste(bbox[1], -bbox[4], bbox[2]-bbox[1], bbox[4]-bbox[3], sep=',')

## choropleths

theme_clean <- function (base_size=12) {
  theme_grey(base_size) %+replace% theme(
      axis.title = element_blank (),
      axis.text = element_blank (),
      panel.background = element_blank (),
      panel.grid = element_blank (),
      axis.ticks.length = unit (0,"cm"),
      axis.ticks.margin = unit (0,"cm"),
      panel.margin = unit (0,"lines"),
      plot.margin = unit(c(0,0,0,0),"lines"),
      complete =TRUE)
}

state_map <- map_data("state")

## use data from cumulative burden analysis to ensure consistency

get_map_data <- function(){
    load('pedv-cum.RData')
    key <- match(mg$State, state.abb)
    ret <- data.frame(label=mg$State, state=tolower(state.name[key]),
                           Cases=mg$cases, Inventory=mg$inventory2012,
                           x=mg$stateLong, y=mg$stateLat)
    ret$Cases <- ifelse(ret$Cases == 0, NA, ret$Cases)
    ret
}

map_data <- get_map_data()
## manually position some state labels
i <- which(map_data$label=='NH')
dy <- 1.4
map_data$x[i] <- -68.08
map_data$y[i] <- map_data$y[map_data$label=='VT'] - 1.3
ii <- which(map_data$label=='MA')
map_data$x[ii] <- map_data$x[i]
map_data$y[ii] <- map_data$y[i] - dy
ii <- which(map_data$label=='RI')
map_data$x[ii] <- map_data$x[i] - 0.3
map_data$y[ii] <- map_data$y[i] - 2*dy
i <- which(map_data$label=='DE')
map_data$x[i] <- map_data$x[i] + 2
map_data$y[i] <- map_data$y[i] - 1
ii <- which(map_data$label=='NJ')
map_data$x[ii] <- map_data$x[i] + 1.2
map_data$y[ii] <- map_data$y[ii] - 1.2
i <- which(map_data$label=='CT')
map_data$x[i] <- map_data$x[ii] + 0.5
map_data$y[i] <- map_data$y[i] - 1.7
i <- which(map_data$label=='MD')
map_data$x[i] <- map_data$x[ii] - 1.8
map_data$y[i] <- map_data$y[i] - 3

get_endpoints <- function(state){
    i1 <- which(map_data$label==state)
    i2 <- which(state.abb==state)
    data.frame(x=state.center$x[i2], y=state.center$y[i2],
               xend=map_data$x[i1] - 1.5, yend=map_data$y[i1] + 0.3)
}

g <- make_choropleth(fill_var='Cases')
ggsave(filename='cases-choropleth.pdf', plot=g, width=8.6/2.54, height=8.6/2.54)
g <- make_choropleth('Inventory')
ggsave(filename='inventory-choropleth.pdf', plot=g, width=8.6/2.54, height=8.6/2.54)

### ggPairs

mats2 <- mats[1:3]
getData <- function(matName, sym, rankv){
    dists <- getMat(matName)
    if(sym){
        dists <- (dists + t(dists))/2
    }
    ret <- as.vector(as.dist(dists))
    if(rankv){
        ret <- rank(ret)
    }
    return(ret)
}
tmpf <- function(x, ...){
    sapply(mats2, getData, sym=x, ...)
}
tmpff <- function(x) {
    lapply(c(undirected=TRUE, directed=FALSE), tmpf, rankv=x)
}
matData <- lapply(c(ranked=TRUE, original=FALSE), tmpff)

theme_set(theme_classic())
theme_update(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             axis.ticks=element_blank(), panel.border=element_blank(),
             axis.line=element_blank())


pm <- makePlotMat(dfl=matData, type='directed', transform='ranked')
pdf(file='ggpairs-spearman-directed.pdf', width=8.6/2.54, height=8.6/2.54)
print(pm)
dev.off()

pm <- makePlotMat(dfl=matData, type='undirected', transform='ranked')
pdf(file='ggpairs-spearman-undirected.pdf', width=8.6/2.54, height=8.6/2.54)
print(pm)
dev.off()

pm <- makePlotMat(dfl=matData, type='directed', transform='original',
                  labelBreaks=TRUE, gridLabelSize=2.0)
plotMatDirectedPearson <- pm
pdf(file='ggpairs-pearson-directed.pdf', width=9.5/2.54, height=6.6/2.54)
print(pm)
dev.off()

pm <- makePlotMat(dfl=matData, type='undirected', transform='original',
                  labelBreaks=TRUE)
pdf(file='ggpairs-pearson-undirected.pdf', width=11.4/2.54, height=1.2*8.6/2.54)
print(pm)
dev.off()

### diagnostics

df <- data.frame(epl=unclass(as.dist(epl)), CC=unclass(as.dist(CC)))
m <- lm(CC~epl, data=df)
pdf('diagnostics.pdf')
plot(m)
dev.off()

## image plots

pdf('cors-heatmap.pdf')
orderings <- heatmap(data.matrix(epl), scale='none')
dev.off()

makeImagePlot <- function(M, orderings, ...){
    x <- data.matrix(M)[orderings$rowInd, orderings$colInd]
    x <- t(x)
    col <- brewer.pal(9, 'Blues')
    nc <- ncol(x)
    nr <- nrow(x)
    xi <- 1:nr
    yi <- 1:nc
    names(xi) <- colnames(x)
    names(yi) <- rownames(x)
    image.plot.ebo(x=xi, y=yi, z=x, horizontal=FALSE,col=col, graphics.reset=TRUE, ...)
}


pdf('matrices.pdf', width=8.7/2.54, height=12/2.54)
layout(matrix(1:2, ncol=1))
par(mar=c(4,4.5,.5,.5))
makeImagePlot(M=data.matrix(epl), orderings=orderings, xlab='Destination', ylab='Source',
              legend.args=list(text=expression(paste(Log[10], '(#swine / y)')),
                  line=2.9, side=4), panelLab='a')
par(mar=c(4,4.5,1,.5))
CCnoDiag <- CC
diag(CCnoDiag) <- NA
makeImagePlot(M=data.matrix(CCnoDiag), orderings=orderings, xlab='Leading state', ylab='Lagging state',
              legend.args=list(text='Cross correlation', line=2.9, side=4), panelLab='b')
dev.off()

### time series

mts <- melt(caseData, id='week')
keepers <- c('MN', 'KS',
             'IL', 'OK',
             'IA', 'NC')
test <- mts$variable %in% keepers
mts <- mts[test, ]
mts$variable <- factor(mts$variable, levels=keepers)
mts$x <- as.Date(as.character(mts$week), format='%m/%d/%Y')
labdf <- ddply(mts, 'variable', summarize, minx=x[1], maxy=max(value))
labdf$minx[labdf$variable == 'OK'] <- as.Date('2013-11-01')

tmpf <- function(df){
    g <- ggplot(df, aes(x=x, y=value, group=variable))
    g <- g + geom_hline(yintercept=0, size=0.25, col='grey')
    g <- g + geom_step(direction="vh", size=0.25)
    g <- g + facet_wrap(~variable, ncol=2, scales='free_y')
    g <- g + scale_x_date()
    g <- g + scale_y_discrete(breaks=pretty_breaks(n=2))
    g <- g + xlab('2013-2014 Date') + ylab('Positive accessions')
    g <- g + theme_classic()
    g <- g + theme(strip.background = element_blank(),
                   strip.text.x = element_blank())
    g <- g + theme(plot.margin=unit(c(0,2,2,0),"mm"))
    g <- g + geom_text(data=labdf, hjust=0, vjust=1, size=10 * 0.35,
                       aes(x=minx, y=maxy, label=variable))
    g <- g + theme(axis.title.x = element_text(vjust=-0.5))
    g <- g + theme(text = element_text(size=10))
    g <- g + theme(line = element_line(size=0.25))
}

g <- tmpf(mts)


ggsave('ts.pdf', width=8.6/2.54, height=6/2.54)

### composite time series and scatterplot matrix

tmpf <- function(){
    grid.newpage()
    lay <- grid.layout(1,2, widths=unit(c(8.6,8.6), 'cm'),
                       heights=unit(c(6), 'cm'))
    pushViewport(viewport(layout=lay))
    pushViewport(viewport(layout.pos.col=1, layout.pos.row=1))
    myPrintGGpairs(plotMatDirectedPearson, newpage=FALSE)
    grid.text(label="a", x=unit(0, "npc") - unit(1.5, "lines"),
              y=unit(1, "npc"), just= "left", gp=gpar(fontsize=10))
    popViewport()
    pushViewport(viewport(layout.pos.col=2, layout.pos.row=1))
    print(g, newpage=FALSE)
    grid.text(label="b", x=unit(0, "npc"), y=unit(1, "npc"),
              just= "left",  gp=gpar(fontsize=10))
    popViewport()
}

cairo_ps(filename = 'plotMatrixWithTimeSeries.eps', width = 19/2.54, height = 6.4/2.54)
tmpf()
dev.off()

cairo_pdf(filename = 'plotMatrixWithTimeSeries.pdf', width = 19/2.54, height = 6.4/2.54)
tmpf()
dev.off()
