#!/usr/bin/Rscript

library(fields)
library(GGally)
library(ggplot2)
library(igraph)
library(mapproj)
library(maps)
library(pander)
library(plyr)
library(RColorBrewer)
library(reshape2)
library(scales)
library(vegan)
source('custom-plot-functions.r')

Sys.setlocale("LC_TIME", "C") #Needed for identical()
Sys.setlocale("LC_COLLATE", "C")

## Get case data

dataDir <- '.'

tmpf <- function(){
  fn <- file.path(dataDir, 'PEDvweeklyreport-state-ts-01-08-14.csv')
  ret <- read.csv(fn)
  ret[1:9, c('CA', 'MD', 'NE', 'WY')] <- 0
  key <- order(colnames(ret)[-c(1:2)])
  ret <- cbind(ret[, c(1:2)], ret[, -c(1:2)][, key])
  ret$Unk <- NULL
  target <- structure(list(week = structure(c(20L, 29L, 32L, 8L, 12L),
.Label = c("10/13/2013", "10/20/2013", "10/27/2013", "10/6/2013",
"11/10/2013", "11/17/2013", "11/24/2013", "11/3/2013", "12/1/2013",
"12/15/2013", "12/22/2013", "12/29/2013", "12/8/2013", "4/15/2013",
"4/22/2013", "4/29/2013", "5/13/2013", "5/20/2013", "5/27/2013",
"5/6/2013", "6/10/2013", "6/16/2013", "6/23/2013", "6/3/2013",
"6/30/2013", "7/14/2013", "7/21/2013", "7/28/2013", "7/7/2013",
"8/11/2013", "8/18/2013", "8/25/2013", "8/4/2013", "9/1/2013",
"9/15/2013", "9/22/2013", "9/29/2013", "9/8/2013"), class = "factor"),
totalNumberSwineAccessions = c(17L, 34L, 26L, 90L, 134L), CA = c(0, 0,
0, 0, 1), CO = c(1L, 1L, 0L, 0L, 0L), IA = c(8L, 6L, 2L, 38L, 54L), IL
= c(0L, 1L, 0L, 1L, 14L), IN = c(3L, 2L, 1L, 1L, 3L), KS = c(0L, 4L,
3L, 6L, 4L), KY = c(0L, 0L, 0L, 0L, 0L), MD = c(0, 0, 0, 0, 0), MI =
c(0L, 0L, 0L, 2L, 0L), MN = c(1L, 2L, 2L, 7L, 20L), MO = c(0L, 0L, 0L,
2L, 4L), NC = c(0L, 3L, 4L, 14L, 18L), NE = c(0, 0, 0, 0, 2), NY =
c(0L, 0L, 0L, 0L, 0L), OH = c(0L, 2L, 1L, 5L, 5L), OK = c(0L, 11L,
10L, 10L, 2L), PA = c(1L, 0L, 3L, 1L, 0L), SD = c(0L, 0L, 0L, 0L, 3L),
TN = c(0L, 0L, 0L, 1L, 1L), TX = c(0L, 0L, 0L, 1L, 1L), WI = c(0L, 0L,
0L, 1L, 0L ), WY = c(0, 0, 0, 0, 1)), .Names = c("week",
"totalNumberSwineAccessions", "CA", "CO", "IA", "IL", "IN", "KS",
"KY", "MD", "MI", "MN", "MO", "NC", "NE", "NY", "OH", "OK", "PA",
"SD", "TN", "TX", "WI", "WY" ), row.names = c(4L, 13L, 20L, 30L, 38L),
class = "data.frame")
  stopifnot(identical(target, ret[c(4,13,20,30,38),]))
  ret
}
caseData <- tmpf()

unwanted <- c('week', 'totalNumberSwineAccessions', 'Unk')
ind <- which(!colnames(caseData) %in% unwanted)
observed <- caseData[, ind]

## Get cross correlations

n <- ncol(observed)
CC <- matrix(nrow=n, ncol=n)

getcc <- function(x,y, lag=1){
    foo <- ccf(x, y, plot=FALSE)
    ind <- which(foo$lag == lag)
    foo$acf[ind]
}

for(i in seq_len(n)){
    for(j in seq_len(n)){
        ## CC[i,j] will be high if deviations from the mean in series i
        ## are shifted to the left of similar deviations to the mean in series j
        ## i.e. i's deviations are indicative of j's future deviations
        CC[i,j] <- getcc(observed[,j], observed[,i])
    }
}

colnames(CC) <- rownames(CC) <- colnames(observed)

check_cc_directionality <- function(){
    t <- 1:100 * .25
    x <- sin(t)
    ## y is shifted to the right
    y <- sin(t-1)

    xylag1 <- getcc(x, y)
    yxlag1 <- getcc(y,x)
    main <- paste('getcc(x, y) ==', round(xylag1,2),
                  '; getcc(y,x) == ', round(yxlag1,2))
    if(xylag1 > yxlag1){
        conc <- 'left arg delayed by lag'
    } else{
        conc <- 'right arg delayed by lag'
    }
    plot(x~t, type='l', ylab='f(t)', main=main, sub=conc)
    lines(y~t, col=2)
    legend('topright', col=1:2, legend=c('x', 'y'), lty=1)
}

png('direction-check.png')
check_cc_directionality()
dev.off()

## Get shared-border neighborhoods

nbEdgelist <- read.csv(file.path(dataDir, 'state_neighbors_fips.txt'), header=FALSE)
data(state.fips)
g <- graph.data.frame(nbEdgelist, directed=FALSE)
g <- simplify(g)
key <- match(V(g)$name, state.fips$fips)
abb <- state.fips$abb[key]
V(g)$name <- as.character(abb)
vids <- which(V(g)$name %in% colnames(caseData))
g2 <- induced.subgraph(g, vids)
nms <- colnames(observed)
nhood <- as.matrix(g2[nms,nms])

## Get shipment flows

epO <- ep <- read.csv(file.path(dataDir, 'shipment-flows-origins-on-rows-dests-on-columns.csv'),row.names=1)
key <- match(colnames(observed), colnames(ep))
ep <- ep[key, key]
rownames(ep) <- colnames(ep)

### assertions based on inspection of original .xls file
stopifnot(ep['CA', 'IL'] == 9415)
stopifnot(ep['MI', 'KS'] == 6)
stopifnot(ep['IL', 'IA'] == 1424813)
stopifnot(ep['KY', 'IA'] == 17658)
stopifnot(ep['MO', 'IA'] == 2389932)
stopifnot(ep['OK', 'CA'] == 16762)
stopifnot(ep['MO', 'CA'] == 645)

epl <- log10(ep +1)
epl <- data.matrix(epl)
## We discard the within-state flows because they do not enter into
## the analysis. Thus it is mostly a matter of how the plots will
## look, and because the original spreadsheet had no data on
## within-state flows (they were added for some other analyses based
## on an assumption of 90% of total flow), there's no good reason to
## include them in the plot.
diag(epl) <- NA

key <- match(colnames(ep), colnames(CC))
CC <- CC[key,key]



## Get great circle distance

key <- match(colnames(CC), state.abb)

cx <- state.center$x[key]
cy <- state.center$y[key]

n <- ncol(CC)
centerDists <- matrix(nrow=n, ncol=n)

# Calculates the geodesic distance between two points specified by radian latitude/longitude using the
# Haversine formula (hf)
# source: http://www.r-bloggers.com/great-circle-distance-calculations-in-r/
gcd.hf <- function(long1, lat1, long2, lat2) {
      R <- 6371 # Earth mean radius [km]
        delta.long <- (long2 - long1)
        delta.lat <- (lat2 - lat1)
        a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
        c <- 2 * asin(min(1,sqrt(a)))
        d = R * c
        return(d) # Distance in km
  }

deg2rad <- function(deg) return(deg*pi/180)

getCenterDist <- function(state1,state2){
    long1 <- deg2rad(cx[state1])
    long2 <- deg2rad(cx[state2])
    lat1 <- deg2rad(cy[state1])
    lat2 <- deg2rad(cy[state2])
    gcd.hf(long1, lat1, long2, lat2)
}

for(i in seq_len(n)){
    for(j in seq_len(n)){
        centerDists[i,j] <- getCenterDist(i, j)
    }
}

colnames(centerDists) <- rownames(centerDists) <- colnames(CC)

## hypothesis tesing

getMat <- function(x) switch(x, 'shipment'=epl, 'cor'=CC,
                             'gcd'=-centerDists, 'sharedBord'=nhood)

doTest <- function(M1, M2, symmetrize=FALSE, ...){
    x <- getMat(M1)
    y <- getMat(M2)
    if(symmetrize){
        x <- (x + t(x))*0.5
        y <- (y + t(y))*0.5
    }
    mantel(x, y, ...)
}

mats <- c('shipment', 'cor', 'gcd', 'sharedBord')
methods <- c('spearman', 'pearson')
symmetrize <- c(TRUE, FALSE)
des <- expand.grid(M1=mats, M2=mats, method=methods, symmetrize=symmetrize, stringsAsFactors=FALSE)
des <- des[des$M1 != des$M2, ]
des$permutations <- 10000

res <- list()
for(i in seq_len(nrow(des))){
    print(i)
    res[[i]] <- do.call(doTest, as.list(des[i,]))
}

des$r <- sapply(res, '[[', 'statistic')
des$pValues <- sapply(res, '[[', 'signif')

sink('mantel-table.txt')
pander(des)
sink()

save.image(file='mantel-testing-checkpoint1.RData')

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
              legend.args=list(text=expression(paste(Log[10], '(transport flow)')),
                  line=2.9, side=4), panelLab='A')
par(mar=c(4,4.5,1,.5))
CCnoDiag <- CC
diag(CCnoDiag) <- NA
makeImagePlot(M=data.matrix(CCnoDiag), orderings=orderings, xlab='Leading state', ylab='Lagging state', 
              legend.args=list(text='Cross correlation', line=2.9, side=4), panelLab='B')
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
    g <- g + geom_hline(yintercept=0, size=0.5, col='grey')
    g <- g + geom_step(direction="vh")
    g <- g + facet_wrap(~variable, ncol=2, scales='free_y')
    g <- g + scale_x_date()
    g <- g + scale_y_discrete(breaks=pretty_breaks(n=2))
    g <- g + xlab('2013-2014 Date') + ylab('Cases')
    g <- g + theme_classic()
    g <- g + theme(strip.background = element_blank(),
                   strip.text.x = element_blank())
    g <- g + theme(plot.margin=unit(c(0,2,2,0),"mm"))
    g <- g + geom_text(data=labdf, hjust=0, vjust=1,
                       aes(x=minx, y=maxy, label=variable))
    g <- g + theme(axis.title.x = element_text(vjust=-0.5))
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
    grid.text(label="A", x=unit(0, "npc") - unit(1.5, "lines"),
              y=unit(1, "npc"), just= "left")
    popViewport()
    pushViewport(viewport(layout.pos.col=2, layout.pos.row=1))
    print(g, newpage=FALSE)
    grid.text(label="B", x=unit(0, "npc"), y=unit(1, "npc"), just= "left")
    popViewport()
}

cairo_ps(filename = 'plotMatrixWithTimeSeries.eps', width = 19/2.54, height = 6.4/2.54)
tmpf()
dev.off()

cairo_pdf(filename = 'plotMatrixWithTimeSeries.pdf', width = 19/2.54, height = 6.4/2.54)
tmpf()
dev.off()



