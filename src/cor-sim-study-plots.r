#!/usr/bin/Rscript

library(ggplot2)
library(grid)
library(gridBase)
library(igraph)

stat_sum_df <- function(fun, geom="crossbar", ...) {
    stat_summary(fun.data=fun, colour="red", geom=geom, width=0.2, ...)
}
load('corSims-layer1.RData')

test <- rdf$expectedIntros == 1 & rdf$couplingFactor == 1
df <- rdf[test,]
test <- df$comSize == 64

colScale <- c("degreeBalanced"="#4477AA", "recipUnbalanced"="#DDCC77", "unbalanced"="#CC6677")
g <- ggplot(data=df[test,], aes(x=cohesion, y=cc1, group=scheme, color=scheme))
g <- g + geom_point(position=position_jitter(width=0.05), alpha=0.75)
g <- g + scale_x_log10()
g <- g + labs(x='Vertex\nconnectivity', y='Cross correlation')
g <- g + scale_color_manual(values=colScale)
g <- g + stat_sum_df("mean_cl_normal", geom = "smooth")
g <- g + scale_y_continuous(limits=c(-0.15,0.85))
g <- g + theme_classic()
## Because label from adjacent plot gets partially covered with opaque background
g <- g + theme(plot.background = element_rect(fill = "transparent",colour = NA))
pCohesion <- g + theme(legend.position='none')

g <- ggplot(data=df[test,], aes(x=nedges, y=cc1, group=scheme, color=scheme))
g <- g + geom_point(position=position_jitter(width=0.05), alpha=0.75)
g <- g + scale_x_log10()
g <- g + labs(x='Edges in\ncut set', y='Cross correlation')
g <- g + stat_sum_df("mean_cl_normal", geom = "smooth")
g <- g + scale_color_manual(values=colScale)
g <- g + scale_y_continuous(limits=c(-0.15,0.85))
g <- g + theme_classic()
#g <- g + theme(axis.title.x=element_text(hjust=1))
pNedges <- g + theme(legend.position='none')

tt <- df$nedges == df$comSize^2
dff <- df[tt,]

g <- ggplot(data=dff, aes(x=R0baseline, y=cc1))
g <- g + geom_point(position=position_jitter(width=0.05), alpha=0.75)
g <- g + labs(x='R0 baseline\n', y='Cross correlation')
g <- g + scale_x_log10(breaks=c(.5,1,2))
g <- g + scale_y_continuous(limits=c(-0.15,0.85))
g <- g + theme_classic()
pR0 <- g

test <- rdf$expectedIntros == 1 & rdf$R0baseline == 1

df <- rdf[test,]
tt <- df$nedges == df$comSize^2
dff <- df[tt,]

g <- ggplot(data=dff, aes(x=couplingFactor, y=cc1))
g <- g + geom_point(position=position_jitter(width=0.05), alpha=0.75)
g <- g + labs(x='Capacity\nfactor', y='Cross correlation')
g <- g + scale_x_log10(breaks=c(.5,1,2))
g <- g + scale_y_continuous(limits=c(-0.15,0.85))
g <- g + theme_classic()
pCoupFac <- g

# make wiring diagrams

gr <- gu <- gbal <- graph.empty(8, directed=FALSE)
lay <- layout.grid(gr, width=2)
gbal <- gbal + edges(1,2, 3,4, 5,6, 7,8)
gbal <- gbal + edges(3,2, 5,4, 7,6)

gu <- gu + edges(7,8, 7,6, 7,4, 7,2)
gu <- gu + edges(1,2, 1,4, 1,6)

gr <- gr + edges(7,8, 7,6, 7,4, 7,2)
gr <- gr + edges(2,1, 2,3, 2,5)

nodeLayout <- layout.grid(gr, width=2)
plotNetwork <- function(x) plot(x, layout=nodeLayout, vertex.label='', vertex.color='black',
                                vertex.size=20, edge.color='black')

## create plot

tmpf <- function(){
    plot.new() # seems necessary for cairo_pdf and cairo_ps devices
    grid.newpage()
    pushViewport(viewport(width=unit(17.8, 'cm'), height=unit(6, 'cm')))
    lay <- grid.layout(1,5, widths=unit(c(1.1,1,1,1,1), 'null'),
                       heights=unit(c(1), 'null'))
    pushViewport(viewport(layout=lay))
    pushViewport(viewport(layout.pos.col=1, layout.pos.row=1))
    print(pR0, newpage=FALSE)
    grid.text(label="A", x=unit(0, "npc") - unit(0.25, "lines"),
              y=unit(1, "npc"), just= "left")
    popViewport()
    pushViewport(viewport(layout.pos.col=2, layout.pos.row=1))
    print(pCoupFac + labs(y=NULL), newpage=FALSE)
    grid.text(label="B", x=unit(0, "npc") - unit(0., "lines"),
              y=unit(1, "npc"), just= "left")
    popViewport()
    pushViewport(viewport(layout.pos.col=3, layout.pos.row=1))
    print(pNedges + labs(y=NULL), newpage=FALSE)
    grid.text(label="C", x=unit(0, "npc") - unit(0., "lines"),
              y=unit(1, "npc"), just= "left")
    popViewport()
    pushViewport(viewport(layout.pos.col=4, layout.pos.row=1))
    print(pCohesion + labs(y=NULL), newpage=FALSE)
    grid.text(label="D", x=unit(0, "npc") - unit(0., "lines"),
              y=unit(1, "npc"), just= "left")
    popViewport()
    pushViewport(viewport(layout.pos.col=5, layout.pos.row=1))
    lay <- grid.layout(5,4, widths=unit(c(.25,.2, 2.2,.5), 'null'),
                       heights=unit(c(.75, .75, .75, .75,1), 'null'))
    pushViewport(viewport(layout=lay))
    pushViewport(viewport(layout.pos.col=4, layout.pos.row=2))
    par(plt=gridPLT(), new=TRUE)
    plotNetwork(gbal)
    popViewport()
    pushViewport(viewport(layout.pos.col=4, layout.pos.row=3))
    par(plt=gridPLT(), new=TRUE)
    plotNetwork(gu)
    popViewport()
    pushViewport(viewport(layout.pos.col=4, layout.pos.row=4))
    par(plt=gridPLT(), new=TRUE)
    plotNetwork(gr)
    popViewport()
    pushViewport(viewport(layout.pos.col=3, layout.pos.row=2))
    grid.text(label='Balanced',  gp=gpar(fontsize=10))
    popViewport()
    pushViewport(viewport(layout.pos.col=3, layout.pos.row=3))
    grid.text(label='Unbalanced',  gp=gpar(fontsize=10))
    popViewport()
    pushViewport(viewport(layout.pos.col=3, layout.pos.row=4))
    grid.text(label='Reciprocally\nunbalanced', gp=gpar(fontsize=10))
    popViewport()
    pushViewport(viewport(layout.pos.col=3, layout.pos.row=1))
    grid.text(label='Wiring scheme', gp=gpar(fontsize=12, fontface=1))
    popViewport()
    pushViewport(viewport(layout.pos.col=2, layout.pos.row=2))
    grid.circle(gp=gpar(fill=colScale["degreeBalanced"], alpha=0.75, col="transparent"))
    popViewport()
    pushViewport(viewport(layout.pos.col=2, layout.pos.row=2))
    grid.circle(gp=gpar(fill=colScale["degreeBalanced"], alpha=0.75, col="transparent"))
    popViewport()
    pushViewport(viewport(layout.pos.col=2, layout.pos.row=3))
    grid.circle(gp=gpar(fill=colScale["unbalanced"], alpha=0.75, col="transparent"))
    popViewport()
    pushViewport(viewport(layout.pos.col=2, layout.pos.row=4))
    grid.circle(gp=gpar(fill=colScale["recipUnbalanced"], alpha=0.75, col="transparent"))
    popViewport()
}

cairo_ps(filename='sim-panel.eps', onefile=TRUE, width=7.1, height=2.45)
tmpf()
dev.off()

cairo_pdf(filename='tmp.pdf', onefile=TRUE, width=7.1, height=2.45)
tmpf()
dev.off()
# First page is blank (due to plot.new?), so:
system('pdfseparate -f 2 -l 2 tmp.pdf sim-panel.pdf && rm tmp.pdf')
