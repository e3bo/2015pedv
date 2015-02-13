#!/usr/bin/Rscript

library(ggplot2)
library(grid)

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

pdf('sim-panel.pdf', width=17.8 / 2.54, height=6/2.54)
lay <- grid.layout(1,5, widths=unit(c(1.1,1,1,1,1), 'null'),
                   heights=unit(c(1), 'null'))
grid.newpage()
pushViewport(viewport(layout=lay))
pushViewport(viewport(layout.pos.col=1, layout.pos.row=1))
print(pR0, newpage=FALSE)
popViewport()
pushViewport(viewport(layout.pos.col=2, layout.pos.row=1))
print(pCoupFac + labs(y=NULL), newpage=FALSE)
popViewport()
pushViewport(viewport(layout.pos.col=3, layout.pos.row=1))
print(pNedges + labs(y=NULL), newpage=FALSE)
popViewport()
pushViewport(viewport(layout.pos.col=4, layout.pos.row=1))
print(pCohesion + labs(y=NULL), newpage=FALSE)
popViewport()
dev.off()
