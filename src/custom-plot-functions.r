## Custom function to suppress X and Y tick labels at corner plots as
## well as gridlines. Also, translate the variable names to those used
## for plotting.
ebo_ggally_diagAxis <-
    function (data, mapping, labelSize = 3, labelXPercent = 0.55,
              labelYPercent = 0.55, labelHJust = 0.5, labelVJust = 0.5,
              gridLabelSize = 3, suppressY=FALSE, suppressX=FALSE, ...)
{
    mapping$y <- NULL
    numer <- !((is.factor(data[, as.character(mapping$x)])) ||
               (is.character(data[, as.character(mapping$x)])))
    if (numer) {
        label <- switch(mapping$x, 'cor'='Cross\ncorrelation',
                        'shipment'='Log10(#swine / y)',
                        'gcd'='-Distance\n(1,000 km)')
        if(mapping$x=='gcd'){
            # to avoid crowding in small figures
            data <- data/1000
        }
        xmin <- min(data[, as.character(mapping$x)], na.rm = TRUE)
        xmax <- max(data[, as.character(mapping$x)], na.rm = TRUE)
        xrange <- c(xmin - 0.01 * (xmax - xmin), xmax + 0.01 *
                    (xmax - xmin))
        p <- ggally_text(label = label, mapping = aes(col = "grey50"),
                         xrange = xrange, yrange = xrange, size = labelSize,
                         xP = labelXPercent, yP = labelYPercent, hjust = labelHJust,
                         vjust = labelVJust)
        p <- p + theme(panel.grid.major=element_blank())
        axisBreaks <- GGally:::get_x_axis_labels(p, xrange)
        if(suppressY){
            test <- axisBreaks$yPos == min(axisBreaks$yPos)
            axisBreaks <- axisBreaks[test, ]
        }
        if(suppressX){
            test <- axisBreaks$xPos == min(axisBreaks$xPos)
            axisBreaks <- axisBreaks[test, ]
        }
        if(!all(suppressX, suppressY)){
            pLabs <- p + geom_text(data = axisBreaks,
                                   mapping = aes_string(x = "xPos", y = "yPos", label = "lab", hjust = "hjust", vjust = "vjust"),
                                   col = "grey50", size = gridLabelSize)
        } else {
            pLabs <- p
        }
    }
    else {
        stop('not implemented')
    }
    pLabs$subType = "internal"
    pLabs$type = "label"
    pLabs
}

myPrintGGpairs <-
    function (x, leftWidthProportion = 0.2, bottomHeightProportion = 0.1,
              spacingProportion = 0.03, showStrips = NULL, newpage=TRUE, ...)
{
    ## adapted from GGally:::print.ggpairs
    vplayout <- function (x, y)
    {
        viewport(layout.pos.row = x, layout.pos.col = y)
    }

    is_blank_plot <- function (p)
    {
        if (!is.null(p$subType) && !is.null(p$type))
            p$subType == "blank" && p$type == "blank"
        else FALSE
    }
    plotObj <- x
    if (identical(plotObj$axisLabels, "internal")) {
        v1 <- viewport(y = unit(0.5, "npc") + unit(1, "lines"),
                       x = unit(0.5, "npc") - unit(1, "lines"),
                       width = unit(1, "npc") - unit(1, "lines"), height = unit(1,
                                                                      "npc") - unit(2, "lines"))
    }
    else {
        v1 <- viewport(width = unit(1, "npc") - unit(3, "lines"),
                       height = unit(1, "npc") - unit(3, "lines"))
    }
    numCol <- length(plotObj$columns)
    if (identical(plotObj$axisLabels, "show")) {
        showLabels <- TRUE
        viewPortWidths <- c(leftWidthProportion, 1, rep(c(spacingProportion,
                                                          1), numCol - 1))
        viewPortHeights <- c(rep(c(1, spacingProportion), numCol -
                                 1), 1, bottomHeightProportion)
    }
    else {
        showLabels <- FALSE
        viewPortWidths <- c(1, rep(c(spacingProportion, 1), numCol -
                                   1))
        viewPortHeights <- c(rep(c(1, spacingProportion), numCol -
                                 1), 1)
    }
    viewPortCount <- length(viewPortWidths)
    v2 <- viewport(layout = grid.layout(viewPortCount, viewPortCount,
                       widths = viewPortWidths, heights = viewPortHeights))
    if(newpage) grid.newpage()
    if (plotObj$title != "") {
        pushViewport(viewport(height = unit(1, "npc") - unit(0.4,
                                  "lines")))
        grid.text(plotObj$title, x = 0.5, y = 1, just = c(0.5,
                                                     1), gp = gpar(fontsize = 15))
        popViewport()
    }
    if (!identical(plotObj$axisLabels, "internal")) {
        pushViewport(viewport(width = unit(1, "npc") - unit(2,
                                  "lines"), height = unit(1, "npc") - unit(3, "lines")))
        pushViewport(viewport(layout = grid.layout(viewPortCount,
                                  viewPortCount, widths = viewPortWidths, heights = viewPortHeights)))
        for (i in 1:numCol) {
            grid.text(plotObj$columnLabels[i], 0, 0.5, rot = 90,
                      just = c("centre", "centre"), vp = vplayout(as.numeric(i) *
                                                        2 - 1, 1))
        }
        popViewport()
        popViewport()
        pushViewport(viewport(width = unit(1, "npc") - unit(3,
                                  "lines"), height = unit(1, "npc") - unit(2, "lines")))
        pushViewport(viewport(layout = grid.layout(viewPortCount,
                                  viewPortCount, widths = viewPortWidths, heights = viewPortHeights)))
        for (i in 1:numCol) {
            grid.text(plotObj$columnLabels[i], 0.5, 0, just = c("centre",
                                                           "centre"), vp = vplayout(ifelse(showLabels, 2 *
                                                                          numCol, 2 * numCol - 1), ifelse(showLabels, 2 *
                                                                                                          i, 2 * i - 1)))
        }
        popViewport()
        popViewport()
    }
    pushViewport(v1)
    pushViewport(v2)
    for (rowPos in 1:numCol) {
        for (columnPos in 1:numCol) {
            p <- getPlot(plotObj, rowPos, columnPos)
            if (is_blank_plot(p)) {
                next
            }
            pGtable <- ggplot_gtable(ggplot_build(p))
            if (columnPos == 1 && showLabels) {
                if (identical(plotObj$verbose, TRUE)) {
                    print("trying left axis")
                }
                pAxisLabels <- gtable_filter(pGtable, "axis-l")
                grobLength <- length(pAxisLabels$grobs)
                leftAxisLayoutHeight <- rep(c(0.1, 1), grobLength)[-1]
                leftAxisLayoutHeightUnits <- rep(c("lines", "null"),
                                                 grobLength)[-1]
                vpLAxis <- viewport(layout = grid.layout(nrow = 2 *
                                        grobLength - 1, ncol = 1, widths = unit(1,
                                                                      "null"), heights = unit(leftAxisLayoutHeight,
                                                                                   leftAxisLayoutHeightUnits)))
                pushViewport(vplayout(rowPos * 2 - 1, 1))
                pushViewport(vpLAxis)
                for (lAxisPos in 1:grobLength) {
                    pushViewport(vplayout(lAxisPos * 2 - 1, 1))
                    grid.draw(pAxisLabels$grobs[[lAxisPos]])
                    popViewport()
                }
                popViewport()
                popViewport()
            }
            if (rowPos == numCol && showLabels) {
                if (identical(plotObj$verbose, TRUE)) {
                    print("trying bottom axis")
                }
                pAxisLabels <- gtable_filter(pGtable, "axis-b")
                grobLength <- length(pAxisLabels$grobs)
                botAxisLayoutWidth <- rep(c(0.1, 1), grobLength)[-1]
                botAxisLayoutWidthUnits <- rep(c("lines", "null"),
                                               grobLength)[-1]
                vpBAxis <- viewport(layout = grid.layout(nrow = 1,
                                        ncol = 2 * grobLength - 1, heights = unit(1,
                                                                       "null"), widths = unit(botAxisLayoutWidth,
                                                                                    botAxisLayoutWidthUnits)))
                pushViewport(vplayout(2 * numCol, 2 * columnPos))
                pushViewport(vpBAxis)
                for (bAxisPos in 1:grobLength) {
                    pushViewport(vplayout(1, bAxisPos * 2 - 1))
                    grid.draw(pAxisLabels$grobs[[bAxisPos]])
                    popViewport()
                }
                popViewport()
                popViewport()
            }
            layoutNames <- c("panel")
            allLayoutNames <- c("panel", "strip-right", "strip-top")
            if (is.null(showStrips)) {
                pShowStrips <- (!is.null(p$type)) && (!is.null(p$subType))
                if (pShowStrips) {
                    if (columnPos == numCol) {
                        layoutNames <- c(layoutNames, "strip-right")
                    }
                    if (rowPos == 1) {
                        layoutNames <- c(layoutNames, "strip-top")
                    }
                }
            }
            else if (showStrips) {
                layoutNames <- allLayoutNames
            }
            if (!is.null(p$axisLabels)) {
                if (p$axisLabels %in% c("internal", "none")) {
                    layoutNames <- allLayoutNames
                }
            }
            layoutRows <- pGtable$layout$name %in% layoutNames
            layoutInfo <- pGtable$layout[layoutRows, ]
            layoutTB <- layoutInfo[, c("t", "b")]
            layoutLR <- layoutInfo[, c("l", "r")]
            pPanel <- pGtable[min(layoutTB):max(layoutTB), min(layoutLR):max(layoutLR)]
            pushViewport(vplayout(2 * rowPos - 1, ifelse(showLabels,
                                                         2 * columnPos, 2 * columnPos - 1)))
            suppressMessages(suppressWarnings(grid.draw(pPanel)))
            popViewport()
        }
    }
    popViewport()
    popViewport()
}

image.plot.ebo <- function (..., add = FALSE, nlevel = 64, horizontal = FALSE,
                            legend.shrink = 0.9, legend.width = 1.2,
                            legend.mar = ifelse(horizontal, 3.1, 5.1),
                            legend.lab = NULL, legend.line = 2,
                            graphics.reset = FALSE, bigplot = NULL,
                            smallplot = NULL, legend.only = FALSE,
                            col = tim.colors(nlevel), lab.breaks = NULL,
                            axis.args = NULL, legend.args = NULL,
                            midpoint = FALSE, border = NA, lwd = 1,
                            panelLab=NULL) {
    old.par <- par(no.readonly = TRUE)
    info <- imageplot.info(...)
    if (add) {
        big.plot <- old.par$plt
    }
    if (legend.only) {
        graphics.reset <- TRUE
    }
    if (is.null(legend.mar)) {
        legend.mar <- ifelse(horizontal, 3.1, 5.1)
    }
    temp <- imageplot.setup(add = add, legend.shrink = legend.shrink,
                            legend.width = legend.width, legend.mar = legend.mar,
                            horizontal = horizontal, bigplot = bigplot, smallplot = smallplot)
    smallplot <- temp$smallplot
    bigplot <- temp$bigplot
    if (!legend.only) {
        if (!add) {
            par(plt = bigplot)
        }
        if (!info$poly.grid) {
            image(..., add = add, col = col, axes=FALSE)
            addAxis1 <- function(x, y, z,...){
                labs <- paste(names(y), c('', '      '))
                axis(1, at=y, labels = labs, las = 2, line = -0.5, tick = 0,
        cex.axis = 0.8)
            }
            addAxis1(...)
            addAxis2 <- function(x, y, z,...){
                labs <- paste(names(x), c('', '      '))
                axis(2, at=x, labels = labs, las = 1, line = -0.5, tick = 0,
        cex.axis = 0.8)
            }
            addAxis2(...)
            mtext(panelLab, side=2, las=2, at=par()$usr[2], line=par()$mar[2] - 1, cex=1)
        }
        else {
            poly.image(..., add = add, col = col, midpoint = midpoint,
                       border = border, lwd.poly = lwd)
        }
        big.par <- par(no.readonly = TRUE)
    }
    if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
        par(old.par)
        stop("plot region too small to add legend\n")
    }
    ix <- 1
    minz <- info$zlim[1]
    maxz <- info$zlim[2]
    binwidth <- (maxz - minz)/nlevel
    midpoints <- seq(minz + binwidth/2, maxz - binwidth/2, by = binwidth)
    iy <- midpoints
    iz <- matrix(iy, nrow = 1, ncol = length(iy))
    breaks <- list(...)$breaks
    par(new = TRUE, pty = "m", plt = smallplot, err = -1)
    if (!is.null(breaks) & !is.null(lab.breaks)) {
        axis.args <- c(list(side = ifelse(horizontal, 1, 4),
                            mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2),
                            at = breaks, labels = lab.breaks), axis.args)
    }
    else {
        axis.args <- c(list(side = ifelse(horizontal, 1, 4),
                            mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2)),
                       axis.args)
    }
    if (!horizontal) {
        if (is.null(breaks)) {
            image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "",
                  ylab = "", col = col)
        }
        else {
            image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "",
                  ylab = "", col = col, breaks = breaks)
        }
    }
    else {
        if (is.null(breaks)) {
            image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "",
                  ylab = "", col = col)
        }
        else {
            image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "",
                  ylab = "", col = col, breaks = breaks)
        }
    }
    do.call("axis", axis.args)
    box()
    if (!is.null(legend.lab)) {
        legend.args <- list(text = legend.lab, side = ifelse(horizontal,
                                                   1, 4), line = legend.line)
    }
    if (!is.null(legend.args)) {
        do.call(mtext, legend.args)
    }
    mfg.save <- par()$mfg
    if (graphics.reset | add) {
        par(old.par)
        par(mfg = mfg.save, new = FALSE)
        invisible()
    }
    else {
        par(big.par)
        par(plt = big.par$plt, xpd = FALSE)
        par(mfg = mfg.save, new = FALSE)
        invisible()
    }
}

make_choropleth <- function(fill_var=c('Cases', 'Inventory'), trans='log10'){
    fill <- match.arg(fill_var)
    g <- ggplot(map_data)
    g <- g + geom_map(aes_string(map_id='state', fill=fill), color='dark grey', size=0.2, map=state_map)
    if(fill=='Cases'){
        g <- g + scale_fill_gradient(trans=trans, high="#4E2527", low="#FA644E")
    } else {
        g <- g + scale_fill_gradient(trans=trans, high="#132B43", low="#459CDA")
    }
    # expand_limits necessary to prevent error: argument "env" is missing
    g <- g + expand_limits(x = state_map$long, y = state_map$lat)
    g <- g + coord_map("azequalarea",orientation=c(30.82,-98.57,0))
    g <- g + theme_clean()
    g <- g + theme(legend.position = "top")
    g <- g + guides(fill = guide_colorbar(barwidth = 10, barheight = 0.5))
    g <- g + geom_text(aes(x=x, y=y, label=label), size=1.9, color='dark grey')
    g <- g + geom_segment(data=get_endpoints('NH'), aes(x=x, y=y, xend=xend-.2, yend=yend),
                          size=0.25, color='dark olive green')
    g <- g + geom_segment(data=get_endpoints('MA'), aes(x=x, y=y, xend=xend-.2, yend=yend),
                          size=0.25, color='dark olive green')
    g <- g + geom_segment(data=get_endpoints('RI'),
                          aes(x=x-.1, y=y, xend=xend+.4, yend=yend), size=0.25,
                          color='dark olive green')
    g <- g + geom_segment(data=get_endpoints('DE'),
                          aes(x=x, y=y, xend=xend+.2, yend=yend-.1),
                          size=0.25, color='dark olive green')
    g <- g + geom_segment(data=get_endpoints('CT'),
                          aes(x=x, y=y, xend=xend+.2, yend=yend-.1), size=0.25,
                          color='dark olive green')
    g <- g + geom_segment(data=get_endpoints('NJ'),
                          aes(x=x, y=y, xend=xend+.2, yend=yend-.1), size=0.25,
                          color='dark olive green')
    g <- g + geom_segment(data=get_endpoints('MD'),
                          aes(x=x, y=y, xend=xend+.2, yend=yend-.1), size=0.25,
                          color='dark olive green')
    g
}

makePlotMat <- function(dfl, type=c('directed', 'undirected'),
                        transform=c('ranked', 'original'), labelBreaks=FALSE, gridLabelSize=4, useSymp=FALSE,
                        ...){
    type <- match.arg(type)
    transform <- match.arg(transform)
    df <- dfl[[transform]][[type]]
    symmetrize <- ifelse(type=='directed', FALSE, TRUE)
    method <- ifelse(transform=='ranked', 'spearman', 'pearson')
    pm <- ggpairs(df, alpha=0.4, axisLabels="internal", ...)
    pal <- brewer.pal(n=9, 'Blues')
    for (i in 1:(length(mats2) - 1)){
        for(j in (i+1):length(mats2)){
            sel <- data.frame(M1=mats[i], M2=mats[j],
                              method=method, symmetrize=symmetrize)
            mg <- merge(sel, des, all.y=FALSE)
            symp <- symnum(mg$pValues, corr = FALSE,
                           cutpoints = c(0, .001,.01,.05, .1, 1),
                           symbols = c("***","**","*","."," "))
            r <- round(mg$r, 2)
            colNumber <- round(r*9, 0)
            if (useSymp){
              label <- paste(r, symp, sep='\n')
            } else {
              label <- paste(r, signif(mg$pValues, digits=2), sep='\n')
            }
            plt <- ggplot() + geom_text(label=label, aes(x=0.5, y=0.5), colour='black',
                                        size=10 * 0.3)
            plt <- plt + xlim(c(0,1)) + ylim(c(0,1))
            plt <- plt + theme(panel.background=element_rect(fill=pal[colNumber]),
                               legend.position='none')
            plt <- plt + labs(x=NULL,y=NULL)
            pm <- putPlot(pm, plt, i, j)
        }
    }
    if(labelBreaks){
        for(i in seq_along(mats2)){
            if(i == 1){
                g <- ebo_ggally_diagAxis(pm$data, mapping=list(x=mats2[i], y=mats2[i]),
                                 suppressX=FALSE, suppressY=TRUE, gridLabelSize=gridLabelSize)
            } else if (i == length(mats2)){
                g <- ebo_ggally_diagAxis(pm$data, mapping=list(x=mats2[i], y=mats2[i]),
                                         suppressX=TRUE, suppressY=FALSE, gridLabelSize=gridLabelSize)
            } else {
                g <- ebo_ggally_diagAxis(pm$data, mapping=list(x=mats2[i], y=mats2[i]),
                                         suppressX=FALSE, suppressY=FALSE, gridLabelSize=gridLabelSize)
            }
            pm <- putPlot(pm, g, i, i)
        }
    } else {
        for(i in seq_along(mats2)){
            g <- ebo_ggally_diagAxis(pm$data, mapping=list(x=mats2[i], y=mats2[i]),
                                         suppressX=TRUE, suppressY=TRUE)
            pm <- putPlot(pm, g, i, i)
        }
    }
    pm
}
