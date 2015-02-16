library(reshape2)

setwd('/home/eamon/src/2015pedv/data')
pages <- 5:6
tmpf <- function(x) paste0('pdftops -f ', x, ' -l ', x,
                           ' SECD_Situation_Report_150212.pdf p', x, '.ps')
cmds <- sapply(pages, tmpf)
invisible(lapply(cmds, system))
ps <- sapply(strsplit(cmds, ' '), '[[', 7)
xml <- gsub('\\.ps', '\\.xml', ps)
invisible(mapply(PostScriptTrace, file=ps, outfilename=xml, charpath=FALSE))
pics <- lapply(xml, readPicture)

tmpf <- function(x, what) {
    test <- sapply(x@paths, inherits, what)
    x[test]
}
strks <- lapply(pics, tmpf, what="PictureStroke")
txts <- lapply(pics, tmpf, what="PictureText")

tmpf <- function(strks) {
    isVertical <- function(x) isTRUE(all.equal(as.numeric(diff(x@x)), 0))
    test <- sapply(strks@paths, isVertical)
    tmpf <- function(x) x@x[1]
    xbreaks <- unique(sapply(strks@paths[test], tmpf))
    colBreaks <- c(0, sort(xbreaks), Inf)
    tmpf <- function(x) x@y[1]
    ybreaks <- unique(sapply(strks@paths[!test], tmpf))
    rowBreaks <- c(sort(ybreaks), max(ybreaks) + 140*1:2, Inf)
    list(colBreaks=colBreaks, rowBreaks=rowBreaks)
}
breaks <- lapply(strks, tmpf)

rb2 <- breaks[[2]]$rowBreaks
breaks[[2]]$rowBreaks <- c(min(rb2) - 140*c(3,1), rb2)

tmpf <- function(txt){
    txtp <- txt@paths
    tmpf <- function(x) grepl("^\\s*$", x@string)
    test <- sapply(txtp, tmpf)
    txtp <- txtp[!test]
    data.frame(x=sapply(txtp, function(x) x@x),
               y=sapply(txtp, function(x) x@y),
               string=sapply(txtp, function(x) x@string),
               stringsAsFactors=FALSE)
}
txtdf <- lapply(txts, tmpf)

colLabs <- c('week', 'total', 'AZ', 'CA', 'CO', 'HI', 'IA', 'IL', 'IN', 'KS',
             'KY', 'MI', 'MN', 'MO', 'MT', 'NC', 'NE', 'NV', 'NY', 'OH', 'OK',
             'PA', 'SD', 'TN', 'TX', 'UT', 'VA', 'WI', 'WY')

nwks <- sapply(breaks, function(x) length(x$rowBreaks) - 3)
wknos <- lapply(nwks, function(x) 1:x)
wknos[[2]] <- wknos[[2]] + nwks[1]

tmpf <- function(x) c('Title', 'headers', paste0('week', x))
rowLabs <- sapply(wknos, tmpf)
rowLabs[[2]] <- sub('week37', 'total', rowLabs[[2]])

tmpf <- function(D, brks, rl) {
    D$cols <- cut(D$x, breaks=brks$colBreaks, labels=colLabs)
    D$rows <- cut(D$y, breaks=brks$rowBreaks, labels=rev(rl))
    splt <- split(D, list(D$col, D$row), drop=TRUE)
    tmpf <- function(piece) {
        ord <- order(-piece$y, piece$x)
        string <- piece$string[ord]
        string <- paste(string, collapse='')
        string
    }
    apld <- sapply(splt, tmpf)
    nms <- strsplit(names(splt), '\\.')
    col <- sapply(nms, '[[', 1)
    row <- sapply(nms, '[[', 2)
    mlt <- cbind(col, row, apld)
    mlt
}
mlt <- mapply(tmpf, D=txtdf, brks=breaks, rl=rowLabs)
mlt <- do.call(rbind, mlt)

test <- mlt[, 'col'] != 'week' & !(mlt[, 'row'] %in% c('headers', 'Title', 'total'))
mltt <- mlt[test,]
vals <- strsplit(mltt[, 'apld'], '/')
vals <- t(sapply(vals, as.numeric))
colnames(vals) <- c('confirmed', 'presumptive')
tab3 <- data.frame(mltt, vals[rownames(mltt), ])
tab3$confPresum <- tab3$confirmed + tab3$presumptive

test <- mlt[, 'col'] == 'week' & !(mlt[, 'row'] %in% c('headers', 'Title', 'total'))
mltt <- mlt[test,]
mltt[mltt[, 'apld']=='1/18/15', 'apld'] <- '1/18/2015'
dt <- as.Date(mltt[, 'apld'], format="%m/%d/%Y")
dt <- data.frame(row=mltt[, 'row'], date=dt)
tab3 <- merge(tab3, dt)



dcast(tab3, date~col, value.var="confPresum")


dcast(tab3, row~col, value.var="confPresum")

colBreaks <- c(0, sort(xbreaks), Inf)
df$cols <- cut(df$x, breaks=sort(colBreaks), labels=colLabs)


df$rows <- cut(df$y, breaks=rowBreaks, labels=rev(rowLabs))

splt <- split(df, list(df$col, df$row), drop=TRUE)

tmpf <- function(piece) {
    ord <- order(-piece$y, piece$x)
    string <- piece$string[ord]
    string <- paste(string, collapse='')
    string
}
apld <- sapply(splt, tmpf)
nms <- strsplit(names(splt), '\\.')
col <- sapply(nms, '[[', 1)
row <- sapply(nms, '[[', 2)
mlt <- cbind(col, row, apld)



