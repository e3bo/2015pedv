library(plyr)
library(reshape2)

setwd('/home/eamon/src/2015pedv/data')
pages <- c(5:6,10:11)
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
    tmpf <- function(x) x@y[1]
    ybreaks <- unique(sapply(strks@paths[!test], tmpf))
    list(xbreaks=xbreaks, ybreaks=ybreaks)
}
breaks <- lapply(strks, tmpf)

tmpf <- function(x, extra){
    n <- length(x$xbreaks)
    n <- n - extra
    xbreaks <- x$xbreaks[seq_len(n)]
    x$colBreaks <- c(0, sort(xbreaks), Inf)
    x
}
breaks <- mapply(tmpf, breaks, extra=c(0,0,1,1), SIMPLIFY=FALSE)

tmpf <- function(x, upper, lower) {
    yb <- sort(x$ybreaks)
    n <- length(yb)
    U <- max(yb) 
    yb <- c(yb[1:(n-1)], U + 140*upper)
    n <- length(yb)
    L <- min(yb)
    yb <- c(L - 140*lower, yb[2:n])
    x$rowBreaks <- c(sort(yb), Inf)
    x    
}
breaks <- mapply(tmpf, breaks, upper=list(0:2,0:2,c(0,1,4),c(0,3)),
                 lower=list(0,c(0,1,3),0,0:2), SIMPLIFY=FALSE)

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

colLabsT3 <- c('week', 'total', 'AZ', 'CA', 'CO', 'HI', 'IA', 'IL', 'IN', 'KS',
               'KY', 'MI', 'MN', 'MO', 'MT', 'NC', 'NE', 'NV', 'NY', 'OH', 'OK',
               'PA', 'SD', 'TN', 'TX', 'UT', 'VA', 'WI', 'WY')
colLabsT4 <- c('week', 'total', 'AZ', 'CA', 'CO', 'HI', 'IA', 'ID', 'IL', 'IN', 'KS',
               'KY', 'MI', 'MN', 'MO', 'MT', 'NC', 'ND', 'NE', 'NV', 'NY', 'OH', 'OK',
               'PA', 'SD', 'TN', 'TX', 'UT', 'VA', 'WI', 'WY', 'UNK')
colLabs <- list(colLabsT3, colLabsT3, colLabsT4, colLabsT4)

nwks <- sapply(breaks, function(x) length(x$rowBreaks) - 3)
wknos <- lapply(nwks, function(x) 1:x)
wknos[[2]] <- wknos[[2]] + nwks[1]
wknos[[4]] <- wknos[[4]] + nwks[3]

tmpf <- function(x) c('Title', 'headers', paste0('week', x))
rowLabs <- sapply(wknos, tmpf)
rowLabs[[2]] <- sub('week37', 'total', rowLabs[[2]])
rowLabs[[4]] <- sub('week37', 'total', rowLabs[[4]])

tmpf <- function(D, brks, rl, cl) {
    D$cols <- cut(D$x, breaks=brks$colBreaks, labels=cl)
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
mlt <- mapply(tmpf, D=txtdf, brks=breaks, rl=rowLabs, cl=colLabs)

m3 <- mlt[1:2]
m8 <- mlt[3:4]

getVals <- function(x) {
    M <- do.call(rbind, x)
    test <- M[, 'col'] != 'week' & !(M[, 'row'] %in% c('headers', 'Title', 'total'))
    M[test,]
}

tmpf <- function(m){
    V <- getVals(m)
    vals <- strsplit(V[, 'apld'], '/')
    vals <- t(sapply(vals, as.numeric))
    colnames(vals) <- c('confirmed', 'presumptive')
    tab3 <- data.frame(V, vals[rownames(mltt), ])
    tab3$confPresum <- tab3$confirmed + tab3$presumptive
    tab3
}
tab3 <- tmpf(m3)

tmpf <- function(m){
    V <- getVals(m)
    vals <- as.numeric(V[, 'apld'])
    tab8 <- data.frame(V, accessions=vals)
    tab8
}
tab8 <- tmpf(m8)

getDates <- function(x) {
    M <- do.call(rbind, x)
    test <- M[, 'col'] == 'week' & !(M[, 'row'] %in% c('headers', 'Title', 'total'))
    M <- M[test,]
    test <- M[, 'apld'] == '1/18/15'
    if(any(test)){
        M[test, 'apld'] <- '1/18/2015'
    }
    dt <- as.Date(M[, 'apld'], format="%m/%d/%Y")
    data.frame(row=M[, 'row'], date=dt)
}
dt3 <- getDates(m3)
dt8 <- getDates(m8)

tab3 <- merge(tab3, dt3)
tab8 <- merge(tab8, dt8)

tab3w <- dcast(tab3, date~col, value.var="confPresum")
tab3w[is.na(tab3w)] <- 0
tab8w <- dcast(tab8, date~col, value.var="accessions")
tab8w[is.na(tab8w)] <- 0

m8 <- melt(tab8w, id.var='date', value.name='accessions')
m3 <- melt(tab3w, id.var='date', value.name='confPresum')
mg <- merge(m3, m8)
test <- mg$variable=='total'
mg <- mg[!test, ]

print('Correlation test of weekly counts of accessions vs. confirmed + presumptive premises:')
cor.test(mg$confPresum, mg$accessions, method='spearman')

tmpf <- function(x){
    colSums(x[, c('confPresum', 'accessions')])
}
mgCum <- ddply(mg, c('variable'), tmpf)
print('Correlation test of cumulative values:')
cor.test(mgCum$confPresum, mgCum$accessions, method='spearman')

save.image('table-scraper.r')
