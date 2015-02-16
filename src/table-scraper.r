
setwd('/home/eamon/src/2015pedv/data')
system('pdftops -f 5 -l 5 SECD_Situation_Report_150212.pdf p5.ps')
PostScriptTrace('p5.ps', 'p5.xml', charpath=FALSE)
foo <- readPicture('p5.xml')
grid.newpage()

isStroke <- sapply(foo@paths, inherits, "PictureStroke")
strks <- foo[isStroke]

isVertical <- function(x) isTRUE(all.equal(as.numeric(diff(x@x)), 0))
test <- sapply(strks@paths, isVertical)
tmpf <- function(x) x@x[1]
xbreaks <- unique(sapply(strks@paths[test], tmpf))
tmpf <- function(x) x@y[1]
ybreaks <- unique(sapply(strks@paths[!test], tmpf))

isText <- sapply(foo@paths, inherits, "PictureText")
txt <- foo[isText]
txtp <- txt@paths
tmpf <- function(x) grepl("^\\s*$", x@string)
test <- sapply(txtp, tmpf)
txtp <- txtp[!test]
tmpf <- function(x) list(x=x@x, y=x@y, string=x@string)
txt <- lapply(txtp, tmpf)

df <- data.frame(x=sapply(txtp, function(x) x@x),
                 y=sapply(txtp, function(x) x@y),
                 string=sapply(txtp, function(x) x@string),
                 stringsAsFactors=FALSE)
colBreaks <- c(0, sort(xbreaks), Inf)
colLabs <- c('week', 'total', 'AZ', 'CA', 'CO', 'HI', 'IA', 'IL', 'IN', 'KS',
             'KY', 'MI', 'MN', 'MO', 'MT', 'NC', 'NE', 'NV', 'NY', 'OH', 'OK',
             'PA', 'SD', 'TN', 'TX', 'UT', 'VA', 'WI', 'WY')
df$cols <- cut(df$x, breaks=sort(colBreaks), labels=colLabs)

rowBreaks <- c(sort(ybreaks), max(ybreaks) + 140*1:2, Inf)
rowLabs <- c('Title', 'headers', paste0('week', seq_len(length(ybreaks))))
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
test <- mlt[, 'col'] != 'week' & !(mlt[, 'row'] %in% c('headers', 'Title'))
mltt <- mlt[test,]
vals <- strsplit(mltt[, 'apld'], '/')
vals <- t(sapply(vals, as.numeric))
colnames(vals) <- c('confirmed', 'presumptive')
D <- data.frame(mltt, vals[rownames(mltt), ])
D$confPresum <- D$confirmed + D$presumptive

