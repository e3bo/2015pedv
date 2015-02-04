#!/usr/bin/Rscript

stateShipments <- read.csv('../data/StateShipments.csv', skip=2, na.strings="")
desel <- c("Origin", "Total", "X")
test <- !(colnames(stateShipments) %in% desel)
ss <- stateShipments[, test]

desel <- c("Total", "NA.")
rn <- make.names(stateShipments$Origin)
test <- !(rn %in% desel)
ss <- ss[test, ]
rownames(ss) <- rn[test]
stopifnot(rownames(ss) == colnames(ss))
ss[is.na(ss)] <- 0

key <- match(colnames(ss), c(make.names(state.name), 'District.of.Columbia',
                             'Mexico', 'Canada', 'Other.states'))
nms <- c(state.abb, 'DC', 'Mexico', 'Canada', 'otherStates')[key]

rownames(ss) <- colnames(ss) <- nms
write.csv(data.matrix(ss),
          'shipment-flows-origins-on-rows-dests-on-columns.csv',
          row.names=TRUE)
