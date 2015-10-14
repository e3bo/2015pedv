#!/usr/bin/Rscript

library(igraph)
library(maps) # for state.fips
library(sp)

Sys.setlocale("LC_TIME", "C") #Needed for identical()
Sys.setlocale("LC_COLLATE", "C")

## Get case data

dataDir <- '.'

GetCaseData <- function(){
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
  wk <- as.character(ret$week)
  ret$week <- as.Date(wk, format = "%m/%d/%Y")
  ret
}
real.case.data <- GetCaseData()

GetSharedBorderAdjacencyMatrix <- function(){
    nbEdgelist <- read.csv(file.path(dataDir, 'state_neighbors_fips.txt'),
                           header=FALSE)
    data(state.fips)
    g <- graph.data.frame(nbEdgelist, directed=FALSE)
    g <- simplify(g)
    key <- match(V(g)$name, state.fips$fips)
    abb <- state.fips$abb[key]
    V(g)$name <- as.character(abb)
    as.matrix(g[, ])
}
shared.border.adjacency <- GetSharedBorderAdjacencyMatrix()

GetShipmentFlows <- function(){
    file <- file.path(dataDir, 'shipment-flows-origins-on-rows-dests-on-columns.csv')
    flows.matrix <- read.csv(file, row.names=1)
    flows.matrix <- data.matrix(flows.matrix)

    stopifnot(rownames(flows.matrix) == colnames(flows.matrix))
    diag(flows.matrix) <- NA ## these are not part of original data

    ### assertions based on inspection of original .xls file
    stopifnot(flows.matrix['CA', 'IL'] == 9415)
    stopifnot(flows.matrix['MI', 'KS'] == 6)
    stopifnot(flows.matrix['IL', 'IA'] == 1424813)
    stopifnot(flows.matrix['KY', 'IA'] == 17658)
    stopifnot(flows.matrix['MO', 'IA'] == 2389932)
    stopifnot(flows.matrix['OK', 'CA'] == 16762)
    stopifnot(flows.matrix['MO', 'CA'] == 645)
    flows.matrix
  }
flows.matrix <- GetShipmentFlows()

GetStateDists <- function(){
  x <- do.call(cbind, state.center)
  M <- spDists(x, longlat=TRUE)
  colnames(M) <- rownames(M) <- state.abb
  M
}
state.to.state.dists <- GetStateDists()

save(real.case.data, shared.border.adjacency, flows.matrix, state.to.state.dists,
     file='common-data.RData')
