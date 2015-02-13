generateGraph <-
function(job, static, comSize, ...){
    g1 <- graph.full(comSize)
    g2 <- graph.full(comSize)
    G <- graph.disjoint.union(g1,g2)
    ret <- list()
    ret$labs <- clusters(G)$membership
    aids <- which(ret$labs == 1)
    bids <- which(ret$labs == 2)
    foo <- addCoupling(G=G, aids=aids, bids=bids, ...)
    ret$A <- foo$G[,]
    ret$abConnectivity <- foo$abConnectivity
    ret$cohesion <- graph.cohesion(foo$G)
    ret
}
