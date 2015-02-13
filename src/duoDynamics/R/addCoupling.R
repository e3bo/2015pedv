addCoupling <-
function(G, aids, bids, nedges=1, couplingFactor=1,
                        weightBalanced=TRUE, scheme='degreeBalanced'){
    E(G)$weight <- rep(1, times=ecount(G))
    if(nedges == 0){
        return(G)
    }
    na <- length(aids)
    nb <- length(bids)
    maxEdges <- na*nb
    totalStrength <- maxEdges * couplingFactor
    stopifnot(nedges <= maxEdges)
    if(weightBalanced){
        weight <- totalStrength / nedges
    }else{
        weight <- 1
    }
    ret <- list()
    if(scheme == 'degreeBalanced'){
        for(n in seq_len(nedges)){
            i <- n %% na + 1
            shift <- (n - 1) %/% nb 
            j <- (shift + n) %% nb + 1
            G[aids[i], bids[j]] <- weight
        }
        ret$abConnectivity <- min(n, na, nb)
    } else if (scheme == 'unbalanced') {
        for(n in seq_len(nedges)){
            i <- (n - 1) %/% na + 1
            j <- n %% na + 1
            G[aids[i], bids[j]] <- weight
        }
        ret$abConnectivity <- i
    } else {
        iter <- getRecipUnbalIter(na,nb)
        for(n in seq_len(nedges)){
            vals <- iter()
            i <- vals[1]; j <- vals[2]
            G[aids[i], bids[j]] <- weight
        }
        ret$abConnectivity <- min(i,j)
    }
    ret$G <- G
    ret
}
