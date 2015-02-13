states2ts <-
function(states, nfreq=1, labs, normalize=FALSE){
    sm <- do.call('rbind', states)
    smts <- ts(sm)
    foo <- aggregate(smts, nfrequency=nfreq)
    splt <- split.data.frame(t(foo), f=labs)
    tsbyblock <- lapply(splt, colSums)
    #tsbyblock <- Filter(function(x) !all(x==0), tsbyblock)
    if(normalize){
        tsbyblock <- lapply(tsbyblock, function(x) (x - mean(x))/sd(x))
    }
    ret <- do.call(rbind, tsbyblock)
    ret
}
