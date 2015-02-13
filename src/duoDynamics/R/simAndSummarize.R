simAndSummarize <-
function(job, static, dynamic, ...){
    states <- runSims(A=dynamic$A, ngens=static$ngens, ...)
    ret <- list()
    ret$tss <- states2ts(states=states, labs=dynamic$labs)
    ret$crossCor <- getcc(ret$tss)
    ret$cohesion <- dynamic$cohesion
    ret$abConnectivity <- dynamic$abConnectivity
    ret
}
