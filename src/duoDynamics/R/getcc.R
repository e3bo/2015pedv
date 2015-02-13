getcc <-
function(x, lag=1, burnin=1:500){
    foo <- ccf(x[1,-burnin], x[2,-burnin], plot=FALSE)
    ind <- which(foo$lag == lag)
    foo$acf[ind]
}
