getdist <-
function(x, burnin=1){
    mean((x[1, -burnin] - x[2,-burnin])^2)
}
