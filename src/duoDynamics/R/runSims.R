runSims <-
function(A, R0baseline=0.5, ngens=30, expectedIntros=0.5){

    getNextState <- function(A, state, beta){
        forces <- as.numeric(A%*% state)
        transmissionProbs <- 1 - exp(-forces * beta)
        rbinom(state, 1, transmissionProbs)
    }

    states <- list()
    N <- ncol(A)
    pintro <- expectedIntros / N
    beta <- R0baseline / (N - 1)
    states[[1]] <- rep(0, length=N)    

    for(i in seq_len(ngens)){
        ns <- getNextState(A=A, state=states[[i]], beta=beta)
        ind <- rbinom(ns, size=1, p=pintro)
        ns[ind==1] <- 1
        states[[i + 1]] <- ns
    }
    states
}
