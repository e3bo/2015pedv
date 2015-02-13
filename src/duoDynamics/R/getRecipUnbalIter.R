getRecipUnbalIter <-
function(na, nb){
    ends <- c(na, nb)
    runPosition <- 0
    runCount <- 0
    ord <- 0
    maxInd <- ends[ord + 1]
    runLen <- maxInd - runCount
    function(){
        i <- maxInd - runPosition
        j <- runCount %/% 2 + 1
        if(ord== 0){
            ret <- c(i,j)
        } else {
            ret <- c(j,i)
        }
        if (runPosition + 1 < runLen){
            runPosition <<- runPosition + 1
        } else {
            runCount <<- runCount + 1
            runPosition <<- 0
            ord <<- ord*-1 + 1
            if(ord==1){
                runLen <<- runLen - 1
            }
            maxInd <<- ends[ord + 1]        
        }
        ret
    }    
}
