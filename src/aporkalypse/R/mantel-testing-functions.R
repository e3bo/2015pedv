GetCrossCorrs <- function(lag, obs){
    n <- ncol(obs)
    CC <- matrix(nrow=n, ncol=n)
    getcc <- function(x,y, lag){
        foo <- ccf(x, y, plot=FALSE)
        ind <- which(foo$lag == lag)
        foo$acf[ind]
    }
    for(i in seq_len(n)){
        for(j in seq_len(n)){
            ## CC[i,j] will be high if deviations from the mean in series i
            ## are shifted to the left of similar deviations to the mean in series j
            ## i.e. i's deviations are indicative of j's future deviations
            CC[i, j] <- getcc(obs[, j], obs[, i], lag=lag)
        }
    }
    colnames(CC) <- rownames(CC) <- colnames(obs)
    CC
}

MakeDesign <- function(mat.list1, mat.list2, mat.list3=NULL,
                       methods=c('spearman', 'pearson'),
                       symmetrize=c(TRUE, FALSE),
                       permutations=10000){
  if (is.null(mat.list3)){
    des <- expand.grid(mat1.name=names(mat.list1),
                       mat2.name=names(mat.list2),
                       method=methods, symmetrize=symmetrize,
                       stringsAsFactors=FALSE)
  } else {
    des <- expand.grid(mat1.name=names(mat.list1),
                       mat2.name=names(mat.list2),
                       mat3.name=names(mat.list3),
                       method=methods, symmetrize=symmetrize,
                       stringsAsFactors=FALSE)
  }
  des <- des[des$mat1.name != des$mat2.name, ]
  des$permutations <- permutations
  des
}

DoDesignedTests <- function(des, mat.list1, mat.list2){
  DoTest <- function(x, y, symmetrize=FALSE, ...){
    if(symmetrize){
        x <- (x + t(x)) * 0.5
        y <- (y + t(y)) * 0.5
    }
    vegan::mantel(x, y, ...)
  }
  res <- mapply(DoTest, x=mat.list1[des$mat1.name],
                y=mat.list2[des$mat2.name],
                symmetrize=des$symmetrize,
                permutations=des$permutations,
                method=des$method, SIMPLIFY=FALSE)
  des$r <- sapply(res, '[[', 'statistic')
  des$pValues <- sapply(res, '[[', 'signif')
  des
}

DoDesignedPartialTests <- function(des, mat.list1, mat.list2, mat.list3){
  DoPartialTest <- function(x, y, z, symmetrize=FALSE, ...){
    if(symmetrize){
        x <- (x + t(x)) * 0.5
        y <- (y + t(y)) * 0.5
        z <- (z + t(z)) * 0.5
    }
    vegan::mantel.partial(x, y, z, ...)
  }
  res <- mapply(DoPartialTest, x=mat.list1[des$mat1.name],
                y=mat.list2[des$mat2.name],
                z=mat.list3[des$mat3.name],
                symmetrize=des$symmetrize,
                permutations=des$permutations,
                method=des$method, SIMPLIFY=FALSE)
  des$r <- sapply(res, '[[', 'statistic')
  des$pValues <- sapply(res, '[[', 'signif')
  des
}

DoMantelTests <- function(mat.list1, mat.list2, ...){
  des <- MakeDesign(mat.list1, mat.list2, ...)
  DoDesignedTests(des, mat.list1, mat.list2)
}

DoPartialMantelTests <- function(mat.list1, mat.list2, mat.list3, ...){
  des <- MakeDesign(mat.list1, mat.list2, mat.list3, ...)
  DoDesignedPartialTests(des, mat.list1, mat.list2, mat.list3)
}

MakePopStructMats <- function(observed){
  nms <- colnames(observed)
  nhood <- shared.border.adjacency[nms, nms]
  ep <- flows.matrix[nms, nms]
  epl <- log10(ep + 1)
  centerDists <- state.to.state.dists[nms, nms]
  list('shipment'=epl, 'gcd'=-centerDists, 'sharedBord'=nhood)
}

MakePopDynMats <- function(observed, lags.sel=1){
  CC <- lapply(lags.sel, GetCrossCorrs, obs=observed)
  names(CC) <- paste0('lag', lags.sel)
  CC
}
