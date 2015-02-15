#!/usr/bin/Rscript

library(BatchExperiments)

load('corSims-layer0.RData')

pdes <- do.call(rbind, list(data.frame(comSize=2^4, nedges=2^8),
                            data.frame(comSize=2^5, nedges=2^10)))
prob.design <- makeDesign('prob', design=pdes,
                 exhaustive=list(couplingFactor=c(1), scheme=c('degreeBalanced')))
alg.design <- makeDesign('sim', exhaustive=list(R0baseline=c(1), expectedIntros=c(1)))
addExperiments(reg, prob.designs=prob.design, algo.designs=alg.design, repls=repls)

summarizeExperiments(reg)

newJobs <- findExperiments(reg, prob.pars=(comSize!=64))
chunked <- chunk(newJobs, n.chunk=20, shuffle=TRUE)

submitJobs(reg, chunked)

tmpf <- function(job, res){
    list(cc1=res$crossCor, cohesion=res$cohesion, abConnectivity=res$abConnectivity)
}
rdf <- reduceResultsExperiments(reg, fun=tmpf)
save.image(file='corSims-layer1.RData')
q('no')

