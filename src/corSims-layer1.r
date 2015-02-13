#!/usr/bin/Rscript

library(BatchExperiments)

load('corSims-layer0.RData')

pdes <- data.frame(comSize=64, nedges=2^(0:12))
prob.design <- makeDesign('prob', design=pdes,
 exhaustive=list(
     couplingFactor=c(1),
     scheme=c('degreeBalanced', 'unbalanced', 'recipUnbalanced')))

alg.design <- makeDesign('sim',
                         exhaustive=list(
                             R0baseline=c(1),
                             expectedIntros=c(1)))

addExperiments(reg, prob.designs=prob.design, algo.designs=alg.design, repls=repls)

summarizeExperiments(reg)

newJobs <- findExperiments(reg, prob.pars=(comSize==64))
chunked <- chunk(newJobs, n.chunk=1, shuffle=TRUE)

submitJobs(reg, chunked)

tmpf <- function(job, res){
    list(cc1=res$crossCor, cohesion=res$cohesion, abConnectivity=res$abConnectivity)
}
rdf <- reduceResultsExperiments(reg, fun=tmpf)
save.image(file='corSims-layer1.RData')
q('no')

