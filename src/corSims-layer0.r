#!/usr/bin/Rscript

library(BatchExperiments)

reg <- makeExperimentRegistry(id="corSims", packages=c("igraph", "duoDynamics"))

repls <- 20
static <- data.frame(ngens=1000)
addProblem(reg, id='prob', static=static, dynamic=generateGraph)
addAlgorithm(reg, id='sim', fun=simAndSummarize)

pdes <- do.call(rbind,
                list(data.frame(comSize=16, nedges=2^(0:8)),
                     data.frame(comSize=32, nedges=2^(0:10))))
#                     data.frame(comSize=64, nedges=2^(0:12))))
prob.design <- makeDesign('prob', design=pdes,
 exhaustive=list(
     couplingFactor=c(0.5,1,2),
     scheme=c('degreeBalanced', 'unbalanced', 'recipUnbalanced')))

alg.design <- makeDesign('sim',
                         exhaustive=list(
                             R0baseline=c(0.5,1,2),
                             expectedIntros=c(0.01,1,10)))

addExperiments(reg, prob.designs=prob.design, algo.designs=alg.design, repls=repls)

summarizeExperiments(reg)

#chunked <- chunk(getJobIds(reg), chunk.size=20, shuffle=TRUE)

submitJobs(reg, 53000:53002)

tmpf <- function(job, res){
    list(cc1=res$crossCor, cohesion=res$cohesion, abConnectivity=res$abConnectivity)
}
rdf <- reduceResultsExperiments(reg, fun=tmpf)
save.image(file='corSims-layer0.RData')
q('no')

