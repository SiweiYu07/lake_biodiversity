#### beta functional dissimilarities


comm.test<-read.csv("CHabundance.csv",row.names = 1)
traits.test<-read.csv("CHtrait.csv",row.names = 1)

library(betapart)
test.pair<-functional.beta.pair(x=comm.test, traits=traits.test, index.family = "jaccard")
lapply(test.pair,round,3)

#### with functional.betapart.core.pairwise
test1 <- functional.betapart.core.pairwise(comm.test, traits.test)
test.pair <- functional.beta.pair(test1)

## Not run:
#### if internal parallelisation would be interesting (large community matrix)
test1 <- functional.betapart.core.pairwise(comm.test, traits.test, parallel = TRUE,
                                           opt.parallel = list(nc = 5))
test.pair <- functional.beta.pair(test1)
test.pair

test.multi<-functional.beta.multi(x=comm.test, traits=traits.test, index.family = "jaccard" )
test.multi
## End(Not run)
