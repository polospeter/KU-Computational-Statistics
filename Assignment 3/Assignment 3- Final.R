
cores <- detectCores()
cluster <- makePSOCKcluster(cores)
system.time(parLapply(cluster, 1:10, function(i) Sys.sleep(i)))

clusterEvalQ()

x <- 10
psock <- parallel::makePSOCKcluster(1L)
clusterEvalQ(psock, x)

clusterEvalQ()

# Parallel sampling  ------------------------------------------------------

library(parallel)
detectCores()  ## Is eight on my laptop
## Setting up a cluster with eight nodes
cl <- makeCluster(8)

## Then you have to copy relevant data and functions to the nodes
clusterExport(cl, c("kernDensEpan3","simul.data"))

## Correct generation of independent (pseudo)random number sequences in parallel 
## processes is a bit tricky. Using the L'Ecuyer et al. generator makes it possible to control 
## this process and makes it reproducible. This is possible by distributing a correct initialization
## of the random number generator to the nodes.
clusterSetRNGStream(cl, iseed = 10)

kern.parallel <- parLapply(
  cl,
  1:8,   ## Eight because I have eight cores
  function(n) kernDensEpan3(simul.data,0.25*sqrt(5))
)


