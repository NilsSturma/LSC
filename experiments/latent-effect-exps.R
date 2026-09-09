source("R/semi-direct-id.R", chdir = TRUE)
source("R/direct-id.R", chdir = TRUE)
source("R/L2O-IV-id.R", chdir = TRUE)
library(rjson)
library(foreach)
library(doParallel)
library(doRNG)

nObs = 10
nLat = 3
nInd = 2        # indicators per latent node
pList = seq(0.02, 0.20, 0.02)
ngraphs = 1000
seed = 100
nCores = 7


# Random acyclic latent digraph in which every latent node has nInd indicators
# and the latent nodes form a chain, so that there always are latent-to-latent
# effects to identify. Every other edge is drawn independently with probability
# p, in particular the edges from a latent node to an indicator of another
# latent node, so that the indicators are not pure by construction.
#
# Acyclicity is ensured by drawing a level for every node and orienting all
# edges from the lower to the higher level. An indicator is drawn after its own
# latent node, all other nodes are placed uniformly at random, so that free
# observed nodes and indicators can appear anywhere in the topological order.
rIndicatorAdjMatrix <- function(nObs, nLat, nInd, p) {
  nTot = nObs + nLat
  observedNodes = seq(nObs)
  latentNodes = (nObs+1):nTot
  indicators = split(seq(nInd*nLat), rep(latentNodes, each = nInd))

  levels = numeric(nTot)
  levels[observedNodes] = runif(nObs)
  levels[latentNodes] = sort(runif(nLat))
  for (h in latentNodes) {
    levels[indicators[[as.character(h)]]] = runif(nInd, levels[h], 1)
  }

  # All pairs consistent with the levels are candidate edges
  candidates = outer(levels, levels, "<")

  L = 1 * (candidates & matrix(runif(nTot^2) < p, nTot, nTot))

  # The indicators and the chain of latent nodes are always present
  for (h in latentNodes) {
    L[h, indicators[[as.character(h)]]] = 1
  }
  for (i in seq(nLat-1)) {
    L[latentNodes[i], latentNodes[i+1]] = 1
  }

  return(L)
}


# Number of pure children of every latent node, needed to see how often the
# methods are applicable at all. identifyL2OIV needs at least one pure
# observed child per latent node, whereas directID needs at least two pure
# children, which may also be latent nodes.
nPureObsChildren <- function(g, latentNodes) {
  return(vapply(latentNodes, function(h) {
    length(validScalingIndicators(g, h))
  }, integer(1)))
}

nPureAnyChildren <- function(g, latentNodes) {
  return(vapply(latentNodes, function(h) {
    length(pureChildren(g, h))
  }, integer(1)))
}


#############################################
# Generate Graphs and check Identifiability #
#############################################

# Both methods are applied in two ways:
#
#   1. Directly to the original graph.
#   2. To the latent subgraph, that is the graph without the edges outgoing
#      from observed nodes, which is what latentCovGraph returns. This is only
#      justified if the semi-direct effects are identified beforehand, which is
#      checked by the latent subgraph criterion.

observedNodes = seq(nObs)
latentNodes = (nObs+1):(nObs+nLat)

# One task per graph, ordered by decreasing p, since the methods take longer on
# denser graphs (longest job first)
tasks = expand.grid(graph = 1:ngraphs, p = pList)
tasks = tasks[order(-tasks$p), ]

cl <- makeCluster(nCores, outfile = "")
registerDoParallel(cl)
set.seed(seed)

results <- foreach(t = 1:nrow(tasks),
                   .combine = 'c',
                   .multicombine=TRUE,
                   .errorhandling="remove",
                   .packages=c("igraph", "sna", "SEMID", "lpSolve")) %dorng% {

p = tasks$p[t]

if ((t%%500)==0){
  print(t)
}

obj <- tryCatch(
  {
    L <- rIndicatorAdjMatrix(nObs, nLat, nInd, p)
    g = LatentDigraph(L, observedNodes, latentNodes)
    gLatent = latentCovGraph(g)

    list("p"=p,
         "lsc"=LSCID(g, subsetSizeControl=Inf)$id,
         "directOrig"=directID(g)$id,
         "l2oOrig"=identifyL2OIV(g)$allIdentified,
         "directLatent"=directID(gLatent)$id,
         "l2oLatent"=identifyL2OIV(gLatent)$allIdentified,
         "pureObsOrig"=nPureObsChildren(g, latentNodes),
         "pureAnyOrig"=nPureAnyChildren(g, latentNodes),
         "pureObsLatent"=nPureObsChildren(gLatent, latentNodes),
         "pureAnyLatent"=nPureAnyChildren(gLatent, latentNodes),
         "adjMatrix"=c(t(L)))  # rowwise
  },
  error = function(e){
    print(e$message)
    list("p"=p,
         "lsc"=NA, "directOrig"=NA, "l2oOrig"=NA,
         "directLatent"=NA, "l2oLatent"=NA,
         "pureObsOrig"=NA, "pureAnyOrig"=NA,
         "pureObsLatent"=NA, "pureAnyLatent"=NA, "adjMatrix"=NA)
  }
)

list(obj)
}

stopCluster(cl)


########
# Save #
########

print(length(results))

jsonList = list()
for (k in 1:length(results)){
  newObj = list(results[[k]])
  names(newObj) <- k
  jsonList = c(jsonList, newObj)
}

name=paste("experiments/O", nObs, "L", nLat, "-latent.json", sep="")
jsonData = toJSON(jsonList)
write(jsonData, name)


######################
# Compute Statistics #
######################

table = matrix(0,length(pList),12)
rowMatching = as.list(1:length(pList))
names(rowMatching) = pList

for (k in 1:length(results)){
  obj = results[[k]]
  row = rowMatching[[as.character(obj$p)]]
  if (any(is.na(obj$lsc))){
    next
  }

  table[row, 1] <- table[row, 1]+1
  if (obj$lsc){
    table[row, 2] <- table[row, 2]+1
  }

  # Methods applied to the original graph
  if (obj$directOrig){
    table[row, 3] <- table[row, 3]+1
  }
  if (obj$l2oOrig){
    table[row, 4] <- table[row, 4]+1
  }

  # Latent subgraph criterion first, then the methods on the latent subgraph
  if (obj$lsc && obj$directLatent){
    table[row, 5] <- table[row, 5]+1
  }
  if (obj$lsc && obj$l2oLatent){
    table[row, 6] <- table[row, 6]+1
  }

  # Methods on the latent subgraph without the latent subgraph criterion,
  # reported to see which of the two steps is binding
  if (obj$directLatent){
    table[row, 7] <- table[row, 7]+1
  }
  if (obj$l2oLatent){
    table[row, 8] <- table[row, 8]+1
  }

  # How often the methods are applicable at all, i.e. how often every latent
  # node has the pure children the methods need
  if (all(obj$pureAnyOrig >= 2)){
    table[row, 9] <- table[row, 9]+1
  }
  if (all(obj$pureObsOrig >= 1)){
    table[row, 10] <- table[row, 10]+1
  }
  if (all(obj$pureAnyLatent >= 2)){
    table[row, 11] <- table[row, 11]+1
  }
  if (all(obj$pureObsLatent >= 1)){
    table[row, 12] <- table[row, 12]+1
  }
}

colnames(table) <- c("nGraphs", "nLSC",
                     "nDirectOrig", "nL2OOrig",
                     "nLSCDirect", "nLSCL2O",
                     "nDirectLatent", "nL2OLatent",
                     "nPure2Orig", "nPure1Orig",
                     "nPure2Latent", "nPure1Latent")
rownames(table) <- pList

print(table)
print(colSums(table))
write.table(table, file = paste("experiments/O", nObs, "L", nLat, "-latent.txt", sep=""),
            sep = "\t", row.names = TRUE, col.names = TRUE)
