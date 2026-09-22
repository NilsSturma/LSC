source("R/semi-direct-id.R", chdir = TRUE)
source("R/canonicalization.R", chdir = TRUE)
library(sna)  # rperm
library(rjson)
library(foreach)
library(doParallel)
library(doRNG)

rAcyclicDirectedAdjMatrix <- function(n, p) {
  return(1 * (upper.tri(matrix(0, n, n)) & matrix(sample(c(T, F), n^2, replace = T,
                                                         prob = c(p, 1 - p)), ncol = n)))
}

nNodes=15
nLat=5
pErdosList = seq(0.15, 0.45, 0.05)
ngraphs = 1000
seed = 100
nCores = 7   # number of cores used for the parallelization


#############################################
# Generate Graphs and check Identifiability #
#############################################

observedNodes = seq(nNodes-nLat)
latentNodes = (nNodes-nLat+1):nNodes

tasks = expand.grid(graph = 1:ngraphs, pErdos = pErdosList)
tasks = tasks[order(-tasks$pErdos), ]

cl <- makeCluster(nCores, outfile = "")
registerDoParallel(cl)
set.seed(seed)

results <- foreach(t = 1:nrow(tasks),
                   .combine = 'c',
                   .multicombine=TRUE,
                   .errorhandling="remove",
                   .packages=c("igraph", "sna", "SEMID", "lpSolve")) %dorng% {

pErdos = tasks$pErdos[t]

if ((t%%100)==0){
  print(t)
}

obj <- tryCatch(
  {
    L <- rAcyclicDirectedAdjMatrix(nNodes,pErdos)
    P = rperm(rep(0,nNodes))
    L <- L[P,P]
    g = LatentDigraph(L, observedNodes, latentNodes)
    gCan <- canonicalization(g)

    tStart <- proc.time()[["elapsed"]]
    idRes <- LSCID(g, subsetSizeControl=Inf)
    time <- proc.time()[["elapsed"]] - tStart

    tStart <- proc.time()[["elapsed"]]
    idResCan <- LSCID(gCan, subsetSizeControl=Inf)
    timeCan <- proc.time()[["elapsed"]] - tStart

    list("g"=g, "pErdos"=pErdos,
         "res"=idRes, "resCan"=idResCan,
         "time"=time, "timeCan"=timeCan)
  },
  error = function(e){
    print(e$message)
    print(L)
    list("g"=NA, "pErdos"=pErdos,
         "res"=NA, "resCan"=NA,
         "time"=NA, "timeCan"=NA)
  }
)

list(obj)
}

stopCluster(cl)

########
# Save #
########

# Change format of list
jsonList = list()
for (k in 1:length(results)){
  oldObj = results[[k]]
  if (!identical(oldObj$g, NA)){
    adjMat = oldObj$g$L()
    newObj <- list(list("pErdos"=oldObj$pErdos,
                        "res"=oldObj$res, "resCan"=oldObj$resCan,
                        "time"=oldObj$time, "timeCan"=oldObj$timeCan,
                        "adjMatrix" = c(t(adjMat))))  # rowwise
  } else {
    newObj <- list(list("pErdos"=oldObj$pErdos,
                        "res"=oldObj$res, "resCan"=oldObj$resCan,
                        "time"=oldObj$time, "timeCan"=oldObj$timeCan,
                        "adjMatrix" = NA))
  }

  names(newObj) <- k
  jsonList = c(jsonList, newObj)
}

# Save as json
name=paste("experiments/O", nNodes-nLat, "L", nLat, ".json", sep="")
jsonData = toJSON(jsonList)
write(jsonData, name)



######################
# Compute Statistics #
######################

table = matrix(0,length(pErdosList),2)
rowMatching = as.list(1:length(pErdosList))
names(rowMatching) = pErdosList

timeSums = matrix(0,length(pErdosList),2)
timeCounts = matrix(0,length(pErdosList),2)
fracSums = matrix(0,length(pErdosList),2)
fracCounts = matrix(0,length(pErdosList),2)
diffSums = matrix(0,length(pErdosList),2)
nILPSums = matrix(0,length(pErdosList),2)

for (k in 1:length(results)){
  obj = results[[k]]
  row = rowMatching[[as.character(obj$pErdos)]]
  if (!any(is.na(obj$res))){
    if (obj$res$id){
      table[row, 1] <- table[row, 1]+1
    }
    nILPSums[row, 1] <- nILPSums[row, 1] + obj$res$nILP
    if (obj$res$nLP > 0){
      fracSums[row, 1] <- fracSums[row, 1] + obj$res$nILP/obj$res$nLP
      diffSums[row, 1] <- diffSums[row, 1] + obj$res$nILPdiff/obj$res$nLP
      fracCounts[row, 1] <- fracCounts[row, 1] + 1
    }
  }
  if (!any(is.na(obj$resCan))){
    if (obj$resCan$id){
      table[row, 2] <- table[row, 2]+1
    }
    nILPSums[row, 2] <- nILPSums[row, 2] + obj$resCan$nILP
    if (obj$resCan$nLP > 0){
      fracSums[row, 2] <- fracSums[row, 2] + obj$resCan$nILP/obj$resCan$nLP
      diffSums[row, 2] <- diffSums[row, 2] + obj$resCan$nILPdiff/obj$resCan$nLP
      fracCounts[row, 2] <- fracCounts[row, 2] + 1
    }
  }
  if (!is.na(obj$time)){
    timeSums[row, 1] <- timeSums[row, 1] + obj$time
    timeCounts[row, 1] <- timeCounts[row, 1] + 1
  }
  if (!is.na(obj$timeCan)){
    timeSums[row, 2] <- timeSums[row, 2] + obj$timeCan
    timeCounts[row, 2] <- timeCounts[row, 2] + 1
  }
}

avgTimes = timeSums / timeCounts
avgFracs = fracSums / fracCounts
avgDiffs = diffSums / fracCounts
table = cbind(table, round(avgTimes, 4), round(avgFracs, 4), round(avgDiffs, 4),
              nILPSums)

colnames(table) <- c("nLSC", "nCanLSC", "timeLSC", "timeCanLSC",
                     "fracILPLSC", "fracILPCanLSC",
                     "fracDiffLSC", "fracDiffCanLSC",
                     "nILPLSC", "nILPCanLSC")
rownames(table) <- pErdosList

print(table)
write.table(table, file = paste("experiments/O", nNodes-nLat, "L", nLat, ".txt", sep=""),
            sep = "\t", row.names = TRUE, col.names = TRUE)
