source("semi-direct-id.R")
source("canonicalization.R")
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


#############################################
# Generate Graphs and check Identifiability #
#############################################

observedNodes = seq(nNodes-nLat)
latentNodes = (nNodes-nLat+1):nNodes

cores = min(detectCores(), length(pErdosList))
cl <- makeCluster(cores, outfile = "")
registerDoParallel(cl)
set.seed(seed)

results <- foreach(k = 1:length(pErdosList),
                   .combine = 'c',
                   .multicombine=TRUE,
                   .errorhandling="remove",
                   .packages=c("igraph", "sna", "SEMID", "lpSolve")) %dorng% {

pErdos = pErdosList[k]
res <- list()

for (i in 1:ngraphs){
  if((i%%10)==0){
    print(i)
  }
  obj <- tryCatch(
    {
      L <- rAcyclicDirectedAdjMatrix(nNodes,pErdos)
      P = rperm(rep(0,nNodes))
      L <- L[P,P]
      g = LatentDigraph(L, observedNodes, latentNodes)
      gCan <- canonicalization(g)

      idRes1 <- checkID(g, subsetSizeControl=1)
      idResCan1 <- checkID(gCan, subsetSizeControl=1)
      idRes2 <- checkID(g, subsetSizeControl=2)
      idResCan2 <- checkID(gCan, subsetSizeControl=2)
      idRes3 <- checkID(g, subsetSizeControl=3)
      idResCan3 <- checkID(gCan, subsetSizeControl=3)

      list("g"=g, "pErdos"=pErdos,
           "res1"=idRes1, "resCan1"=idResCan1,
           "res2"=idRes2, "resCan2"=idResCan2,
           "res3"=idRes3, "resCan3"=idResCan3)
    },
    error = function(e){
      print(e$message)
      print(L)
      list("g"=NA, "pErdos"=pErdos,
           "res1"=NA, "resCan1"=NA,
           "res2"=NA, "resCan2"=NA,
           "res3"=NA, "resCan3"=NA)
    }
  )

  res[[i]] <- obj
  }
res
}

stopCluster(cl)

########
# Save #
########

print(length(results))
# Save as RData file
#name=paste("O", nNodes-nLat, "L", nLat, ".RData", sep="")
#saveRDS(results, file=name)

# Change format of list
jsonList = list()
for (k in 1:length(results)){
  oldObj = results[[k]]
  if (!identical(oldObj$g, NA)){
    adjMat = oldObj$g$L()
    newObj <- list(list("pErdos"=oldObj$pErdos,
                        "res1"=oldObj$res1, "resCan1"=oldObj$resCan1,
                        "res2"=oldObj$res2, "resCan2"=oldObj$resCan2,
                        "res3"=oldObj$res3, "resCan3"=oldObj$resCan3,
                        "adjMatrix" = c(t(adjMat))))  # rowwise
  } else {
    newObj <- list(list("pErdos"=oldObj$pErdos,
                        "res1"=oldObj$res1, "resCan1"=oldObj$resCan1,
                        "res2"=oldObj$res2, "resCan2"=oldObj$resCan2,
                        "res3"=oldObj$res3, "resCan3"=oldObj$resCan3,
                        "adjMatrix" = NA))
  }

  names(newObj) <- k
  jsonList = c(jsonList, newObj)
}

# Save as json
name=paste("O", nNodes-nLat, "L", nLat, ".json", sep="")
jsonData = toJSON(jsonList)
write(jsonData, name)



######################
# Compute Statistics #
######################

table = matrix(0,length(pErdosList),6)
rowMatching = as.list(1:length(pErdosList))
names(rowMatching) = pErdosList

for (k in 1:length(results)){
  obj = results[[k]]
  row = rowMatching[[as.character(obj$pErdos)]]
  if (!any(is.na(obj$res1))){
    if (obj$res1$id){
      table[row, 1] <- table[row, 1]+1
    }
  }
  if (!any(is.na(obj$resCan1))){
    if (obj$resCan1$id){
      table[row, 2] <- table[row, 2]+1
    }
  }
  if (!any(is.na(obj$res2))){
    if (obj$res2$id){
      table[row, 3] <- table[row, 3]+1
    }
  }
  if (!any(is.na(obj$resCan2))){
    if (obj$resCan2$id){
      table[row, 4] <- table[row, 4]+1
    }
  }
  if (!any(is.na(obj$res3))){
    if (obj$res3$id){
      table[row, 5] <- table[row, 5]+1
    }
  }
  if (!any(is.na(obj$resCan3))){
    if (obj$resCan3$id){
      table[row, 6] <- table[row, 6]+1
    }
  }
}

colnames(table) <- c("nLSC1", "nCanLSC1", "nLSC2", "nCanLSC2", "nLSC3", "nCanLSC3")
rownames(table) <- pErdosList

print(table)
write.table(table, file = paste("O", nNodes-nLat, "L", nLat, ".txt", sep=""),
            sep = "\t", row.names = TRUE, col.names = TRUE)
