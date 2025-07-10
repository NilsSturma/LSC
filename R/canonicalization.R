library(SEMID)

canonicalization <- function(g){
  nObs <- g$numObserved()
  nLat <- g$numLatents()
  L <- g$.L
  nTot <- nrow(L)
  
  semiDirectAdjMat <- L[1:nObs,1:nObs] +  
    L[1:nObs,(nObs+1):nTot] %*% 
    solve(diag(nLat) - L[(nObs+1):nTot,(nObs+1):nTot]) %*% 
    L[(nObs+1):nTot,1:nObs]
  semiDirectAdjMat <- 1*(semiDirectAdjMat > 0)
  
  observedNodes <- g$.observedNodes
  latentNodes <- g$.latentNodes
  latentL <- g$.L
  latentL[observedNodes, c(observedNodes, latentNodes)] <- 0
  latentPaths <- solve(diag(nTot) - latentL)
  latentPaths <- 1*(latentPaths > 0)
  
  canL <- matrix(0, nTot,nTot)
  canL[observedNodes,observedNodes] <- semiDirectAdjMat
  canL[latentNodes, observedNodes] <- latentPaths[latentNodes, observedNodes]
  
  gCan <-  LatentDigraph(canL, observedNodes, latentNodes)
  return(gCan)
}
