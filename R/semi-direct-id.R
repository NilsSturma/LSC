library(SEMID)
source("joint-flow.R")

################################
##### Function Definitions #####
################################

subsetsOfSize <- function(x, k) {
  if (k > length(x)) {
    return(list())
  }
  if (k == 0) {
    return(list(numeric(0)))
  }
  if (length(x) == k) {
    return(list(x))
  }
  return(combn(x, k, simplify = F))
}

semiDirectEffectGraph <- function(g){
  nObs <- g$numObserved()
  nLat <- g$numLatents()
  L <- g$.L
  nTot <- nrow(L)
  semiDirectAdjMat <- L[1:nObs,1:nObs] +  
    L[1:nObs,(nObs+1):nTot] %*% 
    solve(diag(nLat) - L[(nObs+1):nTot,(nObs+1):nTot]) %*% 
    L[(nObs+1):nTot,1:nObs]
  semiDirectAdjMat <- 1*(semiDirectAdjMat > 0)
  gSemiDirect <-  LatentDigraph(semiDirectAdjMat, seq(nObs), c())
  return(gSemiDirect)
}


latentCovGraph <- function(g){
  observedNodes <- g$.observedNodes
  latentNodes <- g$.latentNodes
  
  latentL <- g$.L
  latentL[observedNodes, c(observedNodes, latentNodes)] <- 0
  latentGraph <-  LatentDigraph(latentL, observedNodes, latentNodes)
  return(latentGraph)
}


allowedNodesForZ <- function(g, v, S, H1, H2){
  
  semiParentsOfV <- g$.semiDirect$parents(v)
  
  trFromH1 <- g$.latentCovGraph$trFrom(H1, includeLatents = FALSE)
  desH2 <- g$.latentCovGraph$descendants(H2, includeLatents = FALSE)
  
  allowedNodes <- intersect(union(trFromH1, desH2), S)
  allowedNodes <- setdiff(allowedNodes, union(v, semiParentsOfV))
  
  return(allowedNodes)
}


allowedNodesForY <- function(g, v, S, Z, H1, H2){
  latentTrFromZandVAvoidingH2H1 <- g$.latentCovGraph$trFrom(c(Z,v), 
                                         includeLatents = FALSE, 
                                         avoidLeftNodes = H2,
                                         avoidRightNodes = H1)
  # note: this set contains Z and v
  extLatentTrFromZandVAvoidingH2H1 <- intersect(g$.observedNodes, 
                                        g$descendants(latentTrFromZandVAvoidingH2H1))
  notAllowed = setdiff(extLatentTrFromZandVAvoidingH2H1, S)
  notAllowed = union(notAllowed, latentTrFromZandVAvoidingH2H1)
  allowed <- setdiff(g$.observedNodes, notAllowed)
  return(allowed)
}


getLpGraphs <- function(g){
  observedNodes <- g$.observedNodes
  latentNodes <- g$.latentNodes
  nObs <- g$numObserved()
  nLat <- g$numLatents()
  m = nObs + nLat
  
  # Initialize adjacency matrices
  adjMat <- matrix(0, 2 * m, 2 * m)
  adjMatSub <- matrix(0, 2 * m, 2 * m)
  
  # If i -> j with i latent then left j points to left i (opposite direction)
  latentL <- g$.L
  latentL[observedNodes, c(observedNodes, latentNodes)] <- 0
  adjMat[1:m, 1:m] <- t(latentL)
  adjMatSub[1:m, 1:m] <- t(latentL)
  
  # Left nodes point to their corresponding right nodes
  adjMat[cbind(1:m, m + 1:m)] <- 1
  adjMatSub[cbind(1:m, m + 1:m)] <- 1
  
  # If i -> j then right i points to right j
  adjMat[m + 1:m, m + 1:m] <- g$.L
  adjMatSub[m + 1:m, m + 1:m] <- latentL
  
  gLp <- igraph::graph_from_adjacency_matrix(adjMat, mode="directed")
  gLpSub <- igraph::graph_from_adjacency_matrix(adjMatSub, mode="directed")
  
  return(list(gLp,gLpSub))
}


constructTrekSystem <- function(res, flowGraph, flowSubGraph, s, t, m){
  pathSystem <- constructPathSystem(res, flowGraph, flowSubGraph, s, t)
  TrekSystem <- list()
  startNodes <- c()
  for (i in 1:length(pathSystem)){
    path <- pathSystem[[i]]
    startNodes[i] <- path[2,1]
    path <- path[-1,]
    path <- path[-nrow(path),]
    path[] <- vapply(path, function(x){if(x>m) x-m else x}, numeric(1))
    if (is.vector(path)){
      path <- matrix(path, nrow=1)
    } else {
      rows_to_keep <- apply(path, 1, function(row) length(unique(row)) > 1)
      path <- path[rows_to_keep, ]
    }
    TrekSystem[[i]] <- path
  }
  return(list("TrekSystem"=TrekSystem, "startNodes"=startNodes))
}


checkTrekSystem <- function(g, Z, v, Ya){
  
  # Define variables
  semiParentsOfV <- g$.semiDirect$parents(v)
  observedNodes <- g$.observedNodes
  latentNodes <- g$.latentNodes
  nObs <- g$numObserved()
  nLat <- g$numLatents()
  m = nObs + nLat
  
  # Create adjacency matrices
  flowAdjMat <- matrix(0, 2*m+2, 2*m+2)
  flowAdjMatSub <- matrix(0, 2*m+2, 2*m+2)
  s = 2*m+1
  t = 2*m+2
  
  flowAdjMat[1:(2*m), 1:(2*m)] <- as_adjacency_matrix(g$.LpGraph, sparse=FALSE)
  flowAdjMat[s,Ya] <- 1
  flowAdjMat[m + c(semiParentsOfV, Z),t] <- 1
  
  flowAdjMatSub[1:(2*m), 1:(2*m)] <- as_adjacency_matrix(g$.LpSubGraph, sparse=FALSE)
  flowAdjMatSub[s,Ya] <- 1
  flowAdjMatSub[m + Z,t] <- 1
  
  # Define graphs
  flowGraph <- igraph::graph_from_adjacency_matrix(flowAdjMat, mode="directed")
  flowSubGraph <- igraph::graph_from_adjacency_matrix(flowAdjMatSub, mode="directed")

  # Run linear program
  res <- jointFlow(flowGraph, flowSubGraph, s, t)
  objval <- res$objval
  
  if (!all(res$solution - round(res$solution) == 0)){
    #print("Linear program did not return integer solution.")
    res <- jointFlow(flowGraph, flowSubGraph, s, t, int=TRUE)
    if ((res$objval != objval) && ((length(Z)+length(semiParentsOfV))==objval)){
      print("Maximal value of integer program differs to maximal value |Z|+|P| of linear program.")
    }
  } 
  objval <- res$objval
  if(objval==(length(Z)+length(semiParentsOfV))){
    TrekSystem <- constructTrekSystem(res, flowGraph, flowSubGraph, s, t, m)
    return(list("objval" = res$objval,  
                "trekSystem" = TrekSystem$TrekSystem, 
                "Y"=TrekSystem$startNodes))
  } else {
    return(list("objval" = res$objval))
  }
}

checkID <- function(g, subsetSizeControl=Inf){
  
  g$.semiDirect <- semiDirectEffectGraph(g)
  g$.latentCovGraph <-latentCovGraph(g)
  graphs <- getLpGraphs(g)
  g$.LpGraph <- graphs[[1]]
  g$.LpSubGraph <- graphs[[2]]
  observedNodes <- g$.observedNodes
  latentNodes <- g$.latentNodes
  nObs <- g$numObserved()
  nLat <- g$numLatents()
  
  semiDirectParents <- lapply(g$.observedNodes, g$.semiDirect$parents)
  S <- which(sapply(semiDirectParents, function(x) {length(x)==0}))
  
  # For saving results
  Ys <- rep(list(numeric(0)), nObs)
  Zs <- rep(list(numeric(0)), nObs)
  H1s <- rep(list(numeric(0)), nObs)
  H2s <- rep(list(numeric(0)), nObs)
  trekSystems <- rep(list(numeric(0)), nObs)
  
  if (length(S)!=length(observedNodes)){
    changeFlag <- TRUE
  } else {
    changeFlag = FALSE
  }
  
  while(changeFlag){
    
    changeFlag=FALSE
    
    # Loop over all unsolved nodes
    for (v in setdiff(observedNodes, S)) {
      
      # Collect basic info of unsolved node v
      semiParentsOfV <- g$.semiDirect$parents(v)
      
      maxK <- min(subsetSizeControl, 
                  nLat, 
                  floor((nObs - 1 - length(semiParentsOfV)) / 2))
      
      # Loop over possible cardinalities of |H1|+|H2|
      for (k in seq(0, length = 1 + maxK)){
        
        # Loop over all sets H1 and H2 such that |H1|+|H2|=k
        for (l in 0:k){
          for (H1 in subsetsOfSize(latentNodes, l)){
            for (H2 in subsetsOfSize(latentNodes, k-l)){
              
              
              # Loop over all possible sets Z
              Za = allowedNodesForZ(g, v, S, H1, H2)
              for (Z in subsetsOfSize(Za, k)){
                
                # Define the set of allowed nodes for Y
                Ya <- allowedNodesForY(g, v, S, Z, H1, H2)
                if (length(Ya) >= (length(semiParentsOfV) + length(Z))){
                  res <- checkTrekSystem(g, Z, v, Ya)
                  
                  # If trek system exists, we know that v is identified
                  if(res$objval==(length(semiParentsOfV) + length(Z))){
                    Ys[[v]] <- res$Y
                    Zs[[v]] <- Z
                    H1s[[v]] <- H1
                    H2s[[v]] <- H2
                    trekSystems[[v]] <- res$trekSystem
                    S <- c(v, S)
                    if (length(S)!=length(observedNodes)){
                      changeFlag <- TRUE
                    } else {
                      changeFlag = FALSE
                    }
                    break
                  }
                }
              }
              if (v %in% S){break}
            }
            if (v %in% S){break}
          }
          if (v %in% S){break}
        }
        if (v %in% S){break}
      }
      if (v %in% S){break}
    }
  }
  if (length(S)==length(observedNodes)){
    identifiable = TRUE
  } else {
    identifiable = FALSE
  }
  result <- list("S"=S,
                 "Ys"=Ys,
                 "Zs"=Zs,
                 "H1s"=H1s,
                 "H2s"=H2s,
                 "trekSystems"=trekSystems,
                 "id"=identifiable)
  return(result)
}
