library(igraph)
library(lpSolve)

################################
##### Function Definitions #####
################################

isSubGraph <- function(g, gSub){
  adjMat <- as_adjacency_matrix(g)
  adjMatSub <- as_adjacency_matrix(gSub)
  if (all(dim(adjMat)==dim(adjMatSub))){
    return(all((adjMat-adjMatSub)>=0))
  } else {
    return(FALSE)
  }
}


getIncidenceMat <- function(g){
  edges <- as_edgelist(g)
  numVertices <- vcount(g)
  numEdges <- ecount(g)
  
  incMat <- matrix(0, numVertices, numEdges)
  for (e in 1:numEdges) {
    tail <- edges[e, 1]
    head <- edges[e, 2]
    
    incMat[tail, e] <- -1 
    incMat[head, e] <- 1
  }
  return(incMat)
}


getIncomingMat <- function(g){
  edges <- as_edgelist(g)
  numVertices <- vcount(g)
  numEdges <- ecount(g)
  
  incomingMat <- matrix(0, numVertices, numEdges)
  for (e in 1:numEdges) {
    tail <- edges[e, 1]
    head <- edges[e, 2]
    
    incomingMat[head, e] <- 1
  }
  return(incomingMat)
}


jointFlow <- function(g, gSub, s, t, int=FALSE){
  
  # Check
  if (!isSubGraph(g, gSub)){
    print("ERROR, gSub is no a subgraph of g")
    return(FALSE)
  }
  
  # Define variables
  numVertices <- vcount(g)
  numEdgesG <- ecount(g)
  numEdgesGSub <- ecount(gSub)
  edgesG <- as_edgelist(g)
  edgesGSub <- as_edgelist(gSub)
  
  # Objective
  obj <- rep(0, numEdgesG+numEdgesGSub)
  for (e in 1:numEdgesG){
    edge <- edgesG[e,]
    if (edge[2]==t){
      edgeInSubgraph <- apply(edgesGSub, 1, function(row) all(row == edge))
      if (any(edgeInSubgraph)){
        obj[numEdgesG+which(edgeInSubgraph)] <- 1
      } else {
        obj[e] <- 1
      }
    }
  }
  
  # Equality constraints ("flow conservation")
  numEq <- 2*numVertices-4
  matEq <- matrix(0, numEq, numEdgesG + numEdgesGSub)
  matEq[1:(numVertices-2),1:numEdgesG] <- getIncidenceMat(g)[c(-s,-t),]
  matEq[(numVertices-1):numEq,(numEdgesG+1):(numEdgesG + numEdgesGSub)] <- getIncidenceMat(gSub)[c(-s,-t),]
  
  # Inequality constraints
  numInEq <- numVertices-2
  matInEq <- matrix(0, numInEq, numEdgesG + numEdgesGSub)
  matInEq[,1:numEdgesG] <- getIncomingMat(g)[c(-s,-t),]
  matInEq[,(numEdgesG+1):(numEdgesG + numEdgesGSub)] <- getIncomingMat(gSub)[c(-s,-t),]
  
  # Run linear program
  mat <- rbind(matEq,matInEq)
  dir <- c(rep("=",numEq), rep("<=",numInEq))
  rhs <- c(rep(0,numEq), rep(1,numInEq))
  
  res <- lp("max", obj, mat, dir, rhs, all.int=int)
  
  #return(list("res"=res, "mat"=mat, "obj"=obj, "rhs"=rhs))
  return(res)
}  

constructPathSystem <- function(lpRes, g, gSub, s, t){
  
  # Check whether the solution is integer
  if (!all(lpRes$solution - round(lpRes$solution) == 0)){
    print("ERROR, the solution is not integer-valued")
    return(integer(0))
  }
  
  # Define variables
  numEdgesG <- ecount(g)
  numEdgesGSub <- ecount(gSub)
  edgesG <- as_edgelist(g)
  edgesGSub <- as_edgelist(gSub)

  # Paths in subgraph
  activeEdgesSub <- edgesGSub[(lpRes$solution[(numEdgesG+1):(numEdgesG+numEdgesGSub)]==1),]
  if (nrow(activeEdgesSub)==0){
    pathsSub <- numeric(0)
  } else {
    starting_edges <- activeEdgesSub[activeEdgesSub[,1]==s,]
    if (is.vector(starting_edges)){
      starting_edges = matrix(starting_edges,nrow=1)
    }
    pathsSub <- list()
    for (i in 1:nrow(starting_edges)){
      p = matrix(starting_edges[i,],nrow=1) 
      while ( p[nrow(p),2]!=t){
        p = rbind(p,   activeEdgesSub[activeEdgesSub[,1]==p[nrow(p),2],])
      }
      pathsSub[[i]] <- p
    }
  }
  
  # Paths in main graph
  activeEdges <- edgesG[(lpRes$solution[1:numEdgesG]==1),]
  if (nrow(activeEdges)==0){
    pathsMain <- numeric(0)
  } else {
    starting_edges <- activeEdges[activeEdges[,1]==s,]
    if (is.vector(starting_edges)){
      starting_edges = matrix(starting_edges,nrow=1)
    }
    pathsMain <- list()
    for (i in 1:nrow(starting_edges)){
      p = matrix(starting_edges[i,],nrow=1) 
      while ( p[nrow(p),2]!=t){
        p = rbind(p,   activeEdges[activeEdges[,1]==p[nrow(p),2],])
      }
      pathsMain[[i]] <- p
    }
  }
  
  # All paths
  paths <- c(pathsMain, pathsSub)
  return(paths)
}



