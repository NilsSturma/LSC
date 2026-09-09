library(SEMID)


L2OIVequation <- function(graph, outcome, indicators) {
  indicators <- setNames(as.integer(unlist(indicators)), names(indicators))
  
  observed <- observedNodes(graph)
  latents  <- latentNodes(graph)
  nodes    <- c(observed, latents)
  
  allParents      <- parents(graph, outcome)
  latentParents   <- intersect(allParents, latents)
  observedParents <- intersect(allParents, observed)
  
  outcomeIndicator <- as.integer(
    indicators[[as.character(outcome)]]
  )
  parentIndicators <- as.integer(unlist(
    indicators[as.character(latentParents)],
    use.names = FALSE
  ))
  
  # Complete L2O regressor vector
  regressors <- c(parentIndicators, observedParents)
  
  A <- L(graph)
  dimnames(A) <- list(nodes, nodes)
  
  # Explicit error nodes encode the transformed composite error.
  # No error variance is assigned or fixed.
  error <- setNames(max(nodes) + seq_along(nodes), nodes)
  
  transformedNodes <- c(observed, latents, unname(error))
  
  B <- matrix(
    0,
    length(transformedNodes),
    length(transformedNodes),
    dimnames = list(transformedNodes, transformedNodes)
  )
  
  # Original graph and independent error sources
  B[as.character(nodes), as.character(nodes)] <- A
  B[cbind(as.character(error), as.character(nodes))] <- 1
  
  # Whole-equation L2O transformation.
  #
  # All transformed regressor arrows are omitted because the
  # instrumental-set criterion is applied to the edge-deleted graph.
  B[
    as.character(allParents),
    as.character(outcome)
  ] <- 0
  
  B[
    as.character(outcome),
    as.character(outcomeIndicator)
  ] <- 0
  
  # Outcome disturbance enters the transformed outcome.
  B[
    as.character(error[as.character(outcome)]),
    as.character(outcomeIndicator)
  ] <- 1
  
  # Measurement errors of latent-parent indicators enter the
  # transformed composite error.
  B[
    as.character(error[as.character(parentIndicators)]),
    as.character(outcomeIndicator)
  ] <- 1
  
  transformed <- LatentDigraph(
    B,
    observedNodes = observed,
    latentNodes = c(latents, unname(error))
  )
  
  # Valid candidate instruments have no trek to the transformed
  # outcome in the edge-deleted L2O graph.
  candidates <- setdiff(
    observed,
    trFrom(
      transformed,
      outcomeIndicator,
      includeLatents = FALSE
    )
  )
  

  # Find a full-rank instrumental set for the complete equation.
  fit <- if (length(candidates) >= length(regressors)) {
    getTrekSystem(
      transformed,
      fromNodes = candidates,
      toNodes = regressors
    )
  } else {
    list(
      systemExists = FALSE,
      activeFrom = integer()
    )
  }
  
  list(
    outcome = outcome,
    latentParents = latentParents,
    observedParents = observedParents,
    latentEffects = paste(
      latentParents,
      outcome,
      sep = " -> "
    ),
    regressors = regressors,
    candidateInstruments = candidates,
    instruments = fit$activeFrom,
    identified = fit$systemExists,
    transformedGraph = transformed
  )
}


# All valid scaling indicators of a latent node, i.e. all pure observed
# children: observed children whose only parent is that latent node.
validScalingIndicators <- function(graph, h) {
  observedChildren <- children(
    graph,
    h,
    includeObserved = TRUE,
    includeLatents = FALSE
  )

  isPure <- vapply(observedChildren, function(s) {
    pa <- parents(graph, s)
    length(pa) == 1 && pa == h
  }, logical(1))

  return(observedChildren[isPure])
}


# For every latent node, checks that the indicators of its latent parents are
# not also observed parents of it. Returns one logical value per latent node.
disjointParentIndicators <- function(graph, indicators) {
  observed <- observedNodes(graph)
  latents  <- latentNodes(graph)

  return(vapply(latents, function(y) {
    pa <- parents(graph, y)
    latentParents <- intersect(pa, latents)
    observedParents <- intersect(pa, observed)
    parentIndicators <- indicators[as.character(latentParents)]

    length(intersect(parentIndicators, observedParents)) == 0
  }, logical(1)))
}


# Applies the L2O instrumental set criterion for one given choice of scaling
# indicators, one per latent node.
identifyL2OIVforIndicators <- function(graph, indicators) {
  latents  <- latentNodes(graph)

  if (!setequal(names(indicators), as.character(latents)) ||
      any(lengths(indicators) != 1)) {
    stop("Assign exactly one scaling indicator to every latent node.")
  }
  
  indicators <- setNames(as.integer(unlist(indicators)), names(indicators))
  
  if (anyDuplicated(indicators)) {
    stop("Scaling indicators must be distinct.")
  }
  
  # Check that every scaling indicator is a pure observed child.
  valid <- vapply(latents, function(h) {
    indicators[as.character(h)] %in% validScalingIndicators(graph, h)
  }, logical(1))
  
  if (!all(valid)) {
    stop(
      "Invalid scaling indicator for latent node(s): ",
      paste(latents[!valid], collapse = ", ")
    )
  }

  # An indicator of a latent parent cannot also be an observed parent of the
  # same latent node.
  disjointParents <- disjointParentIndicators(graph, indicators)

  if (!all(disjointParents)) {
    stop(
      "Indicators of latent parents must not also be observed parents ",
      "of latent node(s): ",
      paste(latents[!disjointParents], collapse = ", ")
    )
  }
  
  # Check precisely those latent outcomes having latent parents.
  outcomes <- latents[vapply(latents, function(y) {
    any(parents(graph, y) %in% latents)
  }, logical(1))]
  
  equations <- setNames(
    lapply(
      outcomes,
      function(y) {
        L2OIVequation(
          graph,
          outcome = y,
          indicators = indicators
        )
      }
    ),
    outcomes
  )
  
  list(
    allIdentified = all(vapply(
      equations,
      `[[`,
      logical(1),
      "identified"
    )),
    scalingIndicators = indicators,
    equations = equations
  )
}


# One valid choice of scaling indicators, i.e. one pure observed child per
# latent node such that no indicator of a latent parent is also an observed
# parent of the same latent node. The criterion does not depend on which valid
# choice is used, hence the search stops at the first one. A pure observed child
# has exactly one parent, so the candidates of different latent nodes are
# automatically distinct. Returns NULL if there is no valid choice.
firstValidIndicators <- function(graph) {
  latents <- latentNodes(graph)

  if (length(latents) == 0) {
    return(setNames(integer(0), character(0)))
  }

  candidates <- setNames(
    lapply(latents, function(h) validScalingIndicators(graph, h)),
    latents
  )

  if (any(lengths(candidates) == 0)) {
    return(NULL)
  }

  grid <- expand.grid(candidates, KEEP.OUT.ATTRS = FALSE)

  for (i in seq_len(nrow(grid))) {
    indicators <- setNames(as.integer(grid[i, ]), latents)

    if (all(disjointParentIndicators(graph, indicators))) {
      return(indicators)
    }
  }

  return(NULL)
}


# Applies the L2O instrumental set criterion. The scaling indicators are not
# passed but chosen automatically, since the criterion does not depend on which
# valid choice of one pure observed child per latent node is used.
identifyL2OIV <- function(graph) {
  indicators <- firstValidIndicators(graph)

  if (is.null(indicators)) {
    return(list(
      allIdentified = FALSE,
      scalingIndicators = integer(0),
      equations = list()
    ))
  }

  return(identifyL2OIVforIndicators(graph, indicators))
}




###############
### Example ###
###############
A <- matrix(0, 10, 10)

# Pure indicators of latent 8
A[8, 1] <- 1
A[8, 2] <- 1

# Pure indicators of latent 9
A[9, 3] <- 1
A[9, 4] <- 1

# Pure indicators of latent 10
A[10, 5] <- 1
A[10, 6] <- 1

# Observed and latent parents of latent 10
A[7, 10] <- 1
A[8, 10] <- 1
A[9, 10] <- 1

graph <- LatentDigraph(
  A,
  observedNodes = 1:7,
  latentNodes = 8:10
)
plot(graph)

result <- identifyL2OIV(graph)

result$allIdentified

# The first valid choice of scaling indicators, one pure observed child per
# latent node. Every latent node has two of them here, and any other valid
# choice gives the same result.
result$scalingIndicators
identifyL2OIVforIndicators(graph, c(`8` = 2, `9` = 4, `10` = 6))$allIdentified

result$equations[["10"]]
plot(result$equations[["10"]]$transformedGraph)


############################
### More complex example ###
############################
# Observed nodes: 1, ..., 9; latent nodes: 10, 11, 12.
A2 <- matrix(0, 12, 12)

# Each latent has exactly one pure child, used as its scaling indicator.
A2[10, 1] <- 1
A2[11, 2] <- 1
A2[12, 3] <- 1

# Additional, non-pure indicators. Observed node 9 is also their parent.
A2[10, 4] <- 1
A2[11, 5] <- 1
A2[12, 6] <- 1
A2[9, 4:6] <- 1

# Latent variables 11 and 12 both have a latent parent.
A2[10, 11] <- 1
A2[11, 12] <- 1

# Directed edges from observed to latent variables.
A2[7, 11] <- 1
A2[8, 12] <- 1

# Directed edges between observed variables.
A2[7, 8] <- 1
A2[8, 9] <- 1

graph2 <- LatentDigraph(
  A2,
  observedNodes = 1:9,
  latentNodes = 10:12
)
plot(graph2)

# Each latent node has exactly one pure observed child, hence there is only
# one valid choice of scaling indicators.
result2 <- identifyL2OIV(graph2)

result2$scalingIndicators

# There is one equation for each latent outcome with a latent parent.
result2$allIdentified
result2$equations[["11"]]
plot(result2$equations[["11"]]$transformedGraph)
result2$equations[["12"]]
plot(result2$equations[["12"]]$transformedGraph)


####################################
### Sparse non-identified example ###
####################################
# Observed nodes: 1, ..., 5; latent nodes: 6, 7, 8.
A3 <- matrix(0, 8, 8)

# Each latent has exactly one pure child.
A3[6, 1] <- 1
A3[7, 2] <- 1
A3[8, 3] <- 1

# Latent variables 7 and 8 both have a latent parent.
A3[6, 7] <- 1
A3[7, 8] <- 1

# One observed parent per latent outcome and one observed-to-observed edge.
A3[4, 7] <- 1
A3[5, 8] <- 1
A3[4, 5] <- 1

graph3 <- LatentDigraph(
  A3,
  observedNodes = 1:5,
  latentNodes = 6:8
)
plot(graph3)

# Again only one valid choice of scaling indicators, and the criterion fails.
result3 <- identifyL2OIV(graph3)

result3$scalingIndicators

# The equation for latent outcome 7 has no full-rank instrumental set.
result3$allIdentified
result3$equations[["7"]]
plot(result3$equations[["7"]]$transformedGraph)
result3$equations[["8"]]
plot(result3$equations[["8"]]$transformedGraph)
