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
  
  # Explicit error nodes
  error <- setNames(max(nodes) + seq_along(nodes), nodes)
  
  transformedNodes <- c(observed, latents, unname(error))
  
  B <- matrix(
    0,
    length(transformedNodes),
    length(transformedNodes),
    dimnames = list(transformedNodes, transformedNodes)
  )
  
  # Original graph and independent error variables
  B[as.character(nodes), as.character(nodes)] <- A
  B[cbind(as.character(error), as.character(nodes))] <- 1
  
  # L2O transformation.
  B[
    as.character(allParents),
    as.character(outcome)
  ] <- 0
  
  B[
    as.character(outcome),
    as.character(outcomeIndicator)
  ] <- 0
  
  B[
    as.character(error[as.character(outcome)]),
    as.character(outcomeIndicator)
  ] <- 1
  
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
  # outcome in the L2O graph.
  candidates <- setdiff(
    observed,
    trFrom(
      transformed,
      outcomeIndicator,
      includeLatents = FALSE
    )
  )

  # Find a MIIV
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


# All valid scaling indicators of a latent node, i.e. all pure observed children
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


# Find one valid assignment of scaling indicators
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


# Applies the MIIV criterion via a L2O transformation.
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




# The examples are only run in an interactive session, so that this file
# can be sourced from a script.
if (interactive()) {

###############
### Example ###
###############
A <- matrix(0, 10, 10)

A[8, 1] <- 1
A[8, 2] <- 1

A[9, 3] <- 1
A[9, 4] <- 1

A[10, 5] <- 1
A[10, 6] <- 1

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

result$scalingIndicators
identifyL2OIVforIndicators(graph, c(`8` = 2, `9` = 4, `10` = 6))$allIdentified

result$equations[["10"]]
plot(result$equations[["10"]]$transformedGraph)


############################
### More complex example ###
############################
A2 <- matrix(0, 12, 12)

A2[10, 1] <- 1
A2[11, 2] <- 1
A2[12, 3] <- 1

A2[10, 4] <- 1
A2[11, 5] <- 1
A2[12, 6] <- 1
A2[9, 4:6] <- 1

A2[10, 11] <- 1
A2[11, 12] <- 1

A2[7, 11] <- 1
A2[8, 12] <- 1

A2[7, 8] <- 1
A2[8, 9] <- 1

graph2 <- LatentDigraph(
  A2,
  observedNodes = 1:9,
  latentNodes = 10:12
)
plot(graph2)

result2 <- identifyL2OIV(graph2)

result2$scalingIndicators

result2$allIdentified
result2$equations[["11"]]
plot(result2$equations[["11"]]$transformedGraph)
result2$equations[["12"]]
plot(result2$equations[["12"]]$transformedGraph)


####################################
### Sparse non-identified example ###
####################################
A3 <- matrix(0, 8, 8)

A3[6, 1] <- 1
A3[7, 2] <- 1
A3[8, 3] <- 1

A3[6, 7] <- 1
A3[7, 8] <- 1

A3[4, 7] <- 1
A3[5, 8] <- 1
A3[4, 5] <- 1

graph3 <- LatentDigraph(
  A3,
  observedNodes = 1:5,
  latentNodes = 6:8
)
plot(graph3)

result3 <- identifyL2OIV(graph3)

result3$scalingIndicators

result3$allIdentified
result3$equations[["7"]]
plot(result3$equations[["7"]]$transformedGraph)
result3$equations[["8"]]
plot(result3$equations[["8"]]$transformedGraph)
}
