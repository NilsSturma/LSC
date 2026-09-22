library(SEMID)

################################
##### Function Definitions #####
################################

# All subsets of x with at least two elements, largest sets first.
subsetsOfSizeAtLeastTwo <- function(x) {
  if (length(x) < 2) {
    return(list())
  }
  subsets <- list()
  for (k in seq(length(x), 2)) {
    subsets <- c(subsets, combn(x, k, simplify = FALSE))
  }
  return(subsets)
}


# A node is a pure child of l if l is its only parent.
pureChildren <- function(g, l) {
  childrenOfL <- g$children(l)
  isPure <- vapply(childrenOfL, function(c) {
    pa <- g$parents(c)
    length(pa) == 1 && pa == l
  }, logical(1))
  return(childrenOfL[isPure])
}


neighborsOf <- function(g, l) {
  return(union(g$parents(l), g$children(l)))
}


# (emptyset, {l}) trek separates v from c if and only if every trek between
# v and c contains l on its right hand side, i.e., on the side ending in c.
# Note that c is never equal to l, since c is a child of l.
trekSeparatedByLatent <- function(g, v, c, l) {
  trFromV <- g$trFrom(v, avoidRightNodes = l)
  return(!(c %in% trFromV))
}


# Checks the conditions of the theorem for a single latent node l, searching
# over all candidate sets C_l of pure children of l.
checkLatentNode <- function(g, l) {
  allNodes <- c(g$observedNodes(), g$latentNodes())

  for (C in subsetsOfSizeAtLeastTwo(pureChildren(g, l))) {

    # A neighbor n_l of l outside of C that is not a descendant of C.
    # Descendants include C itself, hence n_l is automatically not in C.
    candidateNeighbors <- setdiff(neighborsOf(g, l), g$descendants(C))
    if (length(candidateNeighbors) == 0) {
      next
    }
    n <- candidateNeighbors[1]

    # Every remaining node v has to be trek separated from some c in C
    # by (emptyset, {l}).
    nodesToSeparate <- setdiff(allNodes, c(l, C))
    separators <- integer(0)
    for (v in nodesToSeparate) {
      separatingC <- C[vapply(C, function(c) {
        trekSeparatedByLatent(g, v, c, l)
      }, logical(1))]
      if (length(separatingC) == 0) {
        break
      }
      separators[as.character(v)] <- separatingC[1]
    }

    if (length(separators) == length(nodesToSeparate)) {
      return(list("satisfied" = TRUE,
                  "C" = C,
                  "n" = n,
                  "separators" = separators))
    }
  }

  return(list("satisfied" = FALSE))
}


# Checks whether all direct effects of the latent digraph g are generically
# sign-identifiable by Theorem 5.5
directID <- function(g) {
  latents <- g$latentNodes()

  # One certificate per latent node, named by that node
  certificates <- setNames(lapply(latents, function(l) checkLatentNode(g, l)),
                           latents)
  satisfied <- vapply(certificates, `[[`, logical(1), "satisfied")

  # Latent nodes satisfying the conditions of the theorem
  S <- latents[satisfied]

  result <- list("S" = S,
                 "certificates" = certificates,
                 "id" = all(satisfied))
  return(result)
}


# The examples are only run in an interactive session, so that this file
# can be sourced from a script.
if (interactive()) {

###############
### Example ###
###############
A <- matrix(0, 6, 6)
A[5, 1] <- 1
A[5, 2] <- 1
A[5, 6] <- 1
A[6, 3] <- 1
A[6, 4] <- 1

g <- LatentDigraph(A, observedNodes = 1:4, latentNodes = 5:6)
plot(g)

res <- directID(g)

res$id
res$S

res$certificates[["5"]]$C
res$certificates[["5"]]$n
res$certificates[["5"]]$separators


###########################
### Non-identified case ###
###########################
B <- matrix(0, 6, 6)
B[5, 1] <- 1
B[5, 2] <- 1
B[5, 6] <- 1
B[6, 3] <- 1
B[6, 4] <- 1
B[1, 4] <- 1

gNotID <- LatentDigraph(B, observedNodes = 1:4, latentNodes = 5:6)
plot(gNotID)

resNotID <- directID(gNotID)

resNotID$id
resNotID$S
resNotID$certificates[["5"]]
resNotID$certificates[["6"]]

}
