source("R/semi-direct-id.R", chdir = TRUE)
source("R/canonicalization.R", chdir = TRUE)
library(SEMID)

################################
##### Function Definitions #####
################################

# Observed covariance matrix
observedCov <- function(L, omega, observedNodes) {
  IL <- solve(diag(nrow(L)) - L)
  return((t(IL) %*% diag(omega) %*% IL)[observedNodes, observedNodes])
}


# Semi-direct effect matrix
semiDirectEffects <- function(L, observedNodes, latentNodes) {
  return(L[observedNodes, observedNodes] + L[observedNodes, latentNodes] %*%
    solve(diag(length(latentNodes)) - L[latentNodes, latentNodes]) %*%
    L[latentNodes, observedNodes])
}


# Random edge weights on the support of supp
randomWeights <- function(supp, min = 0.5, max = 1.5) {
  L <- supp
  L[supp == 1] <- runif(sum(supp), min, max)
  return(L)
}


fiberRanks <- function(supp, observedNodes, latentNodes, tol = 1e-6) {
  n <- nrow(supp)
  E <- which(supp == 1, arr.ind = TRUE)
  m <- nrow(E)
  build <- function(th) {
    L <- matrix(0, n, n); L[E] <- th[1:m]
    return(L)
  }
  sigma <- function(th) {
    S <- observedCov(build(th), th[m + (1:n)], observedNodes)
    return(S[upper.tri(S, diag = TRUE)])
  }
  effects <- function(th) {
    return(as.vector(semiDirectEffects(build(th), observedNodes, latentNodes)))
  }
  jacobian <- function(f, th, eps = 1e-6) {
    return(sapply(seq_along(th), function(i) {
      tp <- th; tm <- th; tp[i] <- tp[i] + eps; tm[i] <- tm[i] - eps
      return((f(tp) - f(tm)) / (2 * eps))
    }))
  }
  rank <- function(M) {
    d <- svd(M)$d
    return(sum(d > tol * max(d, 1)))
  }
  th <- c(runif(m, 0.5, 1.5), runif(n, 0.5, 1.5))
  J <- jacobian(sigma, th)
  return(c(rank(J), rank(rbind(J, jacobian(effects, th)))))
}


###############################################################
### (i) The LSC succeeds where the LF-HTC cannot succeed    ###
###############################################################
A <- matrix(0, 7, 7)
A[4, 1] <- 1                      # the semi-direct effect of interest
A[5, 4] <- 1; A[5, 7] <- 1
A[7, 3] <- 1; A[7, 6] <- 1
A[6, 1] <- 1; A[6, 2] <- 1

g <- LatentDigraph(A, observedNodes = 1:4, latentNodes = 5:7)
plot(g)

# The LSC identifies the semi-direct effect matrix.
res <- LSCID(g)
res$id


res$Ys[[1]]; res$Zs[[1]]; res$H2s[[1]]

gCan <- canonicalization(g)
plot(gCan)
lfhtcID(gCan)
LSCID(gCan)$id

# Formula to recover the effect 4 -> 1, implied by proof of LSC
set.seed(1)
lambda <- 0.7
L <- randomWeights(A)
L[4, 1] <- lambda
S <- observedCov(L, runif(7, 0.5, 1.5), 1:4)

solve(rbind(c(S[4, 4], S[2, 4]),
            c(S[4, 3], S[2, 3])),
      c(S[1, 4], S[1, 3]))[1]     # equals lambda


########################################################################################
### (ii) The LSC fails, but the semi-direct effect matrix is rationally identifiable ###
########################################################################################
B <- matrix(0, 5, 5)
B[5, 1] <- 1; B[5, 2] <- 1; B[5, 3] <- 1; B[5, 4] <- 1
B[2, 4] <- 1; B[3, 4] <- 1

gFactor <- LatentDigraph(B, observedNodes = 1:4, latentNodes = 5)
plot(gFactor)

# Neither criterion identifies anything.
LSCID(gFactor)$id
lfhtcID(gFactor)

# Both semi-direct effects are nevertheless rational functions of Sigma. 
set.seed(2)
lambda24 <- 0.9; lambda34 <- -0.4
M <- randomWeights(B)
M[2, 4] <- lambda24; M[3, 4] <- lambda34
S <- observedCov(M, runif(5, 0.5, 1.5), 1:4)

v <- c(S[1, 2] * S[1, 3], S[1, 2] * S[2, 3], S[1, 3] * S[2, 3])
solve(cbind(S[2, 1:3], S[3, 1:3], v), S[4, 1:3])[1:2]   # equals the two effects


###########################################################
### (iii) A semi-direct effect that is not identifiable ###
###########################################################
C <- matrix(0, 3, 3)
C[1, 2] <- 1
C[3, 1] <- 1; C[3, 2] <- 1

gBow <- LatentDigraph(C, observedNodes = 1:2, latentNodes = 3)
plot(gBow)

LSCID(gBow)$id

# Here it is known that the effect 1 -> is not rationally identifiable.
bowL <- function(lambda, a, b) {
  L <- matrix(0, 3, 3)
  L[1, 2] <- lambda; L[3, 1] <- a; L[3, 2] <- b
  return(L)
}

# Two different parameter configurations that give the same covariance matrix
observedCov(bowL(0.5, 1, 1.0), c(1, 1.00, 1), 1:2)
observedCov(bowL(0.8, 1, 0.4), c(1, 1.42, 1), 1:2)


########################################
### (iv) Figure 2 (a) from the paper ###
########################################
# Observed nodes v1, ..., v5 and latent nodes h1 = 6, h2 = 7.
D <- matrix(0, 7, 7)
D[1, 6] <- 1                                          # v1 -> h1
D[6, 7] <- 1                                          # h1 -> h2
D[7, 2] <- 1; D[7, 3] <- 1; D[7, 4] <- 1; D[7, 5] <- 1
D[3, 4] <- 1

gChain <- LatentDigraph(D, observedNodes = 1:5, latentNodes = 6:7)
plot(gChain)

LSCID(gChain)$id

set.seed(3)
fiberRanks(D, 1:5, 6:7) # equal ranks, semi-direct effects are finite-to-one


###############################################################
### (v) Two added edges out of h1 destroy identifiability   ###
###############################################################
E <- D
E[6, 3] <- 1; E[6, 4] <- 1

gTwoFactors <- LatentDigraph(E, observedNodes = 1:5, latentNodes = 6:7)
plot(gTwoFactors)

LSCID(gTwoFactors)$id

# The second rank is larger, so the semi-direct effects vary along the fibers
# of Sigma and are not a function of Sigma, let alone a rational one.
set.seed(4)
fiberRanks(E, 1:5, 6:7)


###############################################################
### (vi) One added edge breaks the LSC but not rational     ###
###      identifiability                                    ###
###############################################################
G <- D
G[2, 4] <- 1

gExtraEdge <- LatentDigraph(G, observedNodes = 1:5, latentNodes = 6:7)
plot(gExtraEdge)

LSCID(gExtraEdge)$id

set.seed(5)
fiberRanks(G, 1:5, 6:7)                               # equal ranks

# The effects into v4 are rational in Sigma. Conditioning on the exogenous node
# v1 removes the latent chain and leaves a one factor model on v2, ..., v5 in
# which v2, v3, v5 are pure indicators, so the rank one argument of (ii)
# applies to the conditional covariance matrix C.
set.seed(6)
L <- randomWeights(G)
S <- observedCov(L, runif(7, 0.5, 1.5), 1:5)

semiDirectEffects(L, 1:5, 6:7)                        # the effects to be recovered

C <- S[2:5, 2:5] - outer(S[1, 2:5], S[1, 2:5]) / S[1, 1]
v <- c(C[1, 2] * C[1, 4], C[1, 2] * C[2, 4], C[1, 4] * C[2, 4])
lambda <- solve(cbind(C[c(1, 2, 4), 1], C[c(1, 2, 4), 2], v), C[c(1, 2, 4), 3])[1:2]

# The effects of v1 follow from its total effects, since v1 is exogenous.
tau <- S[1, 2:5] / S[1, 1]

Lambda <- matrix(0, 5, 5)
Lambda[2, 4] <- lambda[1]
Lambda[3, 4] <- lambda[2]
Lambda[1, 2:5] <- c(tau[1], tau[2], tau[3] - lambda[1] * tau[1] - lambda[2] * tau[2], tau[4])

Lambda
all.equal(Lambda, semiDirectEffects(L, 1:5, 6:7))
