# Latent Subgraph Criterion
This repository contains an implementation of the latent subgraph criterion  in `R`. It builds on the `SEMID` library available on `CRAN`. The function `LSCID(g)` checks the latent subgraph criterion by solving integer linear programs. Here, `g` is a `LatentDigraph` object representing a graph, for more details see the `SEMID`documentation.

```         
> # Latent digraphs are specified by their directed adjacency matrix L
> library(SEMID)
> L = matrix(c(0, 0, 0, 0, 0, 1, 0, 0,
+              0, 0, 0, 0, 0, 0, 0, 0,
+              0, 0, 0, 0, 0, 0, 0, 0,
+              0, 0, 0, 0, 0, 0, 0, 0,
+              0, 0, 0, 0, 0, 0, 0, 0,
+              0, 0, 1, 1, 0, 0, 0, 0,
+              1, 1, 0, 0, 0, 0, 0, 1,
+              0, 1, 1, 1, 1, 0, 0, 0), 8, 8, byrow=TRUE)
> observedNodes = seq(5)
> latentNodes = c(6,7,8)
>
> # Create the latent digraph object corresponding to L
> g = LatentDigraph(L, observedNodes, latentNodes)
>
> # Plot latent digraph
> plot(g)
>
> # Check the latent subgraph criterion
> res <- LSCID(g)
> res$id
```


To reproduce the experimental results in the paper, run:
```
Rscript experiments/random-exps.R
Rscript experiments/compute-stats.R
Rscript experiments/latent-effect-exps.R
```
