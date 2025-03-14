# BAPmultiNMF Package

author: Blake Hansen

Methods to estimate [BAP Multi-NMF](https://doi.org/10.48550/arXiv.2502.01468) via Coordinate Ascent Variational Inference.

## Installation

BAPmultiNMF can be installed using the `devtools` package:

```{r example_install, echo = TRUE, results = TRUE, tidy = TRUE}
library(devtools)
install_github("blhansen/BAPmultiNMF")

library(BAPmultiNMF)
```

## Discovery Only Version
```{r example_1, echo = TRUE, results = TRUE, tidy = TRUE}
# Generate some data
S <- 3
R <- 5
K <- 20
N_s <- c(50, 60, 70)
P <- rdirichlet(R, rep(1/K, K))
E_s <- lapply(1:S, function(s){t(rdirichlet(N_s[s], rgamma(K, shape = 1, rate=0.5)))})
M_s <- lapply(1:S, function(s) matrix(rpois(N_s[s]*K, 500*(P %*% E_s[[s]])), nrow=K, ncol=N_s[s]))

# Run CAVI
approx_1 <- BAPmultiNMF(M_s=M_s, R=5)
```

## Recovery Discovery Version

Simply provide signatures to recover as an argument:

```{r example_2, echo = TRUE, results = TRUE, tidy = TRUE}
approx_2 <- BAPmultiNMF(M_s=M_s, R=5, P_recover = P)
```
## Accessing Signatures

The variational dirichlet concentrations will be returned, so we should normalize the columns of P before examining them:

```{r example_3, echo = TRUE, results = TRUE, tidy = TRUE}
P_estimate <- approx_1$p_conc %*% diag(1/colSums(approx_1$p_conc))
```
