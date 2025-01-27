# BAPmultiNMF

author: Blake Hansen

Methods to estimate BAP Multi-NMF via Coordinate Ascent Variational Inference.

## Discovery Only Version
# Generate some data
S <- 3
R <- 5
K <- 20
N_s <- c(50, 60, 70)
P <- rdirichlet(R, rep(1/K, K))
E_s <- lapply(1:S, function(s){t(rdirichlet(N_s[s], rgamma(K, shape = 1, rate=0.5)))})
M_s <- lapply(1:S, function(s) matrix(rpois(N_s[s]*K, 500*(P %*% E_s[[s]])), nrow=K, ncol=N_s[s]))

# Run CAVI
BAPmultiNMF(M_s=M_s, R=5)

## Recovery Discovery Version

Simply provide signatures to recover as an argument:

BAPmultiNMF(M_s=M_s, R=5, P_recover = P)

