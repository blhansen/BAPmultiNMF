
#' Coordinate Ascent Variational Inference for Bayesian Probit Multi-Study NMF
#'
#' The function implements the CAVI algorithm described in Hansen et al. (2025+)
#' The code will return the relevant parameters of the mean-field
#' approximation described in the main paper.
#'
#' @param M_s A list of length \eqn{S}, containing counts matrices of dimensions \eqn{K \times N_s}.
#' @param x_s A list of length \eqn{S}, containing covariate matrices of dimensions \eqn{N_s \times D_s}, where \eqn{D_s} is the number of covariates in study s. If no covariates are provided, an intercept will be used for each study.
#' @param R The number of discovered signatures. Default is 5.
#' @param hyperparameters A list containing hyperparameter values. The components are:
#' \describe{
#' \item{alpha_p}{Dirichlet concentration parameter of dimension K for p(P_r). Default is 1.}
#' \item{e_conc}{Dirichlet concentration value for exposures when a=1. Default is 5.}
#' \item{e_conc_null}{Dirichlet concentration value for exposures when a=0. Default is 1.}
#' \item{tau_s_shape}{Gamma shape parameter for p(tau_s). Default is 0.1.}
#' \item{taue_s_rate}{Gamma rate parameter for p(tau_s). Default is 0.1.}
#' \item{beta_prior}{Prior mean of beta_s for the intercept term. Default corresponds to 0.1 baseline probability.}
#' }
#' @param P_recover A matrix of dimensions \eqn{K \times R_{recover}}, containing previously known signatures to be recovered.
#' @param recover_weights Either a single value or list of values of length \eqn{R_{recover}}, controlling the strength of the recovery prior. Higher values correspond to a more informative prior.
#' @param tol Convergence criteria. Default is 1e-5.
#' @param min_iter Minimum number of CAVI iterations. Default is 50.
#' @param max_iter Maximum number of CAVI iterations. Default is 1000.
#' @param verbose Logical value determining if the convergence criteria should be printed at each iteration. Default is TRUE.
#'
#' @return A list  containing the the relevant parameters of the mean-field approximation.  The components of the list are:
#' \item{\code{p_conc}}{A matrix containing the variational dirichlet concentration for the signatues matrix.}
#' \item{\code{e_conc}}{A list of length S, where each entry is a matrix containing the variational dirichlet concentration for the exposures matrix.}
#' \item{\code{a_prob}}{A list of length S, where each entry is a nested list of lengths N_s, R containing the variational inclusion probabilities for each subject and signature.}
#' \item{\code{a_star_mean}}{A list of length S, where each entry is a matrix of dimension \eqn{R \times N_s}, containing the variational posterior means for a_star.}
#' \item{\code{mean_beta}}{A list of length S, where each entry is a list of length R, containing the variational mean for beta_s.}
#' \item{\code{var_beta}}{A list of length S, where each entry is a list of length R, containing the variational variance for beta_s.}
#' \item{\code{W_s}}{A list of length S, where each entry is a list of length N_s, containing the total mutation counts.}
#' \item{\code{priors}}{The list of hyperparameter values used.}
#' @importFrom LaplacesDemon rdirichlet
#' @import NMF
#' @export
#' @references Hansen, B., Grabski, I. N., Parmigiani, G.,  De Vito, R. (2025+).
#' Bayesian Probit Multi-Study Non-negative Matrix Factorization for Mutational Signatures.
#' Submitted manuscript.
#' @examples
#' # Generate Data
#' S <- 3
#' R <- 5
#' K <- 20
#' N_s <- c(50, 60, 70)
#'
#' P <- rdirichlet(R, rep(1/K, K))
#' E_s <- lapply(1:S, function(s){t(rdirichlet(N_s[s], rgamma(K, shape = 1, rate=0.5)))})
#' M_s <- lapply(1:S, function(s) matrix(rpois(N_s[s]*K, 500*(P %*% E_s[[s]])), nrow=K, ncol=N_s[s]))
#'
#' # Run CAVI
#' BAPmultiNMF(M_s=M_s, R=5)
#'
BAPmultiNMF <- function(M_s,
                        x_s = NULL,
                        R = 5,
                        hyperparameters = NULL,
                        P_recover = NULL,
                        recover_weights = 1e8,
                        tol = 1e-5,
                        min_iter = 50,
                        max_iter = 1e4,
                        verbose = TRUE) {
  S <- length(M_s)
  W_s <- lapply(M_s, colSums)

  # Dimensions of Data
  if (length(unique(sapply(M_s, function(x)
    dim(x)[1]))) != 1) {
    warning("K is not the same across all studies!")
  } else {
    K <- unique(sapply(M_s, function(x)
      dim(x)[1]))
  }
  N_s <- sapply(M_s, function(x)
    dim(x)[2])

  if(is.null(x_s)){x_s=lapply(1:length(M_s), function(s) matrix(1, nrow=ncol(M_s[[s]]), ncol=1))}

  # Check for specified priors, otherwise use default values
  default_hyperparameters <- list(
    "alpha_p" = rep(1, K),
    "e_conc" = 5,
    "e_conc_null" = 1,
    "tau_s_shape" = 0.1,
    "tau_s_rate" = 0.1,
    "beta_prior" = qnorm(0.1)
  )
  for (hyperparameter in names(default_hyperparameters)) {
    if (is.null(hyperparameters[[hyperparameter]])) {
      hyperparameters[[hyperparameter]] <-
        default_hyperparameters[[hyperparameter]]
    }
  }

  if (is.null(P_recover)) {
    cosineDist <- function(x, y) {
      (x %*% t(y)) / (sqrt(rowSums(x^2) %*% t(rowSums(y^2))))
    }
    cos_threshold <- 0.6
    data_stacked <- do.call(cbind, M_s)

    R_initial <- R_selected <- R

    nmf_initial = nmf(data_stacked, R_selected)
    P_initial <-
      nmf_initial@fit@W %*% diag(1 / colSums(nmf_initial@fit@W))

    r <- 1
    while (TRUE) {
      cos_mat <- cosineDist(t(P_initial), t(P_initial))
      mat_tri <- cos_mat * lower.tri(cos_mat)

      if (r > ncol(mat_tri)) {
        break
      }
      drop <- (mat_tri[, r] >= cos_threshold)
      if (length(drop > 0)) {
        P_initial <- P_initial[, !drop]
      } else if ((length(drop) == 0) &
                 (max(mat_tri) < cos_threshold)) {
        break
      }

      r <- r + 1
    }

    R_selected <- ncol(P_initial)
    if (verbose == TRUE) {
      print(paste(R_selected, "signatures initialized via NMF,", R - R_selected, "signatures initialized with random values."))
    }

    nmf_initial = nmf(data_stacked, R_selected)
    P_initial <-
      nmf_initial@fit@W %*% diag(1 / colSums(nmf_initial@fit@W))


    P_initial <- cbind(P_initial, t(rdirichlet(ceiling((R - R_selected) /
                                                         2
    ), rep(10, K))), t(rdirichlet(floor((R - R_selected) / 2
    ), rep(0.01, K))))

    Eguess <- array(1, dim = c(ncol(P_initial), ncol(data_stacked)))
    for (i in 1:1000) {
      Eguess <- nmf_update.euclidean.h(
        v = as.matrix(data_stacked),
        w = as.matrix(P_initial),
        h = Eguess
      )
    }

    E_initial <- lapply(1:S, function(s, end) {
      Eguess[, ifelse(s > 1, end[s - 1] + 1, 1):end[s]]
    }, end = cumsum(N_s))
    P_prior <- matrix(rep(hyperparameters$alpha_p, R),
                      nrow = K,
                      ncol = R)
  } else {
    if (dim(P_recover)[1] != K) {
      warning("Dimension of P_recovery does not match M_s!")
    }

    # initialize P,E matrices
    if (length(recover_weights) == 1) {
      recover_weights = rep(recover_weights, ncol(P_recover))
    }
    P_prior <- cbind((P_recover %*% diag(recover_weights))+1e-10,
                     matrix(
                       rep(hyperparameters$alpha_p, R),
                       nrow = K,
                       ncol = R
                     ))
    R <- ncol(P_prior)
    P_initial <- sapply(1:R, function(r)
      rdirichlet(1, P_prior[, r]))


    data_stacked <- do.call(cbind, M_s)
    E_initial <- lapply(1:S, function(s)
      t(rdirichlet(N_s[s], rep(
        hyperparameters$e_conc, R
      ))))
  }

  # Initialize latent sig counts
  z_prob <- lapply(1:S, function(s) {
    array(sapply(1:N_s[s], function(j) {
      sapply(1:K, function(i) {
        probs <-  as.vector(P_initial[i, ], mode = "numeric") * as.vector(E_initial[[s]][, j], mode = "numeric")
        probs <- probs / (sum(probs) + 1e-7)
        return(as.array(probs))
      })
    }, simplify = "array"),
    dim = c(R, K, N_s[s]))
  })
  z_prob <- lapply(1:S, function(s)
    aperm(z_prob[[s]], c(2, 1, 3)))

  # Update P, E based on initial Z
  p_conc <- sapply(1:R, function(r) {
    P_prior[, r] + rowSums(sapply(1:S, function(s)
      apply(M_s[[s]] * z_prob[[s]][, r, ], 1, sum)))
  })
  p_mean <- p_conc %*% diag(1 / colSums(p_conc))

  e_conc <- lapply(1:S, function(s) {
    sapply(1:N_s[s], function(j) {
      e_conc_inc <- hyperparameters$e_conc
      e_conc_null <- hyperparameters$e_conc_null
      (1 / 2) * e_conc_inc +
        (1 / 2) * e_conc_null +
        apply(diag(M_s[[s]][, j]) %*% z_prob[[s]][, , j], 2, sum)
    })
  })

  e_mean <- lapply(1:S, function(s)
    e_conc[[s]] %*% diag(1 / (colSums(e_conc[[s]])) + 1e-7))

  # Initialize probit model
  tau_s_initial <- lapply(1:S, function(s) {
    rgamma(R,
           shape = hyperparameters$tau_s_shape,
           rate = hyperparameters$tau_s_rate)
  })

  beta_s_initial <- lapply(1:S, function(s) {
    matrix(
      rnorm(
        ncol(x_s[[s]]) * R,
        mean = c(rep(hyperparameters$beta_prior, R), rep(0, ncol(x_s[[s]] - 1) *
                                                           (R))),
        sd = sqrt(1 / 10)
      ),
      nrow = ncol(x_s[[s]]),
      ncol = R,
      byrow = TRUE
    )
  })

  inc_probs <- lapply(1:S, function(s) {
    sapply(1:N_s[s], function(j) {
      sapply(1:R, function(r) {
        e_conc_inc <- hyperparameters$e_conc
        e_conc_null <- hyperparameters$e_conc_null

        e_alpha_0 = rep(e_conc_inc * (1 / 2) + e_conc_null * (1 / 2), R)
        e_alpha_0 = sum(e_alpha_0[-r])
        v_alpha_0 = rep((e_conc_inc - e_conc_null)^2 * (1 / 2) * (1 - 1 /
                                                                    2), R)
        v_alpha_0 = sum(v_alpha_0[-r])

        e_alpha_0_inc = e_alpha_0 + e_conc_inc
        e_alpha_0_null = e_alpha_0 + e_conc_null

        expect_beta_rel <- function(e, v)
          lgamma(e) + 1 / 2 * psigamma(e, deriv = 1) * v
        expect_log_exposure = digamma(e_conc[[s]][r, j]) - digamma(sum(e_conc[[s]][, j]))


        rho_1 <- -lgamma(e_conc_inc) + expect_beta_rel(e_alpha_0_inc, v_alpha_0) +
          (e_conc_inc - 1) * expect_log_exposure
        rho_0 <- -lgamma(e_conc_null) + expect_beta_rel(e_alpha_0_null, v_alpha_0) +
          (e_conc_null - 1) * expect_log_exposure

        a <- 0.07056
        b <- 1.5976

        mean_xbeta <- c(matrix(x_s[[s]][j, ], nrow = 1) %*% beta_s_initial[[s]][, r])
        var_xbeta <- c(t(as.matrix(x_s[[s]][j, ])) %*% (tau_s_initial[[s]][r]^-1 * diag(ncol(x_s[[s]]))) %*% as.matrix(x_s[[s]][j, ]))

        probit_0 <- log(pnorm(-mean_xbeta)) + (1/2) * var_xbeta * (pnorm(-mean_xbeta)*(-mean_xbeta/sqrt(2*pi)*exp(-(mean_xbeta)^2/2)) - dnorm(-mean_xbeta)^2)/pnorm(-mean_xbeta)^2

        probit_1 <- log(1-pnorm(-mean_xbeta)) + (1/2) * var_xbeta * ( (1-pnorm(-mean_xbeta))*(-mean_xbeta/sqrt(2*pi)*exp(-(mean_xbeta)^2/2))-dnorm(-mean_xbeta)^2)/(1-pnorm(-mean_xbeta))^2

        log_ratio <- (rho_1 + probit_1) - (rho_0  + probit_0)

        if (is.infinite(log_ratio)) {
          log_ratio <- sign(log_ratio) * 1e3
        } else if (is.nan(log_ratio)){
          log_ratio <- sign(-mean_xbeta) * 1e3
        }
        1 / (1 + exp(-log_ratio))
      })
    })
  })

  a_star <- lapply(1:S, function(s) {
    sapply(1:N_s[s], function(j) {
      sapply(1:R, function(r) {
        mean_xbeta <- c(t(x_s[[s]][j, ]) %*% beta_s_initial[[s]][, r])

        pos_expect <- mean_xbeta + dnorm(-mean_xbeta) / (1 - pnorm(-mean_xbeta))
        neg_expect <- mean_xbeta - dnorm(-mean_xbeta) / (pnorm(-mean_xbeta))

        if (is.infinite(pos_expect) | is.nan(pos_expect)) {
          if ((1 - pnorm(-mean_xbeta)) == 0) {
            pos_expect <- 0
          }
        }

        if (is.infinite(neg_expect) | is.nan(neg_expect)) {
          if ((pnorm(-mean_xbeta)) == 0) {
            neg_expect <- 0
          }
        }

        inc_probs[[s]][r, j] * pos_expect + (1 - inc_probs[[s]][r, j]) *
          neg_expect
      })
    })
  })

  tau_s_shape <- lapply(1:S, function(s) {
    hyperparameters$tau_s_shape + ncol(x_s[[s]]) / 2
  })
  tau_s_rate <- lapply(1:S, function(s) {
    sapply(1:R, function(r) {
      hyperparameters$tau_s_rate + 1 / 2 * crossprod(beta_s_initial[[s]][, r])
    })
  })

  var_beta_s <- lapply(1:S, function(s) {
    lapply(1:R, function(r) {
      solve(tau_s_shape[[s]] / tau_s_rate[[s]][r] * diag(ncol(x_s[[s]])) + crossprod(x_s[[s]]))
    })
  })

  mean_beta_s <- lapply(1:S, function(s) {
    lapply(1:R, function(r) {
      var_beta_s[[s]][[r]] %*% (
        t(x_s[[s]]) %*% a_star[[s]][r, ] + tau_s_shape[[s]] / tau_s_rate[[s]][r] *
          diag(ncol(x_s[[s]])) %*% matrix(
            c(hyperparameters$beta_prior, rep(0, ncol(x_s[[s]]) - 1)),
            nrow = ncol(x_s[[s]]),
            ncol = 1
          )
      )
    })
  })

  # Calculate initial likelihood
  M_s_hat <-
    lapply(1:S, function(s)
      (p_mean %*% e_mean[[s]]) %*% diag(W_s[[s]]))
  loglike <-
    sum(mapply(function(x,y) dpois(x = x, lambda = y, log = TRUE) ,unlist(M_s),unlist(M_s_hat)))

  # Perform iterative CAVI updates
  iter <- 0
  delta <- 1
  trace <- vector(length = max_iter)
  while ((delta > tol) | iter < min_iter) {
    iter <- iter + 1
    if (iter > max_iter)
      break

    z_prob <- lapply(1:S, function(s){
      p_mat <- digamma(p_conc+1e-10)-t(replicate(K, digamma(apply(p_conc, 2, sum)+1e-10)))
      e_mat <- digamma(e_conc[[s]]+1e-10)+t(replicate(R, digamma(apply(e_conc[[s]], 2, sum)+1e-10)))
      out <- exp(replicate(N_s[s], p_mat) + aperm(replicate(K, e_mat), c(3,1,2)))
      out <- sweep(out, c(1,3), apply(out, c(1,3), sum), "/")
      return(out)
    })

    p_conc <- sapply(1:R, function(r) {
      P_prior[, r] + apply(
        sapply(1:S, function(s)apply(M_s[[s]] * z_prob[[s]][, r, ], 1, sum)), 1, sum
      )
    })

    p_mean <- sweep(p_conc, 2, apply(p_conc, 2, sum), "/")

    e_conc <- lapply(1:S, function(s) {
      sapply(1:N_s[s], function(j) {
        e_conc_inc <- hyperparameters$e_conc
        e_conc_null <- hyperparameters$e_conc_null
        e_conc_inc * inc_probs[[s]][, j] + e_conc_null * (1 - inc_probs[[s]][, j]) + apply(diag(M_s[[s]][, j]) %*% z_prob[[s]][, , j], 2, sum)
      })
    })
    e_mean <- lapply(1:S, function(s)
      e_conc[[s]] %*% diag(1 / colSums(e_conc[[s]])))

    inc_probs <- lapply(1:S, function(s) {
      sapply(1:N_s[s], function(j) {
        sapply(1:R, function(r) {
          e_conc_inc <- hyperparameters$e_conc
          e_conc_null <- hyperparameters$e_conc_null
          e_alpha_0 = e_conc_inc * inc_probs[[s]][, j] + e_conc_null * (1 -
                                                                          inc_probs[[s]][, j])
          e_alpha_0 = sum(e_alpha_0[-r])
          v_alpha_0 = (e_conc_inc - e_conc_null)^2 * (inc_probs[[s]][, j]) *
            (1 - inc_probs[[s]][, j])
          v_alpha_0 = sum(v_alpha_0[-r])

          e_alpha_0_inc = e_alpha_0 + e_conc_inc
          e_alpha_0_null = e_alpha_0 + e_conc_null

          expect_beta_rel <- function(e, v)
            lgamma(e) + 1 / 2 * psigamma(e, deriv = 1) * v
          expect_log_exposure = digamma(e_conc[[s]][r, j]) - digamma(sum(e_conc[[s]][, j]))


          rho_1 <- -lgamma(e_conc_inc) + expect_beta_rel(e_alpha_0_inc, v_alpha_0) +
            (e_conc_inc - 1) * expect_log_exposure
          rho_0 <- -lgamma(e_conc_null) + expect_beta_rel(e_alpha_0_null, v_alpha_0) +
            (e_conc_null - 1) * expect_log_exposure

          a <- 0.07056
          b <- 1.5976

          mean_xbeta <- c(t(x_s[[s]][j, ]) %*% mean_beta_s[[s]][[r]])
          var_xbeta <- c(t(x_s[[s]][j, ]) %*% (var_beta_s[[s]][[r]]) %*% x_s[[s]][j, ])

          probit_0 <- log(pnorm(-mean_xbeta)) + (1/2) * var_xbeta * (pnorm(-mean_xbeta)*(-mean_xbeta/sqrt(2*pi)*exp(-(mean_xbeta)^2/2)) - dnorm(-mean_xbeta)^2)/pnorm(-mean_xbeta)^2

          probit_1 <- log(1-pnorm(-mean_xbeta)) + (1/2) * var_xbeta * ( (1-pnorm(-mean_xbeta))*(-mean_xbeta/sqrt(2*pi)*exp(-(mean_xbeta)^2/2))-dnorm(-mean_xbeta)^2)/(1-pnorm(-mean_xbeta))^2


          log_ratio <- (rho_1 + probit_1) - (rho_0  + probit_0)

          if (is.infinite(log_ratio)) {
            log_ratio <- sign(log_ratio) * 1e3
          } else if (is.nan(log_ratio)){
            log_ratio <- sign(-mean_xbeta) * 1e3
          }


          1 / (1 + exp(-log_ratio))
        })
      })
    })

    a_star <- lapply(1:S, function(s) {
      sapply(1:N_s[s], function(j) {
        sapply(1:R, function(r) {
          mean_xbeta <- c(x_s[[s]][j, ] %*% mean_beta_s[[s]][[r]])

          pos_expect <- mean_xbeta + dnorm(-mean_xbeta) / (1 - pnorm(-mean_xbeta) + 1e-5)
          neg_expect <- mean_xbeta - dnorm(-mean_xbeta) / (pnorm(-mean_xbeta) + 1e-5)

          if (is.infinite(pos_expect) | is.nan(pos_expect)) {
            if ((1 - pnorm(-mean_xbeta)) == 0) {
              pos_expect <- 0
            }
          }

          if (is.infinite(neg_expect) | is.nan(neg_expect)) {
            if ((pnorm(-mean_xbeta)) == 0) {
              neg_expect <- 0
            }
          }

          inc_probs[[s]][r, j] * pos_expect + (1 - inc_probs[[s]][r, j]) *
            neg_expect
        })
      })
    })

    tau_s_shape <- lapply(1:S, function(s) {
      hyperparameters$tau_s_shape + ncol(x_s[[s]]) / 2
    })
    tau_s_rate <- lapply(1:S, function(s) {
      sapply(1:R, function(r) {
        hyperparameters$tau_s_rate + 1 / 2 * (crossprod(beta_s_initial[[s]][, r] - matrix(c(hyperparameters$beta_prior, rep(0, ncol(x_s[[s]]) - 1)),
                                                                                          nrow = ncol(x_s[[s]]),
                                                                                          ncol = 1)
        ) + sum(diag(var_beta_s[[s]][[r]]))
        )
      })
    })

    var_beta_s <- lapply(1:S, function(s) {
      lapply(1:R, function(r) {
        solve(tau_s_shape[[s]] / tau_s_rate[[s]][r] * diag(ncol(x_s[[s]])) + crossprod(x_s[[s]]))
      })
    })

    mean_beta_s <- lapply(1:S, function(s) {
      lapply(1:R, function(r) {
        var_beta_s[[s]][[r]] %*% (
          t(x_s[[s]]) %*% a_star[[s]][r, ] +
            tau_s_shape[[s]] / tau_s_rate[[s]][r] *
            diag(ncol(x_s[[s]])) %*% matrix(c(hyperparameters$beta_prior, rep(0, ncol(x_s[[s]]) - 1)),
                                            nrow = ncol(x_s[[s]]),
                                            ncol = 1
            )
        )
      })
    })

    M_s_hat <-
      lapply(1:S, function(s)
        (p_mean %*% e_mean[[s]]) %*% diag(W_s[[s]]))
    loglike_new <-
      sum(mapply(function(x,y) dpois(x = x, lambda = y, log = TRUE) ,unlist(M_s),unlist(M_s_hat)))

    delta <- abs((loglike - loglike_new) / loglike)
    loglike <- loglike_new
    trace[iter] <- loglike

    if (verbose) {
      print(paste0("Iteration ", iter, " Delta: ", signif(delta, 3)))
    }
  }

  # Return relevant parameters
  ests <- list(
    "p_conc" = p_conc,
    "e_conc" = e_conc,
    "inc_prob" = inc_probs,
    "a_star" = a_star,
    "mean_beta" = mean_beta_s,
    "var_beta" = var_beta_s,
    "W_s" = W_s,
    "priors" = hyperparameters
  )
  return(ests)
}
