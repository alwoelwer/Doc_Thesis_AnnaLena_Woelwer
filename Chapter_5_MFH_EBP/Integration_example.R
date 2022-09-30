library(fOptions)
library(mvtnorm)


source(paste0(getwd(), "/Chapter_5_MFH_EBP/MFH_gen_dat_m3.R"))
source(paste0(getwd(), "/Chapter_6_MMFH/MMFH_fitting.R"))

# load Gauss-Hermite Weights and Abscissas
GH_weights <- readRDS(file = paste0(getwd(), "/Chapter_5_MFH_EBP/GH_weights_norm_until_n200.RData"))
names(GH_weights)



## Generate example data

# Generate the fixed quantities of a MMFH model including randomly generated auxiliary information (documentation of required inputs in `MFH_gen_dat_m3.R`).

# Set input quantities

m <- 3 # total number of variables of interest, number of dependent variables
D <- 100 # total number of domains
v_ref <- c(2, 3, 4) # variances of random effects of the 3 variables
cor_ref <- c(.2, .3, .4) # correlations of random effects
v_rer <- c(2.5, 3.5, 4.5) # variances of sampling errors of the 3 variables
cor_rer <- c(.15, .25, .35) # correlations of sampling errors
beta <- list(
  c(1.5, 2.5), # list, for each variable: vector of fixed effects
  c(2.3, 3.3, 4.3),
  c(4.1, 3.1, 2.1, 2.2)
)
range_aux <- c(10, 100) # range of the uniform distribution from which auxiliary information is sampled from

# Generate data

d_fix <- f_prep_MFH_m3(seed = 56)
names(d_fix)

# Use `f_gen_MFH_m3` to generate the model information which typically varies between Monte Carlo iterations. That is, the generation of the dependent variables. The function allows to set certain values as missing, input `perc_mis` determines the number of domains for which the survey information of the three dependent variables is missing. Note that in this code the missing dependent variables are non-overlapping. That is, there is maximum one missing dependent variable per domain. This, however, can easily be changed in the code.

d_var <- f_gen_MFH_m3(
  x = d_fix$x,
  beta = d_fix$beta,
  V_ud_est = d_fix$V_ud_est,
  V_ed = d_fix$V_ed,
  seed = 67,
  verbose = TRUE
)
str(d_var)
y <- d_var$y_obs
V_ed <- d_fix$V_ed

### Fit a MFH model
res_MFH <- f_MMFH(
  y = y,
  x = d_fix$x,
  V_ed = V_ed
)

# From the model results: calculate certain quantities

# Calculate X beta
x_bdiag <- do.call(
  "rbind",
  lapply(1:D, function(d) {
    Matrix::Matrix(do.call(
      Matrix::bdiag,
      lapply(d_fix$x, function(x_d) {
        matrix(x_d[d, ], nrow = 1)
      })
    ), sparse = TRUE)
  })
)
beta_est <- res_MFH$est$fit$estcoef[, 1]
x_beta_est <- matrix(x_bdiag %*% beta_est, ncol = m, byrow = TRUE)
residuals <- y - x_beta_est
V_ud_est <- f_Cov(theta = res_MFH$est$fit$refvar)


# ____________________________________________________________________________
# Integration: Rate (y1/(y1+y2)) ----
# ____________________________________________________________________________

 ### Preparation Gauss-Hermite

  # Choose the number of function evaluations per dimension
  # Product rule: Choosing gauss_n = 5 results in a total of 25 function
  # evaluations
  gauss_n <- 10

  # Get the corresponding points and weights
  w_tmp <- GH_weights[[gauss_n]]

  # product rule: Set up the grid of points and weights
  idx <- as.matrix(expand.grid(rep(list(1:gauss_n), m)))
  pts <- matrix(w_tmp[idx, 1], nrow(idx), m)
  wts <- apply(matrix(w_tmp[idx, 2], nrow(idx), m), 1, prod)
  # scale and recenter the weights and nodes to the specific distribution
  wts <- wts * (pi^{
    -m / 2
  })
  pts <- sqrt(2) * pts


### Different integral approximations

set.seed(648)
res_int <- sapply(1:D, function(d) {

### Calculate domain-specific quantities

V_d_inv <- solve(V_ud_est + V_ed[[d]])

  # Calculate the conditional mean
  mu_ebp <- V_ud_est %*% V_d_inv %*% as.vector(residuals[d, ])

  # Calculate the decomposition of the conditional covariance matrix
  # and ensure it is symmetric
  V_u_V_inv_V_e <- V_ud_est %*% V_d_inv %*% V_ed[[d]]
  V_u_V_inv_V_e[1, 2] <- V_u_V_inv_V_e[2, 1]

  V_u_V_inv_V_e_eig <- eigen(V_u_V_inv_V_e, symmetric = TRUE)
  V_u_V_inv_V_e_eigen_rot <- (V_u_V_inv_V_e_eig$vectors %*% diag(sqrt(V_u_V_inv_V_e_eig$values)))

### Gauss-Hermite

# Adjust points for specific distribution
  pts_distr <- t(V_u_V_inv_V_e_eigen_rot %*% t(pts))
  pts_distr <- sweep(pts_distr, 2, mu_ebp, "+")

  # Calculate the approximation of the rate
  res_GH <- sum(apply(pts_distr, 1, function(x) {
    (x_beta_est[d, ][1] + x[1]) /
      (x_beta_est[d, ][1] + x[1] + x_beta_est[d, ][2] + x[2])
  }) * wts)


### Monte Carlo (without antithetic variate)


# draw random effects
u_tmp <- rmvnorm_rebuilt(n = gauss_n^2, mean = u_bv_BLUP, sigma = V.u_V_inv_V.e)
# apply to original variable (i.e. u) (no antithetic variable)
sgl_MC_values_tmp <- sapply(c(1),
  # s = -1
  function(s) {
    apply(
      u_tmp, 1,
      #       function (x) {
      #  (Xbeta[1] + x[1]) / ( (Xbeta[1] + x[1]) + (Xbeta[2] + x[2]) )
      # }
      f_to_eval_for_f_u_given_y
    )
  },
  simplify = "array"
)
rm(u_tmp)
R_MC_Int_2_wtAntit <- sapply(
  MC_variants,
  # MC_i <- 10
  # MC_i <- MC_variants[1]
  function(MC_i) {
    c(mean(sgl_MC_values_tmp[1:(MC_i), ]))
  }
)


})



### Monte Carlo without antithetic variate


V_ud_est_inv <- solve(V_ud_est)
V_e_inv <- solve(V_e)
V_inv <- solve(V)
solve(V_e_inv + V_u_inv)
V_u_V_inv_V_e <- V_u %*% V_inv %*% V_e

# draw random effects
u_tmp <- rmvnorm_rebuilt(n = gauss_n^2, mean = u_bv_BLUP, sigma = V.u_V_inv_V.e)
# apply to original variable (i.e. u) (no antithetic variable)
sgl_MC_values_tmp <- sapply(c(1),
  # s = -1
  function(s) {
    apply(
      u_tmp, 1,
      #       function (x) {
      #  (Xbeta[1] + x[1]) / ( (Xbeta[1] + x[1]) + (Xbeta[2] + x[2]) )
      # }
      f_to_eval_for_f_u_given_y
    )
  },
  simplify = "array"
)
rm(u_tmp)
R_MC_Int_2_wtAntit <- sapply(
  MC_variants,
  # MC_i <- 10
  # MC_i <- MC_variants[1]
  function(MC_i) {
    c(mean(sgl_MC_values_tmp[1:(MC_i), ]))
  }
)
rm(sgl_MC_values_tmp)
names(R_MC_Int_2_wtAntit) <- paste0("R_MC_Int_2_wtAntit", "_e", total_fct_evals)


### Monte Carlo with antithetic variate

# draw random effects
u_tmp <- rmvnorm_rebuilt(n = MC_variants_antith_max, mean = u_bv_BLUP, sigma = V.u_V_inv_V.e)
# apply to original and antithet. variable (i.e. u and -u)
sgl_MC_values_tmp <- sapply(c(1, -1),
  # s = -1
  function(s) {
    if (s == (-1)) {
      u_tmp <- t(apply(u_tmp, 1, function(x) {
        2 * u_bv_BLUP - x
      }))
    }
    apply(u_tmp, 1, f_to_eval_for_f_u_given_y)
  },
  simplify = "array"
)
R_MC_Int_2 <- sapply(
  MC_variants_antith,
  # MC_i <- 10
  # MC_i <- MC_variants[1]
  function(MC_i) {
    c(mean(sgl_MC_values_tmp[1:MC_i, ]))
  }
)
rm(sgl_MC_values_tmp)
names(R_MC_Int_2) <- paste0("R_MC_Int_2", "_e", total_fct_evals)


### Quasi-random: Halton

fOptions::rnorm.halton(n = gauss_n^2, dimension = 2)
