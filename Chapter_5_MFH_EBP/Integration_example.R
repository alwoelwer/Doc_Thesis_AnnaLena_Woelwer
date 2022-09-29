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
  V_ud = d_fix$V_ud,
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
# Apply Gauss-Hermite quadrature ----
# ____________________________________________________________________________

# The approximation of applied to each of the d=1,...,D domains

res_GH <- sapply(1:D, function(d) {

  ### From the model data: Calculate the conditional mean and covariance

  # Calculate the conditional mean
  mu_int2_ebp <- V_ud_est %*% solve(V_ed[[d]] + V_ud_est) %*%
    as.vector(residuals[d, ])

  # Calculate the decomposition of the conditional covariance matrix
  Sigma_Int2_ebp <- V_ud_est %*% solve(V_ud_est + V_ed[[d]]) %*% V_ed[[d]]
  # ensure its symmetric
  Sigma_Int2_ebp[1, 2] <- Sigma_Int2_ebp[2, 1]
  Sigma_Int2_ebp_eig <- eigen(Sigma_Int2_ebp, symmetric = TRUE)
  Sigma_Int2_ebp_eigen_rot <- (Sigma_Int2_ebp_eig$vectors %*%
    diag(sqrt(Sigma_Int2_ebp_eig$values)))


  ### Gauss-Hermite nodes and weights

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
  pts <- t(Sigma_Int2_ebp_eigen_rot %*% t(pts))
  pts <- sweep(pts, 2, mu_int2_ebp, "+")


  ### Calculate the approximation of the rate
  sum(apply(pts, 1, function(x) {
    (x_beta_est[d, ][1] + x[1]) /
      (x_beta_est[d, ][1] + x[1] + x_beta_est[d, ][2] + x[2])
  }) * wts)
})

res_GH
