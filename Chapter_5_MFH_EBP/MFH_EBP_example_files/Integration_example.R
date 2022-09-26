# ______________________________________________________________________________
# ______________________________________________________________________________
#
# Empirical best predictors under multivariate Fay-Herriot models 
# and their numerical approximation ----
#
# Code by Anna-Lena Wölwer in May 2022
# Examples of integral approximations
#
# ______________________________________________________________________________
# ______________________________________________________________________________


# ____________________________________________________________________________
# Load data ----
# ____________________________________________________________________________

# Set user path
path.user = "C:/Users/alwoe/Seafile/me/Uni/Proj/MFH_nonlinear_comb/Sim/Code/GH_example/"

# Laut Gauss-Hermite Weights and Abscissas
gauss_hermite_abscissas_weights <- readRDS(paste0(path.user, "weight_list_normalized_until_n200", ".RData"))


###
### Load some test data
###

(load(file = paste0(path.user, "test_dat_for_GH", ".RData"))  )

# D           Number of domains
# m           Number of dependent variables
# y           Direct estimates (D x m Matrix)
# X_bdiag     Block-diagonal matrix of auxiliary information, dimension: D*m x p
# V_ed        List: Domain-specific Covariance matrices of sampling errors, each list element has dimension m x m
# V_ud_est    Covariance matrix of random effects, estimated from a bivariate FH model, dimension: m x m
# beta_est    Vector of fixed effects, estimated form a bivariate FH model, length: p

# From the test data: calculate, define certain quantities

# Calculate X beta
X_beta_est  <- matrix(X_bdiag %*% beta_est, ncol = m, byrow = TRUE)

# Calculate residuals
residuals   <- y - X_beta_est



# ____________________________________________________________________________
# Apply Gauss-Hermite quadrature ----
# ____________________________________________________________________________

# The approximation of applied to each of the d=1,...,D domains

results <- sapply(1:D, function (d) {
  
  ###
  ### From the model data: Calculate the conditional mean and covariance
  ###
  
  # Calculate the conditional mean
  mu_int2_ebp <- V_ud_est %*% solve(V_ed[[d]] + V_ud_est) %*% as.vector(residuals[d,])

    # calculate the decomposition of the conditional covariance matrix
  Sigma_Int2_ebp <- V_ud_est %*% solve(V_ud_est + V_ed[[d]]) %*% V_ed[[d]]
  # ensure its symmetric
  Sigma_Int2_ebp[1,2] = Sigma_Int2_ebp[2,1]
  Sigma_Int2_ebp_eig       <- eigen(Sigma_Int2_ebp, symmetric = TRUE)
  Sigma_Int2_ebp_eigen_rot <- (Sigma_Int2_ebp_eig$vectors %*% diag(sqrt(Sigma_Int2_ebp_eig$values)))
  
  
  ### 
  ### Get Gauss-Hermite nodes and weights
  ### 
  
  # Choose the number of function evaluations per dimension
  # Choosing gauss_n = 5 results in a total of 25 function evaluations
  gauss_n = 10
  
  # get the corresponding points and weights
  w_tmp <- gauss_hermite_abscissas_weights[[gauss_n]]
  
  # product rule: Set up the grid of points and weights
  idx <- as.matrix(expand.grid(rep(list(1:gauss_n), m)))
  pts <- matrix(w_tmp[idx,1], nrow(idx), m)
  wts <- apply(matrix(w_tmp[idx,2], nrow(idx), m), 1, prod)
  
  # scale and recenter the weights and nodes to the specific distribution
  wts <- wts * (pi^{- m / 2 })
  pts <- sqrt(2) * pts
  pts <- t(Sigma_Int2_ebp_eigen_rot %*% t(pts))
  pts <- sweep(pts, 2, mu_int2_ebp, "+")
  
  
  ### 
  ### Calculate the integral approximation of the unemployment rate
  ### 

  sum(apply(pts, 1, function(x) {
    (X_beta_est[d,][1] + x[1]) / (X_beta_est[d,][1] + x[1] + X_beta_est[d,][2] + x[2])
  }) * wts) 
  
})

# results
# contains the approximation of the unemployment rate for the 
# D domains
