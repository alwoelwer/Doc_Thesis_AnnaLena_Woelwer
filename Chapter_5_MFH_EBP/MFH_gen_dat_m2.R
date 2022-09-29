
# Necessary packages ----
# ______________________________________________________________________________

library(Matrix)
library(mvtnorm)


# Function for covariance matrix ----
# ______________________________________________________________________________

f_Cov <- function(m = 2, theta) {

  # The function returns an m-times-m covariance matrix

  ###
  ### Input
  ###

  # m
  #  positive integer. The returned matrix is of dimension m X m
  #
  # theta
  #  vector of length m*(m+1)/2 containing the variance parameters.
  #  The first m entries correspond to the variances.
  #  The following entries correspond to the correlations and are in [0,1)


  # check input quantities
  if (length(theta) != m * (m + 1) / 2) {
    stop("theta should be a
         vector of length m*(m+1)/2 containing the variance parameters.
         The first m entries correspond to the variances,
         the next entries to the correlations in [0,1)")
  }

  # Covariance matrix
  v_kl <- cbind(rep(1:m, each = m), rep(1:m, times = m))
  v_kl <- v_kl[v_kl[, 1] <= v_kl[, 2], ]
  v_kl <- rbind(
    v_kl[v_kl[, 1] == v_kl[, 2], ],
    v_kl[v_kl[, 1] != v_kl[, 2], ]
  )
  v_kl <- cbind(v_kl, 1:length(theta))

  V_ud <- matrix(nrow = m, ncol = m)
  diag(V_ud) <- theta[1:m]

  for (a in (m + 1):length(theta)) {
    V_ud[v_kl[a, 1], v_kl[a, 2]] <-
      V_ud[v_kl[a, 2], v_kl[a, 1]] <-
      theta[a] * sqrt(theta[v_kl[a, 1]] * theta[v_kl[a, 2]])
  }

  return(V_ud)
}


###
### Test function
###

if (FALSE) {
  out_test <- f_Cov(
    m = 3,
    theta = c(1, 2, 3)
  )
  out_test
}


# Function for preparation of MFH data generation ----
# ______________________________________________________________________________

f_prep_MFH_m2 <- function(D = 100,
                          v_ref = c(2, 3),
                          cor_ref = .7,
                          v_rer = c(2.5, 3.5),
                          cor_rer = .5,
                          beta = list(
                            c(1.5, 2.5),
                            c(2.3, 3.3, 4.3)
                          ),
                          range_aux = c(10, 100),
                          seed = 56) {

  # The function prepares the data generation of a 
  # multivariate Fay-Herriot (MFH) model with m=2 dependent variables.
  # It generates all quantities which are typically assumed fixed for a small
  # area simulation, like the matrix of auxiliary information x

  ###
  ### Input
  ###

  # D
  #  postive integer. Number of domains
  #
  # v_ref
  #  vector of length 2, positive real values. Variance of random effects
  #  (variable 1, 2)
  #
  # cor_ref
  #  positive real value in [0,1). Correlation of random effects of variable 1 and 2
  #
  # v_rer
  #  vector of length 2, positive real values. Variance of sampling errors
  #  (variable 1, 2), same for all domains
  #
  # cor_rer
  #  positive real value in [0,1). Correlation of sampling errors of variable 1 and 2
  #
  # beta
  #  list of length 2 containing a numeric vector of fixed effects for each
  #  dependent variable. The first value is for the intercept
  #
  # range_aux
  #  vector of length 2, real values in ascending order.
  #  The range of values from which the uniformly distributed auxiliary data are
  #  drawn the dimension of the generated auxiliary information is in
  #  accordance with beta
  #
  # seed
  #  seed for data generation

  m <- 2 # number of dependent variables

  ###
  ### check input quantities
  ###

  if (length(D) != 1 | abs(D - round(D)) != 0) {
    stop("D should be a
         postive integer. Number of domains")
  }

  if (length(v_ref) != 2 | any(v_ref < 0)) {
    stop("v_ref should be a
         vector of length 2, positive real values.
         Variance of random effects (variable 1, 2)")
  }

  if (length(cor_ref) != 2 | any(cor_ref < 0) | any(cor_ref >= 1)) {
    stop("cor_ref should be a
         vector of length 2, positive real values in [0,1).
         Correlation of random effects variable 1 to 2")
  }

  if (length(v_rer) != 2 | any(v_rer < 0)) {
    stop("v_rer should be a
         vector of length 2, positive real values.
         Variance of sampling errors (variable 1, 2), same for all domains")
  }

  if (length(cor_rer) != 2 | any(cor_rer < 0) | any(cor_rer >= 1)) {
    stop("cor_rer should be a
         vector of length 2, positive real values in [0,1).
         Correlation of sampling errors of variable 1 and 2,
         same for all domains")
  }

  if (length(beta) != 2 | !(is.list(beta))) {
    stop("beta should be a
         list of length 2 containing a numeric vector of fixed effects for each
         dependent variable")
  }

  if (length(range_aux) != 2 | any(range_aux < 0) |
    range_aux[1] >= range_aux[2]) {
    stop("range_aux should be a
         vector of length 2, real values in ascending order.
         The range of values from which the uniformly distributed auxiliary
         data are drawn")
  }


  ###
  ### covariance matrices
  ###

  # variance-covariance matrix of sampling errors for the D domains
  # (in this example code the covariance matrix is the same for all domains)
  cov_mat_rer <- f_Cov(m = m, theta = c(v_rer, cor_rer))
  V_ed <- lapply(1:D, function(a) {
    cov_mat_rer
  })

  # variance-covariance matrix of random effects
  V_ud <- f_Cov(m = m, theta = c(v_ref, cor_ref))


  ###
  ### auxiliary information
  ###

  p_k <- sapply(beta, length) # number of fixed effects per dependent variable
  p <- sum(p_k) # total number of fixed effects

  # generate auxiliary information x

  x <- list()
  for (k in 1:m) {
    x[[k]] <- sapply(1:p_k[k],
      function(k_i) {
        if (k_i == 1) {
          rep(1, D) # Intercept as first column
        } else {
          sapply(
            1:D,
            function(d) {
              runif(
                n = 1,
                min = range_aux[1],
                max = range_aux[2]
              )
            }
          )
        }
      },
      simplify = "array"
    )
  }

  ###
  ### Output
  ###

  out <- list(
    m,
    D,
    v_ref,
    cor_ref,
    v_rer,
    cor_rer,
    beta,
    range_aux,
    seed,
    m,
    V_ed,
    V_ud,
    x
  )
  names(out) <- c(
    "m",
    "D",
    "v_ref",
    "cor_ref",
    "v_rer",
    "cor_rer",
    "beta",
    "range_aux",
    "seed",
    "m",
    "V_ed",
    "V_ud",
    "x"
  )
  return(out)
}


###
### Test function
###

if (FALSE) {
  out_test <- f_prep_MFH_m3(
    D = 100,
    v_ref = c(2, 3),
    cor_ref = .7,
    v_rer = c(2.5, 3.5),
    cor_rer = .5,
    beta = list(
      c(1.5, 2.5),
      c(2.3, 3.3, 4.3)
    ),
    range_aux = c(10, 100),
    seed = 56
  )
  str(out_test)
}


# Function for MFH data generation ----
# ______________________________________________________________________________

f_gen_MFH_m3 <- function(D = 100,
                         x = x,
                         beta = list(
                           c(1.5, 2.5),
                           c(2.3, 3.3, 4.3),
                           c(4.1, 3.1, 2.1, 2.2)
                         ),
                         V_ud = V_ud,
                         V_ed = V_ed,
                         seed = 67,
                         verbose = TRUE) {
  # This function generates data according to a 
  # multivariate Fay-Herriot (MFH) model with m=3 dependent variables.
  # Before applying this function, the information which is typically considered
  # fixed in a small area simulation study, can be generated via f_prep_MFH_m3


  ###
  ### Input
  ###

  # D
  #  positive integer. Number of domains
  #
  # beta
  #  list of length 3 containing a numeric vector of fixed effects
  #  for each dependent variable. The first value is for the intercept
  #
  # x
  #  list of of length 3. Each list element contains the auxiliary
  #  data for all three dependent variables and D domains, first column should
  #  be an intercept column (consisting of 1s) the dimensions have to fit those
  #  of beta
  #
  # V_ud
  # 3-times-3 covariance matrix of random effects
  #
  # V_ed
  #  list of length D, each containing a 3-times-3 covariance
  #               matrix of sampling errors
  #
  # seed
  #  seed for data generation
  #
  # verbose
  #  logical. Whether information of data generation is to be returned

  m <- 3 # number of dependent variables

  ###
  ### check input matrices
  ###

  if (length(D) != 1 | abs(D - round(D)) != 0) {
    stop("D should be a
         postive integer. Number of domains")
  }

  if (length(beta) != 3 | !(is.list(beta))) {
    stop("beta should be a
         list of length 3 containing a numeric vector of fixed effects for each
         dependent variable")
  }


  ###
  ### Transform input
  ###

  x_bdiag <- do.call(
    "rbind",
    lapply(
      1:D,
      function(d) {
        Matrix::Matrix(do.call(
          Matrix::bdiag,
          lapply(
            x,
            function(x_d) {
              matrix(x_d[d, ],
                nrow = 1
              )
            }
          )
        ),
        sparse = TRUE
        )
      }
    )
  )
  x_beta <- matrix(x_bdiag %*% unlist(beta), ncol = 3, byrow = TRUE)


  ###
  ### Generate model data
  ###

  set.seed(seed)
  u <- as.array(mvtnorm::rmvnorm(n = D, mean = rep(0, m), sigma = V_ud))
  e <- t(sapply(1:D, # sampling errors
    function(d) {
      as.vector(mvtnorm::rmvnorm(
        n = 1, mean = rep(0, m),
        sigma = V_ed[[d]]
      ))
    },
    simplify = "array"
  ))
  y_true <- x_beta + u
  y_obs <- y_true + e


  ###
  ### Output
  ###

  out <- list(
    y_true,
    y_obs
  )
  names(out) <- c(
    "y_true",
    "y_obs"
  )
  return(out)
}


###
### Test function
###

if (FALSE) {
  out_prep_test <- f_prep_MFH_m3(
    D = 100,
    v_ref = c(2, 3, 4),
    cor_ref = c(.2, .3, .4),
    v_rer = c(2.5, 3.5, 4.5),
    cor_rer = c(.15, .25, .35),
    beta = list(
      c(1.5, 2.5),
      c(2.3, 3.3, 4.3),
      c(4.1, 3.1, 2.1, 2.2)
    ),
    range_aux = c(10, 100),
    seed = 56
  )
  str(out_prep_test)
  names(out_prep_test)

  out_test <- f_gen_MFH_m3(
    D = out_prep_test$D,
    x = out_prep_test$x,
    beta = out_prep_test$beta,
    V_ud = out_prep_test$V_ud,
    V_ed = out_prep_test$V_ed,
    seed = 67,
    verbose = TRUE
  )
  str(out_test)
}
