
# Necessary packages ----
# ______________________________________________________________________________

library(Matrix)
library(mvtnorm)


# Function for covariance matrix ----
# ______________________________________________________________________________

f_Cov <- function(m = 3, theta) {

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
    theta = c(1, 2, 3, .2, .3, .4)
  )
  out_test
}


# Function for preparation of MMFH data generation ----
# ______________________________________________________________________________

f_prep_MMFH_m3 <- function(D = 100,
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
                           seed = 56) {

  # The function prepares the data generation of a MMFH model
  # It generates all quantities which are typically assumed fixed for a small
  # area simulation, like the matrix of auxiliary information x

  ###
  ### Input
  ###

  # D
  #  postive integer. Number of domains
  #
  # v_ref
  #  vector of length 3, positive real values. Variance of random effects
  #  (variable 1, 2, 3)
  #
  # cor_ref
  #  vector of length 3, positive real values in [0,1). Correlation of random
  #  effects (variables 12, 13, 23)
  #
  # v_rer
  #  vector of length 3, positive real values. Variance of sampling errors
  #  (variable 1, 2, 3), same for all domains
  #
  # cor_rer
  #  vector of length 3, positive real values in [0,1). Correlation of sampling
  #  errors (variables 12, 13, 23), same for all domains
  #
  # beta
  #  list of length 3 containing a numeric vector of fixed effects for each
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

  m <- 3 # number of dependent variables

  ###
  ### check input quantities
  ###

  if (length(D) != 1 | abs(D - round(D)) != 0) {
    stop("D should be a
         postive integer. Number of domains")
  }

  if (length(v_ref) != 3 | any(v_ref < 0)) {
    stop("v_ref should be a
         vector of length 3, positive real values.
         Variance of random effects (variable 1, 2, 3)")
  }

  if (length(cor_ref) != 3 | any(cor_ref < 0) | any(cor_ref >= 1)) {
    stop("cor_ref should be a
         vector of length 3, positive real values in [0,1).
         Correlation of random effects (variables 12, 13, 23)")
  }

  if (length(v_rer) != 3 | any(v_rer < 0)) {
    stop("v_rer should be a
         vector of length 3, positive real values.
         Variance of sampling errors (variable 1, 2, 3), same for all domains")
  }

  if (length(cor_rer) != 3 | any(cor_rer < 0) | any(cor_rer >= 1)) {
    stop("cor_rer should be a
         vector of length 3, positive real values in [0,1).
         Correlation of sampling errors (variables 12, 13, 23),
         same for all domains")
  }

  if (length(beta) != 3 | !(is.list(beta))) {
    stop("beta should be a
         list of length 3 containing a numeric vector of fixed effects for each
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
  out_test <- f_prep_MMFH_m3(
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
  str(out_test)
}


# Function for MMFH data generation ----
# ______________________________________________________________________________

f_gen_MMFH_m3 <- function(D = 100,
                          x = x,
                          beta = list(
                            c(1.5, 2.5),
                            c(2.3, 3.3, 4.3),
                            c(4.1, 3.1, 2.1, 2.2)
                          ),
                          V_ud = V_ud,
                          V_ed = V_ed,
                          perc_mis = c(5, 5, 0),
                          seed = 67,
                          verbose = TRUE,
                          return = "2") {
  # This function generates data according to a MMFH model with missing
  # dependent variables (MCAR)
  # Before applying this function, the information which is typically considered
  # fixed in a small area simulation study, can be generated via f_prep_MMFH_m3


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
  # perc_mis
  #  vector of length 3, positive real values in [0,100),
  #  giving the percentages of MCAR values of y for the 3 dependent variables
  #  note that in this version, the missing values are non-overlapping such that
  #  sum(perc_mis) in [0,100)
  #
  # seed
  #  seed for data generation
  #
  # verbose
  #  logical. Whether information of data generation is to be returned
  #
  # return
  #  return = "1" returns all output including the list of sampling error
  #  covariance matrices where the positions of missing dependent variables are
  #  filled by NA,
  #  return = "2" only returns the generated dependent variables

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

  if (length(perc_mis) != 3 | any(perc_mis < 0) | any(perc_mis >= 100) |
    sum(perc_mis) >= 100) {
    stop("perc_mis should be a
       vector of length 3, positive real values in [0,100),
       giving the percentages of MCAR values of y for the three dependent
       variables.
       Note that in this version, the missing values are non-overlapping such
       that sum(perc_mis) in [0,100)")
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
  ### Determine missing values
  ###

  y_mis <- y_obs
  nr_mis <- round(D * perc_mis / 100)
  nr_mis_cumsum <- cumsum(round(D * perc_mis / 100))
  if (max(nr_mis_cumsum) > D) {
    stop("Percentage of missing values is too high")
  }
  if (verbose) {
    cat(
      "\nGenerated number of missing y values (non-overlapping over the", D,
      "domains) for the three dependent variables is", nr_mis, "\n"
    )
  }

  for (k in 1:m) {
    if ((1 + sum(nr_mis_cumsum[k - 1])) < (nr_mis_cumsum[k])) {
      y_mis[(1 + sum(nr_mis_cumsum[k - 1])):(nr_mis_cumsum[k]), k] <- NA
    }
  }

  # Set the corresponding sampling variances missing as well
  V_ed_mis <- V_ed
  for (k in 1:m) {
    y_mis_k_tmp <- which(is.na(y_mis[, k]))
    V_ed_mis[y_mis_k_tmp] <- sapply(y_mis_k_tmp,
      function(d) {
        V_ed_mis[[d]][k, ] <- NA
        V_ed_mis[[d]][, k] <- NA
        V_ed_mis[[d]]
      },
      simplify = FALSE
    )
  }


  ###
  ### Output
  ###

  out <- list(
    y_true,
    y_obs,
    y_mis
  )
  names(out) <- c(
    "y_true",
    "y_obs",
    "y_mis"
  )
  if (return == "1") {
    out["V_ed_mis"] <- V_ed_mis
  }
  return(out)
}


###
### Test function
###

if (FALSE) {
  out_prep_test <- f_prep_MMFH_m3(
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

  out_test <- f_gen_MMFH_m3(
    D = out_prep_test$D,
    x = out_prep_test$x,
    beta = out_prep_test$beta,
    V_ud = out_prep_test$V_ud,
    V_ed = out_prep_test$V_ed,
    perc_mis = c(5, 5, 0),
    seed = 67,
    verbose = TRUE
  )
  str(out_test)
}
