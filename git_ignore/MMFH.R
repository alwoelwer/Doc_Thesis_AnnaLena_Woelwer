f_MMFH <- function(y, 
                   x, 
                   V_ed, 
                   theta = NULL, 
                   theta_start = NULL, 
                   method = "REML", 
                   eps = 1e-8, 
                   maxiter = 500, 
                   verbose = FALSE){
  
  ### Notation
  # d = 1, .., D 1:D of interest
  # m variables of interest
  # k, l commonly used indices for the dependent variables(k <= m)
  
  ### Input
  #
  # y 
  # D X m matrix of observed values for(y_1, , ..., y_m) in each 
  # domain d = 1, ..., D. 
  # Missing values are filled by NA.
  #
  # x 
  # list of length m containing the D X p matrix of covariates for 
  # variable(y_1, ..., y_m) with the same ordering as in y.
  # Missing values are not allowed.
  # The first list element refers to the first column of y and so forth.
  #
  # V_ed 
  # list of length D giving the covariance matrix of(y_1, ..., y_m) 
  # for each domain d = 1, ..., D with the same ordering as in y.
  # Missing values for missing direct estimates can be filled by NA.
  #
  # theta 
  # vector of known variance components \sigma_u1^2, \sigma_u2^2, rho. 
  # If set =NULL(default), the variance components are estimated.
  # The order is:
  # The first m(number variables of interest) elements 
  # correspond to the random effects variances.
  # The(m+1):length(theta) elements correspond to the 
  # random effects correlations.
  #
  # theta_start 
  # can give starting values for theta(if theta is unknown) 
  # which are otherwise determined by univariate estimates.
  #
  # eps 
  # precision needed for assuming convergence.
  #
  # maxiter 
  # maximum number of iterations.
  #
  # verbose 
  # logical. TRUE gives some numerical diagnostic output.
  
  # Note for x:
  # X only works with at least one auxiliary variable for each variable y.
  # If no intercept is given as first column, the intercept is added.
  
  
  # Load Packages ----
  # ____________________________________________________________________________
  
  library(sae)
  library(Matrix)
  
  
  # Check input ----
  # ____________________________________________________________________________
  
  if(any(rowSums(apply(y, 1, is.na) * 1) == dim(y)[1])){
    tmp_v1 <- names(which(rowSums(apply(y, 1, is.na)) == dim(y)[1]))
    stop(paste0("y: All values of variable ", tmp_v1, " are missing.
 Exclude variable and run again"))
  }
  m <- dim(y)[2]
  
  
  ###
  ### Check mode of input
  ###
  
  if(mode(y) != "numeric"){
    y <- as.matrix(y)
    warning(paste0(
      "mode(y) is not numeric.", 
      "y has to be a D X m matrix of observed values for(y_1, ..., y_m) in each 
 domain d=1, ..., D and is therefore tranformed with as.matrix(y).\n"
    ))
  }
  
  if(mode(x) != "list"){
    stop(paste0(
      "mode(x) is not a list.", 
      "x has to be list of length m containing the D X p matrix of ", 
      "covariates for variable(y_1, ..., y_m). ", 
      "It has to be of same order as y and is therefore transformed.\n"
    ))
  }
  
  if(any(sapply(x, function(x_k){ mode(x_k) }) %in% c("list"))){
    warning(paste0(
      "At least one mode of list elements of x is a list instead of numeric. ", 
      "x has to be list of length m containing the D X p matrix of ", 
      "covariates for variable(y_1, ..., y_m). ", 
      "The list elements are therefore transformed via as.matrix().\n"
    ))
    x <- lapply(x, function(x_k){ as.matrix(x_k) })
  }
  
  if(any(unlist(lapply(x, is.matrix)) == FALSE)){
    stop(paste0(
      "The list elements of x have to be a D X p matrix of covariates. ", 
      "They are therefore transformed via as.matrix().\n"
    ))
    x <- lapply(x, function(x_k){ as.matrix(x_k) })
  }
  
  if(mode(V_ed) != "list"){
    stop(paste0(
      "mode(V_ed) is not a list.", 
      "V_ed has to be a list of length D giving the m x m covariance matrix of
 the direct estimates of(y_1, ..., y_m) for each domain d=1, ..., D ", 
      "in same order as y.", 
      "Where direct estimates are missing, V_ed can contain NA values.\n"
    ))
  }
  
  if(any(sapply(V_ed, function(V_ed_d){ mode(V_ed_d) }) %in% c("list"))){
    warning(paste0(
      "The list elements of V_ed are lists.", 
      "They have to be covariance matrices ", 
      "of(y_1, ..., y_m) for each domain d=1, ..., D. ", 
      "The list elements of V_ed are therefore transformed via as.matrix().\n"
    ))
    V_ed <- lapply(V_ed, function(V_ed_d){ as.matrix(V_ed_d) })
  }
  
  
  ###
  ### Check dimensions of input
  ###
  
  if((length(unique(unlist(lapply(x, nrow)))) != 1) | 
     nrow(y) != unlist(lapply(x, nrow))[1]){
    stop(paste0(
      "nrow(y) must equal nrow() of all list elements of x."
    ))
  }
  
  
  if(nrow(y) != length(V_ed)){
    stop(paste0(
      "nrow(y) must equal length(V_ed)."
    ))
  }
  
  if((length(unique(unlist(lapply(V_ed, dim)))) != 1) | 
     (unlist(lapply(V_ed, dim))[1] != ncol(y))){
    stop(paste0(
      "ncol(y), i.e. the number of variables of interest, must equal ", 
      " nrow() and ncol() of each list element of V_ed."
    ))
  }
  
  # check if there are missing values in x where y is observed
  if(anyNA(x)){
    stop(paste0(
      "At least one missing value found in auxiliary variables x.\n"
    ))
  }
  
  # insert intercept if not already in x
  x <- lapply(x, function(x_k){
    
    tmp_v1 <- apply(x_k, 2, function(x_k_k){ range(x_k_k) })
    tmp_v2 <- apply(tmp_v1, 2, function(x_k_k){
      if(x_k_k[1] == x_k_k[2] & x_k_k[1] != 1){
        stop(paste0("There is a constant column in auxiliary ", 
                    "data x which is not an intercept, i.e. it is != 1"))
      }
      x_k_k[1] == x_k_k[2] & x_k_k[1] == 1
    }) 
    if(sum(tmp_v2) > 1){
      stop(paste0("There are two intercept columns in x."))
    }
    if(sum(tmp_v2) != 1){
      x_k <- cbind("(Intercept)" = 1, x_k)
    }
    x_k
  }) 
  
  
  ###
  ### same dim names for y, x, and V_ed
  ###
  
  if(is.null(dimnames(y)[[1]])){ 
    dimnames(y)[[1]] <- 1:nrow(y) 
  }
  
  if(is.null(dimnames(y)[[2]])){ 
    dimnames(y)[[2]] <- paste0("v", 1:m)
  }
  
  x <- lapply(1:m, function(k){
    dimnames(x[[k]])[[1]] <- dimnames(y)[[1]]
    dimnames(x[[k]])[[2]] <- paste0("v", k, ".", "x", 
                                    1:dim(x[[k]])[2])
    x[[k]]
  })
  names(x) <- dimnames(y)[[2]]
  names(V_ed) <- dimnames(y)[[1]]
  
  
  ###
  ### completely missing observations in y are omitted
  ###
  
  obs.comp.mis <- rownames(y)[rowSums(is.na(y)) == dim(y)[2]]
  if(!identical(obs.comp.mis, character(0))){
    warning("There were 1:D with no direct estimate in any variable of interest.
 These 1:D are excluded from the estimation.")
    y <- y[!rownames(y) %in% obs.comp.mis, ]
    x <- lapply(x, function(x_k){
      x_k <- x_k[!rownames(x_k) %in% obs.comp.mis, ]
      x_k
    })
    V_ed <- V_ed[!names(V_ed) %in% obs.comp.mis]
  }
  
  
  ###
  ### method and start values
  ###
  
  if(method != "REML" & method != "ML"){
    stop(paste0("method=\"", method, "\" must be \"REML\" or \"ML\".\n"))
  }
  
  if(!is.null(theta)){
    if(length(theta) !=(m*(m+1)/2)){
      stop(paste0("Length of theta must be equal to ", m*(m+1)/2, 
                  " or =NULL.\n"))
    }
    theta_fix <- TRUE
  } else {
    theta_fix <- FALSE
  }
  
  if(!is.null(theta_start) & theta_fix){
    warning(paste0(
      "theta and starting values for theta are provided at the same time. ", 
      "Only theta will be considered.\n"
    ))
    theta_start <- NULL
  }
  
  if(!theta_fix){
    if(!is.null(theta_start)){
      if(length(theta_start) !=(m*(m+1)/2)){
        warning(paste0(
          "Length of theta_start must be equal to ",(m*(m+1)/2), 
          ". Therefore, theta_start is ignored.\n"
        ))
        theta_start_fix <- FALSE
      } else {
        if(any(theta_start[1:m] <= 0)){
          warning(paste0(
            "Negative or zero starting values for theta[", 1, ":", m, "]. ", 
            "Therefore, theta_start is ignored.\n"
          ))
          theta_start_fix <- FALSE
        } else {
          if(abs(theta_start[(m+1):length(theta_start)]) >= 1){
            warning(paste0("Starting value greater or equal 1 for ", 
                           "theta[",(m+1), ":", length(theta_start), "]. ", 
                           "theta_start is ignored.\n"
            ))
            theta_start_fix <- FALSE
          } else {
            theta_start_fix <- TRUE
          }
        }
      }
    } else {
      theta_start_fix <- FALSE
    }
  }
  
  # Functions ----
  # ____________________________________________________________________________
  
  #if(MMFH.verbose.steps){ print("Define functions") }
  
  ###
  # Create matrix V_ud from theta vector
  f_V_ud <- function(theta, order_ab){
    
    V_ud <- matrix(nrow = m, ncol = m)
    
    # variances
    diag(V_ud) <- theta[1:m]
    
    # covariances
    for(theta_a in(m + 1):length(theta)){ 
      kl_a <- order_ab[theta_a, ]
      
      V_ud[kl_a$k, kl_a$l] <- 
        V_ud[kl_a$l, kl_a$k] <- 
        theta[kl_a$theta] * sqrt(theta[kl_a$k] * 
                                   theta[kl_a$l]) 
    }
    V_ud
  }
  
  # if(FALSE){
  # # test
  # tmp_v1 = abs(rnorm(n = dim(order_ab)[1]))
  # tmp_v1[(m+1):dim(order_ab)[1]] <- sample(c(.2,.3,.4), 
  # length((m+1):dim(order_ab)[1]), 
  # replace = TRUE)
  # f_V_ud(theta = tmp_v1, order_ab)
  # # test
  # tmp_v1[6] * sqrt(tmp_v1[2] * tmp_v1[3]) == f_V_ud(theta = tmp_v1, order_ab)[2, 3]
  # tmp_v1[2] == f_V_ud(theta = tmp_v1, order_ab)[2, 2]
  # }
  
 
  
  # if(FALSE){
  # f_theta_m(theta = theta_iter1, order_ab)
  # theta_iter1
  # theta <- theta_iter1
  # }
  
  
  ###
  # Create covariance matrix of y, given by V_ud + V_ed
  f_V <- function(theta, V_ed, D, order_ab){
    V_ud <- f_V_ud(theta, order_ab)
    unname(Matrix::Matrix(do.call(Matrix::bdiag, lapply(1:D, function(d){ 
      V_ed[[d]] + V_ud 
    })), 
    sparse = TRUE))
  }
  
  # f_V(theta, V_ed, D, order_ab)
  # dimension: D*m x D*m
  # Note:
  # V_inv_obs: solve(V) of current V, take only rows and columns for observed direct estimates
  # V_inv_x_obs: solve(V) %*% X = V_inv_obs %*% X or current V
  if(FALSE){
    f_V(theta = theta_iter1, V_ed, D, order_ab)
    V_inv_obs <- tryCatch(expr = solve(f_V(theta = theta_iter1, V_ed, D, order_ab)[I_obs_vec, I_obs_vec]), 
                          error = function(e){ ginv(f_V(theta = theta_iter1, V_ed, D, order_ab)[I_obs_vec, I_obs_vec]) })
    # dimension: length of observed variables of intrest ^2 
    V_inv_x_obs <- crossprod(V_inv_obs, x_bdiag_obs)
    # dimension: length of observed variables of intrest x p(nr. auxiliary variables)
  }
  
  ###
  # Create matrix Q
  f_Q <- function(x_bdiag_obs, V_inv_x_obs){
    tmp_v1 <- t(x_bdiag_obs) %*% V_inv_x_obs
    tryCatch(expr = solve(tmp_v1), 
             error = function(e){ ginv(tmp_v1) })
  }
  
  if(FALSE){
    f_Q(x_bdiag_obs, V_inv_x_obs)
    Q_obs <- f_Q(x_bdiag_obs, V_inv_x_obs)
  }
  
  ###
  # Create matrix P
  f_P <- function(x_bdiag_obs, Q_obs, V_inv_obs, V_inv_x_obs){
    V_inv_obs - V_inv_x_obs %*% Q_obs %*% t(x_bdiag_obs) %*% V_inv_obs
  }
  
  if(FALSE){
    # are the same
    all.equal(
      (diag(nrow(x_bdiag_obs)) - V_inv_x_obs %*% Q_obs %*% t(x_bdiag_obs)) %*% V_inv_obs, 
      V_inv_obs - V_inv_x_obs %*% Q_obs %*% t(x_bdiag_obs) %*% V_inv_obs
    )
    # which is faster?
    {
      a1_beg_tmp <- Sys.time()
      for(a in 1:50){
        (diag(nrow(x_bdiag_obs)) - V_inv_x_obs %*% Q_obs %*% t(x_bdiag_obs)) %*% V_inv_obs
      }
      a1_end_tmp <- Sys.time()
      a2_beg_tmp <- Sys.time()
      for(a in 1:50){
        V_inv_obs - V_inv_x_obs %*% Q_obs %*% t(x_bdiag_obs) %*% V_inv_obs
      }
      a2_end_tmp <- Sys.time()
      cat(a1_end_tmp - a1_beg_tmp, "\n")
      cat(a2_end_tmp - a2_beg_tmp)
    }
    # second approach is faster
  }
  
  ###
  # First partial derivatives of V_ud w.r.t. theta
  f_V_ud_l <- function(theta, order_ab){
    
    ### 
    # Create matrix of theta vector
    f_theta_m <- function(theta, order_ab){
      
      theta_mat <- matrix(nrow = m, ncol = m)
      # variances
      diag(theta_mat) <- theta[1:m]
      # correlations
      for(theta_a in(m+1):length(theta)){ 
        tmp_v1 <- order_ab[theta_a, ]
        theta_mat[tmp_v1$k, tmp_v1$l] <- 
          theta_mat[tmp_v1$l, tmp_v1$k] <- 
          theta[tmp_v1$theta] 
      }
      theta_mat
    }

    # for all elements of theta
    V_ud_l <- lapply(1:length(theta), function(theta_a){
      V_ud_l_tmp <- mapply(function(k, l){
        tmp_v1 <- 0
        
        # derivatives with respect to variances
        if(theta_a %in% 1:m){
          
          # diagonal
          if(k == l & k == theta_a){ 
            tmp_v1 <- 1 
          } 
          
          # off-diagonal
          if((k != l & k == theta_a)){ 
            tmp_v1 <-(theta_mat[k, l] * 
                        sqrt(theta_mat[l, l])) / 
              (2*sqrt(theta_mat[k, k]))
          }
          if((k != l & l == theta_a)){ 
            tmp_v1 <-(theta_mat[k, l] * 
                        sqrt(theta_mat[k, k])) / 
              (2*sqrt(theta_mat[l, l]))
          }
        } else {
          
          # derivatives with respect to correlations
          if((order_ab$k[theta_a] == k & order_ab$l[theta_a] == l) | 
             (order_ab$k[theta_a] == l & order_ab$l[theta_a] == k)){ 
            tmp_v1 <- sqrt(theta_mat[k, k]) * 
              sqrt(theta_mat[l, l])
          } 
        }
        tmp_v1
      }, order_a, order_b)
      V_ud_l_tmp <- matrix(V_ud_l_tmp, 
                           ncol = m, 
                           nrow = m, 
                           byrow = TRUE)
      V_ud_l_tmp
    })
    V_ud_l
  }
  
  if(FALSE){
    f_V_ud_l(theta = theta_iter1, order_ab)
  }
  
  ###
  # List of sparse block diagonal matrices each of dimension D*m X D*m
  # of derivatives of V_ud w.r.t. theta
  f_V_ud_l_mDmD <- function(theta, D, order_ab){
    
    V_ud_d_tmp_v1 <- f_V_ud_l(theta, order_ab)
    
    lapply(1:length(theta), function(theta_a){ # theta_a=1
      do.call(Matrix::bdiag, lapply(1:D, function(d) {
        V_ud_d_tmp_v1[[theta_a]]
      }))
    })
  }
  
  if(FALSE){
    f_V_ud_l(theta = theta_iter1, order_ab)[[1]]
    f_V_ud_l_mDmD(theta = theta_iter1, 1, order_ab)[[1]]
    # f_V_ud_l_mDmD(theta = theta_iter1, c(1), order_ab)
    # f_V_ud_l_mDmD(theta = theta_iter1, c(150), order_ab)
  }
  
  
  ###
  # ML log-likelihood
  f_ML_ll <- function(y, x_bdiag, beta, theta, V_ed, order_ab, O_d, 
                      pos_y_obs_d){
    
    I_obs_mat <- !is.na(y)
    m_d <- rowSums(I_obs_mat)
    m <- ncol(y)
    resids <- y - matrix(x_bdiag %*% beta, ncol = m, byrow = TRUE)
    #rownames(resids) <- 1:D
    V_ud <- f_V_ud(theta, order_ab)
    
    # for all domains d
    tmp_v1 <- sapply(1:D, function(d){ 
      
      res_d_obs <- resids[d, pos_y_obs_d[[d]]]
      O_d_d <- O_d[[d]]
      
      # if only one variable is observed in domain d
      if(!is.matrix(O_d_d)){
        V_d_obs <- as.numeric(O_d_d %*% (V_ed[[d]] + V_ud) %*% 
                                t(matrix(O_d_d, nrow = 1)))
        V_d_obs_inv <- V_d_obs^(-1)
        
        - m_d[d]/2 * log(2*pi) -
          .5 * log(V_d_obs) -
          .5 * as.numeric(V_d_obs_inv * res_d_obs^2)
        
        # if at least two variables are observed in domain d
      } else {
        V_d_obs <- O_d_d %*% (V_ed[[d]] + V_ud) %*% t(O_d_d) 
        V_d_obs_inv <- tryCatch(expr = solve(V_d_obs), 
                                error = function(e){ ginv(V_d_obs) })
        
        - m_d[d]/2 * log(2*pi) -
          .5 * determinant(V_d_obs, logarithm = TRUE)$modulus[1] -
          .5 * as.numeric(t(res_d_obs) %*% V_d_obs_inv %*% res_d_obs)
      }
    })
    sum(tmp_v1)
  }
  
  if(FALSE){
    beta <- beta_new
    d = "3"
    theta = theta_iter1
    theta = theta_new
    V_inv_obs <- tryCatch(expr = solve(f_V(theta = theta_iter1, V_ed, D, order_ab)[I_obs_vec, I_obs_vec]), 
                          error = function(e){ ginv(f_V(theta = theta_iter1, V_ed, D, order_ab)[I_obs_vec, I_obs_vec]) })
    # dimension: length of observed variables of intrest ^2 
    V_inv_x_obs <- crossprod(V_inv_obs, x_bdiag_obs)
    beta <- as.vector(tryCatch(expr = solve(a =(t(x_bdiag_obs) %*% V_inv_x_obs), b =(t(V_inv_x_obs) %*% Y_obs)), 
                               error = function(e){ ginv(t(x_bdiag_obs) %*% V_inv_x_obs) %*% (t(V_inv_x_obs) %*% Y_obs) })
    )
    f_ML_ll(y, x_bdiag, beta, theta = theta_iter1, V_ed, order_ab, O_d, pos_y_obs_d)
    
    
    # is equal
    all.equal(log(det(V_d_obs)), 
              determinant(V_d_obs, logarithm = TRUE)$modulus[1])
  }
  
  ###
  # ML log-likelihood in dependence of step size lambda
  f_ML_ll_lambda <- function(lambda, y, x_bdiag, beta, theta_old, theta_delta, 
                            V_ed, order_ab, O_d, pos_y_obs_d){
    f_ML_ll(y, x_bdiag, beta, 
            theta = as.vector(theta_old + lambda * theta_delta), 
            V_ed, order_ab, O_d, pos_y_obs_d)
  }
  
  if(FALSE){
    f_ML_ll_lambda(lambda = 1, y, x_bdiag, beta, 
                  theta_old = theta_iter1, 
                  theta_delta = rep(0, length(theta_iter1)), 
                  V_ed, m, order_ab, O_d, pos_y_obs_d)
  }
  
  ###
  # ML score function theta
  f_S_theta_ML <- function(y, x_bdiag, beta, theta, order_ab, O_d, V_ed){
    
    m <- dim(y)[2] 
    resids <- y - matrix(x_bdiag %*% beta, ncol = m, byrow = TRUE)
    #rownames(resids) <- 1:D
    V_ud <- f_V_ud(theta, order_ab)
    V_ud_l <- f_V_ud_l(theta, order_ab)
    
    # for all elements of theta
    sapply(1:length(theta), function(a){
      
      V_ud_d_a <- V_ud_l[[a]]
      
      # for all domains d
      tmp_v1 <- sapply(1:D, function(d){ 
        
        pos_y_obs_d_d <- pos_y_obs_d[[d]]
        theta_d <- order_ab$theta[order_ab$k %in% pos_y_obs_d_d & order_ab$l %in% pos_y_obs_d_d]
        res_d_obs <- resids[d, pos_y_obs_d[[d]]]
        O_d_d <- O_d[[d]]
        
        # if theta_a does not correspond to the observed y values in domain d
        if(!a %in% theta_d){
          0
          
          # if theta_a corresponds to the observed y values in domain d
        } else {
          
          # if only one variable is observed in domain d
          if(!is.matrix(O_d_d)){
            V_d_obs <- as.numeric(O_d_d %*% (V_ed[[d]] + V_ud) %*% 
                                    t(matrix(O_d_d, nrow = 1)))
            V_d_obs_inv <- V_d_obs^(-1)
            
            - .5 * V_d_obs_inv +
              .5 * V_d_obs_inv^2 * res_d_obs^2
            
            # if at least two variables are observed in domain d
          } else {
            V_d_obs <- O_d_d %*% (V_ed[[d]] + V_ud) %*% t(O_d_d) 
            V_d_obs_inv <- tryCatch(expr = solve(V_d_obs), 
                                    error = function(e){ ginv(V_d_obs) })
            V_ud_d_a_obs <- O_d_d %*% V_ud_d_a %*% t(O_d_d)
            
            - .5 * sum(diag(V_d_obs_inv %*% V_ud_d_a_obs)) +
              .5 * t(res_d_obs) %*% V_d_obs_inv %*% V_ud_d_a_obs %*% 
              V_d_obs_inv %*% res_d_obs
          }
        }
        
      })
      sum(tmp_v1)
    })
  }
  
  if(FALSE){
    f_S_theta_ML(y, x_bdiag, beta, theta, order_ab, O_d, V_ed)
  }
  
  ###
  # ML score function beta
  f_S_beta_ML <- function(y, x_bdiag, x_bdiag_obs_list, beta, theta, order_ab, 
                          O_d, V_ed){
    
    m <- dim(y)[2] 
    resids <- y - matrix(x_bdiag %*% beta, ncol = m, byrow = TRUE)
    #rownames(resids) <- 1:D
    V_ud <- f_V_ud(theta, order_ab)
    
    # for all domains d
    tmp_v1 <- sapply(1:D, function(d){ 
      
      res_d_obs <- resids[d, pos_y_obs_d[[d]]]
      O_d_d <- O_d[[d]]
      
      # if only one variable is observed in domain d
      if(!is.matrix(O_d_d)){
        V_d_obs <- O_d_d %*% (V_ed[[d]] + V_ud) %*% t(matrix(O_d_d, nrow = 1))
        V_d_obs_inv <- V_d_obs^(-1)
        as.vector(
          t(matrix(x_bdiag_obs_list[[d]], 1)) %*% V_d_obs_inv %*% res_d_obs)
        
        # if at least two variables are observed in domain d
      } else {
        V_d_obs <- O_d_d %*% (V_ed[[d]] + V_ud) %*% t(O_d_d) 
        V_d_obs_inv <- tryCatch(expr = solve(V_d_obs), 
                                error = function(e){ ginv(V_d_obs) })
        as.vector(t(x_bdiag_obs_list[[d]]) %*% V_d_obs_inv %*% res_d_obs)
      }
    })
    rowSums(tmp_v1)
  }
  
  if(FALSE){
    f_S_beta_ML(y, x_bdiag, beta, theta, order_ab, V_ed)
  }
  
  ###
  # ML Fisher information function theta
  f_F_theta_ML <- function(y, x_bdiag, beta, theta, order_ab, O_d, V_ed){
    
    m <- dim(y)[2] 
    resids <- y - matrix(x_bdiag %*% beta, ncol = m, byrow = TRUE)
    #rownames(resids) <- 1:D
    V_ud <- f_V_ud(theta, order_ab)
    V_ud_l <- f_V_ud_l(theta, order_ab)
    
    theta_l <- length(theta)
    theta_iters_a <- rep(1:theta_l, times = theta_l)
    theta_iters_b <- rep(1:theta_l, each = theta_l)
    
    # for all elements of theta
    Flist <- mapply(function(theta_a, theta_b){
      
      V_ud_d_a <- V_ud_l[[theta_a]]
      V_ud_d_b <- V_ud_l[[theta_b]]
      
      # for all domains d
      tmp_v1 <- sapply(1:D, function(d){ 
        
        O_d_d <- O_d[[d]]
        pos_y_obs_d_d <- pos_y_obs_d[[d]]
        theta_d <- order_ab$theta[order_ab$k %in% pos_y_obs_d_d & 
                                    order_ab$l %in% pos_y_obs_d_d]
        
        # if theta_a does not correspond to the observed y values in domain d
        if(!theta_a %in% theta_d){
          0
          
          # if theta_a corresponds to observed y values in domain d
        } else {
          
          # if only one variable is observed in domain d
          if(!is.matrix(O_d_d)){
            V_d_obs <- O_d_d %*% (V_ed[[d]] + V_ud) %*% t(matrix(O_d_d, 1))
            V_d_obs_inv <- V_d_obs^(-1)
            
            .5 *(V_d_obs_inv^2)
            
            # if at least two variables are observed in domain d
          } else {
            V_d_obs <- O_d_d %*% (V_ed[[d]] + V_ud) %*% t(O_d_d) 
            V_d_obs_inv <- tryCatch(expr = solve(V_d_obs), 
                                    error = function(e){ ginv(V_d_obs) })
            V_ud_d_a_obs <- O_d_d %*% V_ud_d_a %*% t(O_d_d)
            V_ud_d_b_obs <- O_d_d %*% V_ud_d_b %*% t(O_d_d)
            
            .5 * sum(diag(V_d_obs_inv %*% V_ud_d_a_obs %*% V_d_obs_inv %*% 
                            V_ud_d_b_obs))
          }
        }
      })
      sum(tmp_v1)
    }, theta_iters_a, theta_iters_b)
    matrix(Flist, length(theta), length(theta))
  }
  
  if(FALSE){
    f_F_theta_ML(y, x_bdiag, beta, theta, order_ab, V_ed)
  }
  
  ###
  # ML Fisher information function beta
  f_F_beta_ML <- function(y, x_bdiag, x_bdiag_obs_list, beta, theta, order_ab, 
                          O_d, V_ed){
    
    m <- dim(y)[2] 
    resids <- y - matrix(x_bdiag %*% beta, ncol = m, byrow = TRUE)
    #rownames(resids) <- 1:D
    V_ud <- f_V_ud(theta, order_ab)
    
    # for all elements of beta
    # for all domains d
    tmp_v1 <- lapply(1:D, function(d){ 
      O_d_d <- O_d[[d]]
      
      # if only one variable is observed in domain d
      if(!is.matrix(O_d_d)){
        V_d_obs <- O_d_d %*% (V_ed[[d]] + V_ud) %*% t(matrix(O_d_d, nrow = 1))
        V_d_obs_inv <- V_d_obs^(-1)
        
        t(matrix(x_bdiag_obs_list[[d]], nrow = 1)) %*% V_d_obs_inv %*% 
          x_bdiag_obs_list[[d]]
        
        # if at least two variables are observed in domain d
      } else {
        V_d_obs <- O_d_d %*% (V_ed[[d]] + V_ud) %*% t(O_d_d) 
        V_d_obs_inv <- tryCatch(expr = solve(V_d_obs), 
                                error = function(e){ ginv(V_d_obs) })
        
        t(x_bdiag_obs_list[[d]]) %*% V_d_obs_inv %*% x_bdiag_obs_list[[d]]
      }
    })
    Reduce('+', tmp_v1)
  }
  
  if(FALSE){
    f_F_beta_ML(y, x_bdiag, beta, theta, order_ab, V_ed)
  }
  
  
  ###
  # log likelihood
  f_REML_ll <- function(Y_obs, x_bdiag_obs, theta, V_ed, I_mis_vec){
    
    V_obs <- f_V(theta, V_ed, D, order_ab)[I_obs_vec, I_obs_vec]
    V_inv_obs <- tryCatch(expr = solve(V_obs), 
                          error = function(e){ ginv(V_obs) })
    V_inv_x_obs <- crossprod(V_inv_obs, x_bdiag_obs)
    Q_obs <- f_Q(x_bdiag_obs, V_inv_x_obs)
    P_obs <- f_P(x_bdiag_obs, Q_obs, V_inv_obs, V_inv_x_obs)
    V_ud <- f_V_ud(theta, order_ab)
    
    # see https://en.wikipedia.org/wiki/Block_matrix#Block_diagonal_matrices 
    # Let A be a block-diagonal matrix with blocks A_i i=1, ..., m
    # Then det(A) = det(A_1) x ... x det(A_m) = prod_i=1^m(det(A_i))
    # 
    # sum of log(determinants(single blocks)) instead of
    # log(determinant(complete block-diagonal matrix)), 
    # i.e.
    # sum(sapply(1:D, function(d){ determinant(V_ed[[d]] + V_ud, logarithm = TRUE)$modulus[1] }))
    # instead of 
    # determinant(V_obs, logarithm = TRUE)$modulus[1]
    
    # old 
    # 
    # as.vector(.5 * determinant(V_obs, logarithm = TRUE)$modulus[1] -
    # .5 * determinant(t(x_bdiag_obs) %*% V_inv_x_obs, logarithm = TRUE)$modulus[1] -
    # .5 * t(Y_obs) %*% P_obs %*% Y_obs)
    
    tmp_v1 <- sum(sapply(1:D, function(d){
      b_tmp <-(V_ed[[d]] + V_ud)[I_obs_mat[d, ], I_obs_mat[d, ]]
      if(is.null(dim(b_tmp))){
        log(b_tmp)
      } else {
        determinant(b_tmp, logarithm = TRUE)$modulus[1]
      }}))
    
    as.vector(- .5 * tmp_v1 -
                .5 * determinant(t(x_bdiag_obs) %*% V_inv_x_obs, 
                                 logarithm = TRUE)$modulus[1] -
                .5 * t(Y_obs) %*% P_obs %*% Y_obs)
  }
  
  if(FALSE){
    # log(det(block-digonal matrix))
    theta <- theta_iter1
    V_obs <- f_V(theta, V_ed, D, order_ab)[I_obs_vec, I_obs_vec]
    determinant(V_obs, logarithm = TRUE)$modulus[1]
    
    # sum log(det(blocks))
    V_ud <- f_V_ud(theta, order_ab)
    sum(sapply(1:D, function(d){
      b_tmp <-(V_ed[[d]] + V_ud)[I_obs_mat[d, ], I_obs_mat[d, ]]
      if(is.null(dim(b_tmp))){
        log(b_tmp)
      } else {
        determinant(b_tmp, logarithm = TRUE)$modulus[1]
      }
    }))
    
    #
    theta <- theta_new
    V_obs <- f_V(theta, V_ed, D, order_ab)[I_obs_vec, I_obs_vec]
    V_obs[1:6, 1:6]
    V_obs[100:105, 100:105]
    V_inv_obs[1:6, 1:6]
    f_REML_ll(Y_obs, x_bdiag_obs, theta, V_ed, I_mis_vec)
    
    # test if equal
    determinant(t(x_bdiag_obs) %*% V_inv_x_obs, logarithm = TRUE)$modulus[1]
    log(det(t(x_bdiag_obs) %*% V_inv_x_obs)) 
    
    determinant(t(x_bdiag_obs) %*% V_inv_x_obs, logarithm = TRUE)$modulus[1]
    log(abs(det(t(x_bdiag_obs) %*% V_inv_x_obs)))
    
    all.equal(determinant(V_obs, logarithm = TRUE)$modulus[1], 
              log(det(V_obs)))
    
    # calculation of log determinant
    
    A <- unname(t(x_bdiag_obs) %*% V_inv_x_obs)
    log(det(A)) 
    log(prod(eigen(A)$values))
    log(eigen(A)$values)
    nlme::logDet(as.matrix(A)) 
    determinant(A, logarithm = TRUE)$modulus[1]
  }
  
  ###
  # REML log-likelihood in dependence of step size lambda
  f_REML_ll_lambda <- function(lambda, Y_obs, x_bdiag_obs, theta_old, V_ed, 
                              I_mis_vec, theta_delta){
    f_REML_ll(Y_obs, x_bdiag_obs, 
              theta = as.vector(theta_old + lambda * theta_delta), 
              V_ed, I_mis_vec)
  }
  
  if(FALSE){
    theta_old <- theta_iter1
    f_REML_ll_lambda(lambda = 1, Y_obs, x_bdiag_obs, theta_old, V_ed, I_mis_vec, theta_delta = rep(0, length(theta_old)))
  }
  
  ###
  # REML score vector
  f_S_REML <- function(P_obs, Y_obs, V_ud_l_mDmD_obs, theta_l){
    sapply(1:theta_l, function(theta_a){
      V_ud_l_mDmD_obs_theta_a <- V_ud_l_mDmD_obs[[theta_a]]
      as.vector(-0.5 * sum(diag(P_obs %*% V_ud_l_mDmD_obs_theta_a)) +
                  0.5 * t(Y_obs) %*% P_obs %*% V_ud_l_mDmD_obs_theta_a %*% 
                  P_obs %*% Y_obs)
    })
  }
  
  if(FALSE){
    theta <- theta_iter1
    V_inv_obs <- tryCatch(expr = solve(f_V(theta = theta, V_ed, D, order_ab)[I_obs_vec, I_obs_vec]), 
                          error = function(e){ ginv(f_V(theta = theta, V_ed, D, order_ab)[I_obs_vec, I_obs_vec]) })
    V_inv_x_obs <- crossprod(V_inv_obs, x_bdiag_obs)
    Q_obs <- f_Q(x_bdiag_obs, V_inv_x_obs)
    P_obs <- f_P(x_bdiag_obs, Q_obs, V_inv_obs, V_inv_x_obs)
    V_ud_l_mDmD <- f_V_ud_l_mDmD(theta = theta*1.1, D, order_ab)
    V_ud_l_mDmD_obs <- lapply(1:theta_l, function(theta_a){ V_ud_l_mDmD[[theta_a]][I_obs_vec, I_obs_vec] })
    
    f_S_REML(P_obs, Y_obs, V_ud_l_mDmD_obs, theta_l)
  }
  
  ###
  # REML Fisher information matrix
  f_F_REML <- function(P_obs, V_ud_l_mDmD_obs, theta_l){
    
    theta_iters_a <- rep(1:theta_l, times = theta_l)
    theta_iters_b <- rep(1:theta_l, each = theta_l)
    
    Flist <- mapply(function(theta_a, theta_b){
      V_ud_l_mDmD_obs_theta_a <- V_ud_l_mDmD_obs[[theta_a]]
      V_ud_l_mDmD_obs_theta_b <- V_ud_l_mDmD_obs[[theta_b]]
      
      0.5 * sum(diag(P_obs %*% V_ud_l_mDmD_obs_theta_a %*% P_obs %*% 
                       V_ud_l_mDmD_obs_theta_b))
    }, theta_iters_a, theta_iters_b)
    matrix(Flist, theta_l, theta_l)
  }
  
  if(FALSE){
    theta <- theta_iter1
    V_inv_obs <- tryCatch(expr = solve(f_V(theta = theta, V_ed, D, order_ab)[I_obs_vec, I_obs_vec]), 
                          error = function(e){ ginv(f_V(theta = theta, V_ed, D, order_ab)[I_obs_vec, I_obs_vec]) })
    V_inv_x_obs <- crossprod(V_inv_obs, x_bdiag_obs)
    Q_obs <- f_Q(x_bdiag_obs, V_inv_x_obs)
    P_obs <- f_P(x_bdiag_obs, Q_obs, V_inv_obs, V_inv_x_obs)
    V_ud_l_mDmD <- f_V_ud_l_mDmD(theta = theta_new, D, order_ab)
    V_ud_l_mDmD_obs <- lapply(1:theta_l, function(theta_a){ V_ud_l_mDmD[[theta_a]][I_obs_vec, I_obs_vec] })
    
    f_F_REML(P_obs, V_ud_l_mDmD_obs, theta_l)
  }
  
  ###
  # Adjust the Fisher-Scoring step size delta
  # If the algorithm uses many iterations: sequentially reduce step size
  # Furthermore: check whether the parameter estimates under delta are
  # feasible. If not: choose optimal delta for which estimates are feasible
  
  # f_lambda_opt <- function(iter, method, theta_old, theta_delta,
  #                          Y_obs, x_bdiag_obs, V_ed, I_mis_vec,
  #                          y, x_bdiag, beta, m, order_ab, O_d, 
  #                          pos_y_obs_d){
  #   
  #   theta_l <- length(theta_old)
  #   
  #   # The possible step sizes for lambda are more and more restricted
  #   lambda_bounds <- data.frame("names" = levels(cut(20,
  #                                              breaks = c(0, 
  #                                                         seq(20, 200, by = 20)),
  #                                              right = TRUE)),
  #                         "lower" = 0,
  #                         "upper" = c(1, 1, (9:2)/10))
  #   iter_class_i <- cut(iter,
  #                       breaks = c(0, seq(20, 200, by = 20)),
  #                       right = TRUE)
  #   lower <- lambda_bounds[lambda_bounds$names == iter_class_i, "lower"]
  #   upper <- lambda_bounds[lambda_bounds$names == iter_class_i, "upper"]
  #   
  #   # check if resulting values under lambda = upper are feasible
  #   tmp_v3 <- as.vector(theta_old + upper * theta_delta)
  #   
  #   # If theta estimates under lambda = upper are out of border: 
  #   # search for optimal lambda in bounds of feasible values
  #   if(any(tmp_v3[1:m] <= 0) | any(abs(tmp_v3[(m+1):theta_l]) >= 1)){
  #     I_tmp <- TRUE
  #     lambda <- f_lambda_opt_feas(iter, method, m, theta_old, theta_delta)
  #     
  #     # if resulting values under lambda = upper are feasible
  #   } else {
  #     if(method == "REML"){
  #       lambda <- optim(par = upper, 
  #                       fn = f_REML_ll_lambda, 
  #                       lower = lower, 
  #                       upper = upper, 
  #                       method = "L-BFGS-B", 
  #                       Y_obs = Y_obs, 
  #                       x_bdiag_obs = x_bdiag_obs, 
  #                       theta_old = theta_old, 
  #                       V_ed = V_ed, 
  #                       I_mis_vec = I_mis_vec, 
  #                       theta_delta = theta_delta, 
  #                       control = list(fnscale = -1))$par
  #     }
  #     if(method == "ML"){
  #       lambda <- optim(par = upper, 
  #                       fn = f_ML_ll_lambda, 
  #                       lower = lower, 
  #                       upper = upper, 
  #                       method = "L-BFGS-B", 
  #                       y = y, 
  #                       x_bdiag = x_bdiag, 
  #                       beta = beta_new, 
  #                       theta_old = theta_old, 
  #                       theta_delta = theta_delta, 
  #                       V_ed = V_ed, 
  #                       m = m, 
  #                       order_ab = order_ab, 
  #                       O_d = O_d, 
  #                       pos_y_obs_d = pos_y_obs_d, 
  #                       control = list(fnscale = -1))$par
  #     }
  #   }
  #   lambda
  # }
  
  
  ###
  # it may appear that new values are out of boundaries
  # function searches for best lambda within allowed boundaries
  # of parameters
  
  f_lambda_opt_feas <- function(method, m, theta_old, theta_delta,
                                lower, Y_obs, x_bdiag_obs, V_ed, 
                                I_mis_vec, y, x_bdiag, beta, order_ab, 
                                O_d, pos_y_obs_d){
    
    theta_l <- length(theta_old)
    upper <- c(sapply(1:m, function(theta_a){ 
      abs(theta_old[theta_a] / theta_delta[theta_a]) 
    }), 
    sapply((m+1):theta_l, function(theta_a){
      ((1 - theta_old[theta_a]) / theta_delta[theta_a]) 
    }), 
    sapply((m+1):theta_l, function(theta_a){
      ((-1 - theta_old[theta_a]) / theta_delta[theta_a]) 
    }))
    
    # only accept delta values in [0,1]
    upper2 <- upper[upper >= 0 & upper <= 1]
    if(length(upper2) > 1){ upper3 <- sort(upper2, decreasing = FALSE)[1] }
    if(length(upper2) == 0){ upper3 <- 0 }
    if(length(upper2) == 1){ upper3 <- upper2 }
    if(length(upper3) == 0){ lambda <- 0 }
    if((upper3 > 1) |(upper3 == 0)){
      lambda <- 0
    } else {
      upper <- upper3 * 0.99
      if(upper < lower){ lower <- upper * .5 }
      if(method == "REML"){
        lambda <- optim(par = upper, 
                        fn = f_REML_ll_lambda, 
                        lower = lower, 
                        upper = upper, 
                        method = "L-BFGS-B", 
                        Y_obs = Y_obs, 
                        x_bdiag_obs = x_bdiag_obs, 
                        theta_old = theta_old, 
                        V_ed = V_ed, 
                        I_mis_vec = I_mis_vec, 
                        theta_delta = theta_delta, 
                        control = list(fnscale = -1))$par
      }
      if(method == "ML"){
        lambda <- optim(par = upper, 
                        fn = f_ML_ll_lambda, 
                        lower = lower, 
                        upper = upper, 
                        method = "L-BFGS-B", 
                        y = y, 
                        x_bdiag = x_bdiag, 
                        beta = beta_new, 
                        theta_old = theta_old, 
                        theta_delta = theta_delta, 
                        V_ed = V_ed, 
                        order_ab = order_ab, 
                        O_d = O_d, 
                        pos_y_obs_d = pos_y_obs_d, 
                        control = list(fnscale = -1))$par
      }
    }
    lambda
  }
  
  ###
  # Print results of single Fisher-Scoring iterations
  f_FS_print <- function(beta_new, theta_new, ll, lambda, 
                         dig_theta, dig_beta, dig_lambda, dig_ll, dig_iter, 
                         wid_theta, wid_beta, wid_lambda, wid_ll, wid_iter){
    if(iter == 0){
      verb1 <- c(formatC(c("iter"), width = wid_iter), 
                 formatC(c("loglike"), width = wid_ll), 
                 formatC(c("lambda"), width = wid_lambda), 
                 formatC(c(theta_names), width = wid_theta), 
                 formatC(c(paste0("", beta_names)), width = wid_beta))
      cat(verb1, "\n")
    }
    verb2 <- c(formatC(c(iter), format = "f", digits = dig_iter, 
                       width = wid_iter), 
               formatC(c(round(ll, digits = dig_ll)), format = "f", 
                       digits = dig_ll, width = wid_ll), 
               formatC(c(as.character(round(lambda, digits = dig_lambda))), 
                       format = "f", digits = dig_lambda, width = wid_lambda), 
               formatC(c(round(theta_new, digits = dig_theta)), format = "f", 
                       digits = dig_theta, width = wid_theta), 
               formatC(c(round(beta_new, digits = dig_beta)), format = "f", 
                       digits = dig_beta, width = wid_beta))
    cat(verb2, "\n")
  }
  
  if(FALSE){
    dig_theta <- 3
    dig_beta <- 3
    dig_lambda <- 3
    dig_ll <- 3
    dig_iter <- 0
    
    wid_theta <- range(nchar(round(theta_new, digits = dig_theta)))[2]; if(wid_theta < 5){ wid_theta <- 5 }
    wid_beta <- range(nchar(round(beta_new, digits = dig_beta)))[2]; if(wid_beta < 5){ wid_beta <- 5 }
    wid_lambda <- 6
    wid_ll <- range(nchar(round(ll, digits = dig_ll)))[2]; if(wid_ll < 18){ wid_ll <- 18 }
    wid_iter <- 5
    f_FS_print(beta_new, theta_new, ll, lambda, 
               dig_theta, dig_beta, dig_lambda, dig_ll, dig_iter, 
               wid_theta, wid_beta, wid_lambda, wid_ll, wid_iter)
  }
  
  
  # Transformation of input quantities ----
  # ____________________________________________________________________________
  
  ###
  ### variables of interest
  ###
  
  #1:D <- rownames(y)
  D <- dim(y)[1] # number of 1:D
  
  # indicators for missing values
  I_obs_mat <- !is.na(y)
  pos_y_obs_d <- lapply(1:nrow(I_obs_mat), function(a){ 
    unname(which(I_obs_mat[a, ] == TRUE)) 
  })
  pos_y_mis_d <- lapply(1:nrow(!I_obs_mat), function(a){ 
    unname(which(!I_obs_mat[a, ] == TRUE)) 
  })
  # names(pos_y_obs_d) <- rownames(I_obs_mat)
  # names(pos_y_mis_d) <- rownames(!I_obs_mat)
  
  # y as vector
  Y <- as.vector(t(y))
  I_obs_vec <- !is.na(Y)
  I_mis_vec <- is.na(Y)
  Y_obs <- Y[I_obs_vec]
  #names(Y_obs) <- paste0(rep(1:D, each = m), ".", "v", rep(1:m, times = D))[I_obs_vec]
  
  
  # for each domain: 
  # diagonal matrix of indicator whether variable is observed
  # where all rows equal to zero are deleted such that
  # dimension: number of observed variables x m
  O_d <- lapply(1:D, function(d){
    tmp_v1 <- diag((I_obs_mat * 1)[d, ])
    tmp_v1 <- tmp_v1[!rowSums(tmp_v1) == 0, ]
    tmp_v1
  })
  #names(O_d) <- 1:D
  
  
  ###
  ### auxiliary data
  ###
  
  p_k <- unlist(lapply(x, function(x_k){ 
    dim(x_k)[2] 
  }))
  #names(p_k) <- 1:m
  p <- sum(p_k)
  beta_assign_k <- unlist(lapply(1:m, function(k){ 
    rep(k, times = p_k[k]) 
  }))
  x_bdiag <- do.call("rbind", 
                     lapply(1:D, function(d){
                       Matrix::Matrix(do.call(Matrix::bdiag, 
                                              lapply(x, function(x_d){
                                                matrix(x_d[d, ], nrow = 1)
                                              })), sparse = TRUE)
                     }))
  
  
  # x_bdiag <- do.call("rbind", 
  # lapply(1:D, 
  # function(d){ 
  # Matrix::Matrix(do.call(Matrix::bdiag, 
  # lapply(x, 
  # function(x_d){ 
  # matrix(x_d[d, ], nrow = 1) 
  # })), sparse = TRUE) 
  # }))
  
  ### note: bdiag cannot handle lists as arguments on server
  # works
  # lapply(1:D, function(d){ do.call(Matrix::bdiag, lapply(x, function(x_k){ matrix(x_k[d, ], nrow = 1) })) })
  # works not
  # lapply(1:D, function(d){ as.matrix(Matrix::bdiag(lapply(x, function(x_k){ matrix(x_k[d, ], nrow = 1) }))) })
  
  #x_bdiag <- do.call("rbind", x_bdiag)
  
  # dimnames(x_bdiag) <- list("d.v" = paste0("d", rep(1:D, each = m), ".", "v", rep(1:m, times = D)), 
  # "v.x" = paste0(unlist(lapply(x, 
  # function(x_k){
  # colnames(x_k)
  # })))
  #)
  
  beta_names <- unlist(lapply(x, function(x_k){ 
    colnames(x_k) 
  }))
  beta_assign_k_mat <- sapply(1:m, function(k){ 
    beta_assign_k == k 
  })
  #colnames(beta_assign_k_mat) <- 1:m 
  
  # x_bdiag only for observed variables
  x_bdiag_obs_list <- lapply(1:D, function(d){ 
    Matrix::Matrix(do.call(Matrix::bdiag, 
                           lapply(x, function(x_k){ # x_k <- x[[1]]
                             matrix(x_k[d, ], nrow = 1) 
                           })), sparse = TRUE)[pos_y_obs_d[[d]], ]
    # matrix(do.call(Matrix::bdiag, 
    # lapply(x, 
    # function(x_k){ 
    # matrix(x_k[d, ], nrow = 1) 
    # }))[pos_y_obs_d[[d]], ], ncol = p) 
  })
  #names(x_bdiag_obs_list) <- 1:D
  
  # check dimensions
  # sum(matrix(unlist(lapply(x_bdiag_obs_list, dim)), byrow = TRUE, ncol = 2)[, 1]) == dim(x_bdiag_obs)[1]
  
  x_bdiag_obs <- x_bdiag[I_obs_vec, ]
  if(det(t(x_bdiag_obs) %*% x_bdiag_obs) %in% c(-Inf, Inf)
     #| det((t(x_bdiag) %*% x_bdiag)) < 10^(-10) | det((t(x_bdiag) %*% x_bdiag)) > 10^(10)
  ){
    try(expr = solve(t(x_bdiag_obs) %*% x_bdiag_obs))
    stop(paste0("x is not of full rank. Please transform auxiliary data x.\n"))
  }
  
  
  ###
  ### variance components
  ### 
  
  # theta
  # assignment which element of theta belongs to which variable
  order_a <- rep(1:m, each = m)
  order_b <- rep(1:m, times = m)
  order_ab <- cbind(order_a, order_b)
  order_ab <- order_ab[order_ab[, 1] <= order_ab[, 2], ]
  order_ab <- rbind(order_ab[order_ab[, 1] == order_ab[, 2], ], 
                    order_ab[order_ab[, 1] != order_ab[, 2], ])
  colnames(order_ab) <- c("k", "l")
  order_ab <- cbind(order_ab, "theta" = 1:dim(order_ab)[1])
  order_ab <- data.frame(order_ab)
  theta_l <- m *(m + 1) / 2
  theta_names <- c(paste0("var_", 1:m), 
                   paste0("rho_", 
                          apply(order_ab[(m+1):theta_l, ], 1, function(theta_a){ 
                            paste0(theta_a[1], theta_a[2]) 
                          })))
  
  
  ###
  ### If necessary: Adjust covariance matrix of sampling errors
  ### 
  
  #if(MMFH.verbose.steps){ print("Check input: variance covariance matrix of sampling errors") }
  
  # where y is missing, V_ed should be 0
  V_ed <- lapply(1:D, function(d){
    if(anyNA(V_ed[[d]][pos_y_obs_d[[d]], pos_y_obs_d[[d]]])){
      stop(paste0(
        "Sampling variances or covariances for y is missing in V_ed, 
 at least for domain ", d, ".\n"
      ))
    }
    V_ed_tmp <- V_ed[[d]]
    V_ed_tmp[pos_y_mis_d[[d]], ] <- 0
    V_ed_tmp[, pos_y_mis_d[[d]]] <- 0
    V_ed_tmp
  })
  #names(V_ed) <- 1:D
  
  
  # Initial parameter estimates for Fisher-Scoring ----
  # ____________________________________________________________________________
  
  # As starting values: Take the solutions of single (univariate) FH models
  
  if(is.null(theta)){
    if(!theta_start_fix){
      
      theta_iter1 <- vector(length = theta_l)
      
      for(k in 1:m){
        obs_k_d <- which(I_obs_mat[, k] == TRUE)
        tmp_mod <- tryCatch(
          expr = sae::eblupFH(formula = y[obs_k_d, k] ~ - 1 + x[[k]][obs_k_d, ], 
                              vardir = c(unlist(lapply(obs_k_d, function(d){ 
                                V_ed[[d]][k, k] 
                              }))), 
                              method = method, 
                              MAXITER = maxiter, 
                              PRECISION = eps), 
          error = function(e){ "Error in univariate estimation" }) 
        if(identical(tmp_mod, "Error in univariate estimation")){
          warning(paste0(
            "Search for starting values: \n", 
            "Univariate model(function sae::eblupFH()) for variable ", k, 
            " did not converge.\n", 
            "Therefore, the starting value of ", 
            "the respective variance parameter is set ", 
            "to 1.\n"
          ))
          tmp_mod <- list()
          tmp_mod$fit <- list()
          tmp_mod$fit$refvar <- 1
        }
        if(abs(tmp_mod$fit$refvar) > 10^(-5)){
          theta_iter1[k] <- tmp_mod$fit$refvar
        } else {
          theta_iter1[k] <- 1
        }
        rm(tmp_mod, obs_k_d)
      }
      kl_tmp <- order_ab[order_ab$k < order_ab$l, ]
      cor_kl_tmp <- apply(kl_tmp, 1, function(kl_tmp_i){
        obs_k_d <- !is.na(y[, c(kl_tmp_i[1], kl_tmp_i[2])])
        b_tmp <- rownames(obs_k_d[rowSums(obs_k_d) == 2, ])
        if(length(b_tmp) > 2){
          cor(y[b_tmp, kl_tmp_i[1]], y[b_tmp, kl_tmp_i[2]]) 
        } else {
          0 
        }
      })
      theta_iter1[(m+1):length(theta_iter1)] <- cor_kl_tmp
    } else {
      theta_iter1 <- theta_start
    }
  } else {
    theta_iter1 <- theta
  }
  if(exists("cor_kl_tmp")){ rm(cor_kl_tmp) }
  if(exists("kl_tmp")){ rm(kl_tmp) }
  names(theta_iter1) <- theta_names
  
  
  # Preparation of Fisher Scoring ----
  # ____________________________________________________________________________
  
  #if(MMFH.verbose.steps){ print("Preparation of Fisher Scoring") }
  
  # initial values
  iter <- 0
  converged <- TRUE
  lambda <- NA
  I_tmp <- FALSE # indicating if FS step size was restricted due to in-feasible values
  theta_new <- theta_iter1
  #names(theta_new) <- theta_names
  V_obs <- f_V(theta = theta_new, V_ed, D, order_ab)[I_obs_vec, I_obs_vec]
  V_inv_obs <- tryCatch(expr = solve(V_obs), 
                        error = function(e){ ginv(V_obs) })
  V_inv_x_obs <- V_inv_obs %*% x_bdiag_obs
  a_tmp <- t(x_bdiag_obs) %*% V_inv_x_obs
  b_tmp <- t(V_inv_x_obs) %*% Y_obs
  beta_new <- as.vector(tryCatch(
    expr = solve(a = a_tmp, b = b_tmp), 
    error = function(e){ ginv(a_tmp %*% b_tmp) }))
  #names(beta_new) <- beta_names
  
  # log likelihood of estimates
  if(method == "ML"){ ll <- f_ML_ll(y, x_bdiag, beta = beta_new, theta = theta_new, V_ed, order_ab, O_d, pos_y_obs_d) }
  if(method == "REML"){ ll <- f_REML_ll(Y_obs, x_bdiag_obs, theta = theta_new, V_ed, I_mis_vec) }
  if(identical(NaN, ll)){ stop("NaNs poduced for the likelihood. Check matrix of auxiliary data.") }
  
  if(verbose){
    # prepare verbose output
    dig_theta <- 3
    dig_beta <- 3
    dig_lambda <- 3
    dig_ll <- 3
    dig_iter <- 0
    
    wid_theta <- range(nchar(round(theta_new, digits = dig_theta)))[2]; if(wid_theta < 6){ wid_theta <- 6 }
    wid_beta <- range(nchar(round(beta_new, digits = dig_beta)))[2]; if(wid_beta < 6){ wid_beta <- 6 }
    wid_lambda <- 6
    wid_ll <- range(nchar(round(ll, digits = dig_ll)))[2]; if(wid_ll < 18){ wid_ll <- 18 }
    wid_iter <- 5
    
    f_FS_print(beta_new, theta_new, ll, lambda, 
               dig_theta, dig_beta, dig_lambda, dig_ll, dig_iter, 
               wid_theta, wid_beta, wid_lambda, wid_ll, wid_iter)
  }
  
  # set parameters
  theta_old <- .5 * theta_new
  beta_old <- .5 * beta_new
  
  # step size
  lambda <- 1
  lambda_bounds <- data.frame("names" = levels(cut(20,
                                                   breaks = c(0, 
                                                              seq(20, 200, by = 20)),
                                                   right = TRUE)),
                              "lower" = 0,
                              "upper" = c(1, 1, (9:2)/10))
  
  
  # Fisher Scoring ----
  # ____________________________________________________________________________
  
  while((any(abs((beta_new - beta_old) / beta_old) > eps) |
         any(abs((theta_new - theta_old) / theta_old) > eps)) &
        iter <= maxiter){
    
    iter <- iter + 1
    beta_old <- beta_new
    theta_old <- theta_new
    
    if(!theta_fix){
      
      ### REML
      if(method == "REML"){
        V_obs <- f_V(theta = theta_old, V_ed, D, order_ab)[I_obs_vec, I_obs_vec]
        V_inv_obs <- tryCatch(expr = solve(V_obs), 
                              error = function(e){ ginv(V_obs) })
        V_inv_x_obs <- V_inv_obs %*% x_bdiag_obs
        Q_obs <- f_Q(x_bdiag_obs, V_inv_x_obs)
        P_obs <- f_P(x_bdiag_obs, Q_obs, V_inv_obs, V_inv_x_obs)
        V_ud_l_mDmD <- f_V_ud_l_mDmD(theta = theta_new, D, order_ab)
        V_ud_l_mDmD_obs <- lapply(1:theta_l, function(theta_a){ 
          V_ud_l_mDmD[[theta_a]][I_obs_vec, I_obs_vec] 
        })
        F_REML <- f_F_REML(P_obs, V_ud_l_mDmD_obs, theta_l)
        S_REML <- f_S_REML(P_obs, Y_obs, V_ud_l_mDmD_obs, theta_l)
        theta_delta <- as.vector(tryCatch(
          expr = solve(a = F_REML, b = S_REML), 
          error = function(e){ ginv(F_REML %*% S_REML) }) )
      } 
      
      ### ML
      if(method == "ML"){
        F_beta <- f_F_beta_ML(y, x_bdiag, x_bdiag_obs_list, beta = beta_old, 
                              theta = theta_old, order_ab, O_d, V_ed)
        S_beta <- f_S_beta_ML(y, x_bdiag, x_bdiag_obs_list, beta = beta_old, 
                              theta = theta_old, order_ab, O_d, V_ed)
        beta_new <- beta_old + as.vector(tryCatch(
          expr = solve(a = F_beta, b = S_beta), 
          error = function(e){ ginv(F_beta) %*% S_beta }) )
        F_theta <- f_F_theta_ML(y, x_bdiag, beta = beta_new, theta = theta_old, 
                                order_ab, O_d, V_ed)
        S_theta <- f_S_theta_ML(y, x_bdiag, beta = beta_new, theta = theta_old, 
                                order_ab, O_d, V_ed)
        theta_delta <- as.vector(tryCatch(
          expr = solve(a = F_theta, b = S_theta), 
          error = function(e){ ginv(F_theta) %*% S_theta }) )
      } 
      
      ###
      ### choice of step size lambda
      ###
      
      # In first iteration ll usually goes down, but even more up in next
      # iterations.
      # Therefore, if lambda is chosen optimally in first iteration, 
      # it is too often set = zero.
      # Hence, if possible, lambda = 1 in first rounds
      
      # In further iterations, lambda can be chosen optimally, because
      # elsewise often too many iterations
      
      iter_class_i <- cut(iter,
                          breaks = c(0, seq(20, 200, by = 20)),
                          right = TRUE)
      lower <- lambda_bounds[lambda_bounds$names == iter_class_i, "lower"]
      upper <- lambda_bounds[lambda_bounds$names == iter_class_i, "upper"]
      lambda <- upper
      
      # check if resulting values are feasible
      # If estimated variance components are out of border: 
      # search for optimal lambda in bounds of feasible values
      theta_cand <- as.vector(theta_old + lambda * theta_delta)
      if(any(theta_cand[1:m] <= 0) | any(abs(theta_cand[(m+1):theta_l]) >= 1)){
        I_tmp <- TRUE
        lambda <- f_lambda_opt_feas(method, m, theta_old, theta_delta,
                                    lower, Y_obs, x_bdiag_obs, V_ed, 
                                    I_mis_vec, y, x_bdiag, beta = beta_new, 
                                    order_ab, 
                                    O_d, pos_y_obs_d)
      }
      # lambda <- f_lambda_opt(iter, method, theta_old, theta_delta,
      #                        Y_obs, x_bdiag_obs, V_ed, I_mis_vec,
      #                        y, x_bdiag, beta, m, order_ab, O_d, 
      #                        pos_y_obs_d)
      
      theta_new <- as.vector(theta_old + lambda * theta_delta)
      theta_new1 <- theta_new
      
      ###
      ### plot log likelihood depending on possible lambda
      ###
      
      # if(FALSE){
      #   rbind(theta.true, theta_old, theta_delta)
      #   upper <- 1
      #   deltas_tmp <- seq(0, upper, .01)
      #   names(deltas_tmp) <- seq(0, upper, .01)
      #   if(method == "REML"){
      #     ll_lambdas_tmp <- sapply(deltas_tmp, function(a){
      #       f_REML_ll_lambda(lambda = a, 
      #                       Y_obs, x_bdiag_obs, 
      #                       theta_old, V_ed, 
      #                       I_mis_vec, theta_delta) })
      #   }
      #   if(method == "ML"){
      #     ll_lambdas_tmp <- sapply(deltas_tmp, function(a){
      #       f_ML_ll_lambda(lambda = a, 
      #                     y, x_bdiag, beta = beta_new, 
      #                     theta_old, theta_delta, 
      #                     V_ed, m, order_ab, 
      #                     O_d, pos_y_obs_d) })
      #   }
      #   plot(deltas_tmp, ll_lambdas_tmp)
      #   sort(ll_lambdas_tmp, decreasing = TRUE)[1:10]
      #   rm(deltas_tmp, ll_lambdas_tmp)
      # }
      
      ### Check whether parameter estimates are within bounds
      
      if(any(theta_new[1:m] <= 10^(-9))){
        warning(paste0(
          "Negative or zero ", 
          if(method == "REML"){"restricted "}, 
          "maximum-likelihood estimate for sigma_u^2. ", 
          "Taking values of previous iteration.\n"
        ))
        theta_new <- theta_old
        converged <- FALSE
        break;
      }
      
      if(any(abs(theta_new[(m+1):theta_l]) >= 1)){
        warning(paste0(
          if(method == "REML"){"Restricted m"} else {"M"}, 
          "maximum-likelihood estimate greater 1 for rho. ", 
          "Taking values of previous iteration.\n"
        ))
        theta_new <- theta_old
        converged <- FALSE
        break;
      }
      
      
      # ### See whether estimates alternate by comparing last iterations
      # ### if they alternate, choose middle value to break alternation
      # 
      # if(converged){
      #   if(iter == 1){
      #     theta_iters <- matrix(theta_new, 
      #                           byrow = TRUE, 
      #                           ncol = length(theta_old))
      #     colnames(theta_iters) <- theta_names
      #   } else {
      #     theta_iters <- rbind(theta_iters, theta_new)
      #   }
      #   if(iter > 5){
      #     theta_iters <- theta_iters[2:6, ]
      #     #if results alternate, choose middle values
      #     if(prod(abs((theta_iters[5, ] - theta_iters[3, ]) / 
      #                 theta_iters[3, ]) < 1e-8) == 1){
      #       warning(paste0(
      #         "Alternating parameter estimates in iter =. ", iter, ". ", 
      #         "New estimates equal mean of preceding two estimates.\n"
      #       ))
      #       theta_new <- (theta_iters[5, ] + theta_iters[4, ]) / 2
      #     }
      #   }
      # }
    } # end if(!theta_fix){
    
    ### Update estimates
    if(method == "ML"){
      # if theta_new changed, adapt beta
      if(! sum(abs(theta_new1 - theta_new)) < eps){
        F_beta <- f_F_beta_ML(y, x_bdiag, x_bdiag_obs_list, beta = beta_new, 
                              theta = theta_new, order_ab, V_ed)
        S_beta <- f_S_beta_ML(y, x_bdiag, x_bdiag_obs_list, beta = beta_new, 
                              theta = theta_new, order_ab, V_ed)
        beta_new <- beta_old + as.vector(tryCatch(
          expr = solve(a = F_beta, b = S_beta), 
          error = function(e){ ginv(F_beta %*% S_beta) })
        ) 
      }
    }
    
    ### print results of current iteration
    if(verbose){
      if(method == "REML"){ 
        ll <- f_REML_ll(Y_obs, x_bdiag_obs, theta = theta_new, V_ed, I_mis_vec) 
        V_obs <- f_V(theta = theta_new, V_ed, D, order_ab)[I_obs_vec, I_obs_vec]
        V_inv_obs <- tryCatch(expr = solve(V_obs), 
                              error = function(e){ ginv(V_obs) })
        V_inv_x_obs <- V_inv_obs %*% x_bdiag_obs
        tmp_v1 <- t(x_bdiag_obs) %*% V_inv_x_obs
        tmp_v2 <- t(V_inv_x_obs) %*% Y_obs
        beta_new <- as.vector(tryCatch(
          expr = solve(a = tmp_v1, b = tmp_v2), 
          error = function(e){ ginv(tmp_v1 %*% tmp_v2) }))
      }
      if(method == "ML"){
        ll <- f_ML_ll(y, x_bdiag, beta = beta_new, theta = theta_new, V_ed, 
                      order_ab, O_d, pos_y_obs_d)
      }
      f_FS_print(beta_new, theta_new, ll, lambda, 
                 dig_theta, dig_beta, dig_lambda, dig_ll, dig_iter, 
                 wid_theta, wid_beta, wid_lambda, wid_ll, wid_iter)
    }
    
  } # end while
  #rm(theta_old, beta_old)
  #gc(verbose = FALSE)
  
  
  # Results parameter estimation ----
  # ____________________________________________________________________________
  
  #if(MMFH.verbose.steps){ print("Collect results point estimation") }
  if(round(lambda, digit = 4) == 0){
    warning("Fisher scoring step size was reduced towards 0 as elsewise non-feasible", 
            " estimates of random effect variances or covariances", 
            " would have occured.", 
            " Hence, the algorithm did not converge.")
    converged <- FALSE
  }
  beta <- beta_new
  names(beta) <- beta_names
  theta <- theta_new
  names(theta) <- theta_names
  Xbeta <- matrix(x_bdiag %*% beta, ncol = m, byrow = TRUE)
  colnames(Xbeta) <- colnames(y)
  rownames(Xbeta) <- rownames(y)
  
  V_ud <- f_V_ud(theta, order_ab)
  V_ud_inv <- tryCatch(expr = solve(V_ud), 
                       error = function(e){ ginv(V_ud) })
  V_ud_l <- f_V_ud_l(theta, order_ab)
  V_ud_l_mDmD <- f_V_ud_l_mDmD(theta, D, order_ab)
  V_ud_l_mDmD_obs <- lapply(V_ud_l_mDmD, function(V_ud_l_mDmD_theta_a){ 
    V_ud_l_mDmD_theta_a[I_obs_vec, I_obs_vec] 
  })
  V_obs <- f_V(theta, V_ed, D, order_ab)[I_obs_vec, I_obs_vec]
  V_inv_obs <- tryCatch(expr = solve(V_obs), 
                        error = function(e){ ginv(V_obs) })
  V_inv_x_obs <- crossprod(V_inv_obs, x_bdiag_obs)
  Q_obs <- f_Q(x_bdiag_obs, V_inv_x_obs)
  #if(is.null(colnames(Q_obs))){ rownames(Q_obs) <- colnames(Q_obs) <- colnames(x_bdiag_obs) }
  # if system = singular for calculation of Q ginv() is used instead of solve() and Q has no rownames
  P_obs <- f_P(x_bdiag_obs, Q_obs, V_inv_obs, V_inv_x_obs)
  
  # # for cAIC, see Han(2013)
  # P_obs.star <- diag(nrow(x_bdiag_obs)) - x_bdiag_obs %*% Q_obs %*% t(x_bdiag_obs) %*% V_inv_obs
  # S_tmp <- f_V(rep(0, theta_l), V_ed, D, order_ab)[I_obs_vec, I_obs_vec]
  # phi0.approx <- length(I_obs_vec) - sum(diag(S_tmp %*% V_inv_obs %*% P_obs.star))
  # r.s <- V_inv_obs %*% P_obs.star %*% Y_obs
  
  if(method == "REML"){
    F_REML <- f_F_REML(P_obs, V_ud_l_mDmD_obs, theta_l)
    F_inv <- tryCatch(expr = solve(F_REML), 
                      error = function(e){ ginv(F_REML) })
    ll <- f_REML_ll(Y_obs, x_bdiag_obs, theta, V_ed, I_mis_vec)
    
    # for cAIC, see Han(2013)
    # delta_s_tau2 <- sum(diag((V_inv_obs %*% P_obs.star)^2)) - 2 * t(r.s) %*% V_inv_obs %*% P_obs.star %*% r.s
    # delta_s_tau2_inv <- tryCatch(expr = solve(delta_s_tau2), 
    # error = function(e){ ginv(delta_s_tau2) })
    # phi0 <- phi0.approx - 2 * delta_s_tau2_inv %*% t(r.s) %*% V_inv_obs %*% P_obs.star %*% S_tmp %*% V_inv_obs %*% P_obs.star %*% r.s
  }
  
  if(method == "ML"){
    V <- f_V(theta, V_ed, D, order_ab)
    V_diag <- diag(V)
    # resids <- y - Xbeta
    # resid <- as.vector(t(resids))
    #rownames(resids) <- 1:D
    F_beta_ML <- f_F_beta_ML(y, x_bdiag, x_bdiag_obs_list, beta, theta, 
                             order_ab, O_d, V_ed)
    F_theta_ML <- f_F_theta_ML(y, x_bdiag, beta, theta, 
                               order_ab, O_d, V_ed)
    # F_inv_beta_ML <- tryCatch(expr = solve(F_beta_ML), 
    #                          error = function(e){ ginv(F_beta_ML) })
    F_inv_theta_ML <- tryCatch(expr = solve(F_theta_ML), 
                               error = function(e){ ginv(F_theta_ML) })
    # F_inv_tmp1 <- cbind(F_beta_ML, matrix(0, ncol = ncol(F_theta_ML), nrow = nrow(F_beta_ML)))
    # F_inv_tmp2 <- cbind(matrix(0, ncol = ncol(F_beta_ML), nrow = nrow(F_theta_ML)), F_theta_ML)
    # F_inv_tmp <- rbind(F_inv_tmp1, F_inv_tmp2)
    # F_inv <- tryCatch(expr = solve(F_inv_tmp), 
    # error = function(e){ ginv(F_inv_tmp) })
    # rm(F_inv_tmp1, F_inv_tmp2, F_inv_tmp)
    ll <- f_ML_ll(y, x_bdiag, beta, theta, V_ed, order_ab, O_d, pos_y_obs_d)
    
    # for cAIC, see Han(2013)
    # delta_s_tau2 <- sum(diag(V_inv_obs^(-2))) - 2 * t(r.s) %*% V_inv_obs %*% P_obs.star %*% r.s
    # delta_s_tau2_inv <- tryCatch(expr = solve(delta_s_tau2), 
    # error = function(e){ ginv(delta_s_tau2) })
    # phi0 <- phi0.approx - 2 * delta_s_tau2_inv %*% t(r.s) %*% V_inv_obs %*% P_obs.star %*% S_tmp %*% V_inv_obs %*% P_obs.star %*% r.s
  }
  
  # AIC <-(-2) * ll + 2 *(p) # speziell BIC Formel für bivariate ??
  # BIC <-(-2) * ll +(p) * log(length(Y_obs)) # speziell BIC Formel für bivariate ??
  std_err_beta <- sqrt(diag(Q_obs))
  tvalue <- beta / std_err_beta
  pvalue <- 2 * pnorm(abs(tvalue), lower.tail = FALSE)
  coef <- data.frame(beta, std_error = std_err_beta, tvalue, pvalue)
  y_num <- y # replace missing values by 0
  y_num[is.na(y_num)] <- 0
  resid <- y_num - Xbeta
  rm(theta_new, beta_new)
  
  
  # Preparation of MSE estimation ----
  # ____________________________________________________________________________
  
  #if(MMFH.verbose.steps){ print("Compute values") }
  
  # inverse of V_ed for observed variables, then 
  # V_ed in dimension m x m, where places of missing variables are filled by 0
  V_ed_0_inv <- lapply(1:D, function(d){
    I_obs_mat_d <- I_obs_mat[d, ]
    O_d_d <- O_d[[d]]
    V_ed_d_obs <- V_ed[[d]][I_obs_mat_d, I_obs_mat_d]
    if(!is.matrix(O_d_d)){
      # if only one variable observed
      t(matrix(O_d_d, nrow = 1)) %*% (V_ed_d_obs)^(-1) %*% O_d_d
    } else {
      # if at least two variables observed
      V_ed_d_obs_inv <- tryCatch(expr = solve(V_ed_d_obs), 
                                 error = function(e){ ginv(V_ed_d_obs) })
      t(O_d_d) %*% V_ed_d_obs_inv %*% O_d_d
    }
  })
  #names(V_ed_0_inv) <- 1:D
  
  
  # V_ed.0 = V_ed
  # V_ed, but where variables are missing, corresponding entries are 0
  # is equal to V_ed
  
  
  # Phi
  Phi <- lapply(1:D, function(d){
    tryCatch(expr = solve(V_ed_0_inv[[d]] + V_ud_inv), 
             error = function(e){ ginv(V_ed_0_inv[[d]] + V_ud_inv) })
  })
  #names(Phi) <- 1:D
  
  
  # random effects
  ref <- t(sapply(1:D, function(d){
    as.vector(Phi[[d]] %*% V_ed_0_inv[[d]] %*% resid[d, ])
  }))
  colnames(ref) <- colnames(y)
  rownames(ref) <- rownames(y)
  
  
  # ebp
  ebp <- Xbeta + ref
  
  
  # cAIC crude approximation from p. 66 Han(2013)
  # -2 log g_c(Y|theta_hat) + 2 phi_0
  # cAIC.approx <- as.vector(-2 * log(Reduce("+", 
  # mapply(function(res, V){
  # res <- res[[1]]
  # res[is.na(res)] <- 0
  # res <- as.vector(res)
  # V <- as.matrix(V)
  # t(res) %*% V %*% res
  # }, 
  # apply(y - ebp, 1, function(la) list(la)), # res
  # V_ed) # V
  #))) + 2 * phi0.approx
  
  
  # cAIC Theorem 2 in Han(2013)
  # -2 log g_c(Y|theta_hat) + 2 phi_0
  # cAIC <- as.vector(-2 * log(Reduce("+", 
  # mapply(function(res, V){
  # res <- res[[1]]
  # res[is.na(res)] <- 0
  # res <- as.vector(res)
  # V <- as.matrix(V)
  # t(res) %*% V %*% res
  # }, 
  # apply(y - ebp, 1, function(la) list(la)), # res
  # V_ed) # V
  #)) + 2 * phi0)
  
  
  # goodness of fit measures
  goodness <- c(ll = ll#, 
                # AIC = AIC, 
                # BIC = BIC#, 
                #cAIC = cAIC, 
                #cAIC.approx = cAIC.approx
  )
  
  
  # Derivatives of Phi w.r.t. theta_a 
  Phi_d_a <- lapply(1:D, function(d){
    Phi_d <- Phi[[d]]
    res_tmp <- lapply(1:theta_l, function(theta_a){
      Phi_d %*% V_ud_inv %*% V_ud_l[[theta_a]] %*% V_ud_inv %*% Phi[[d]]
    })
    names(res_tmp) <- theta_names
    res_tmp
  })
  #names(Phi_d_a) <- 1:D
  
  # Phi_d_a.11 <- t(sapply(Phi_d_a, function(d){ sapply(d, function(l){ l[1, 1] }) }))
  # Phi_d_a.22 <- t(sapply(Phi_d_a, function(d){ sapply(d, function(l){ l[2, 2] }) }))
  # Phi_d_a.12 <- t(sapply(Phi_d_a, function(d){ sapply(d, function(l){ l[1, 2] }) }))
  
  
  if(method == "ML"){ var_theta <- F_inv_theta_ML }
  if(method == "REML"){ var_theta <- F_inv }
  
  
  # MSE estimation ----
  # ____________________________________________________________________________
  
  #if(MMFH.verbose.steps){ print("mse estimation") }
  
  mse <- t(sapply(1:D, function(d){ 
    
    x_d <- lapply(x, function(x_k){ x_k[d, ] })
    Phi_d <- Phi[[d]]
    V_ed_0_inv_d <- V_ed_0_inv[[d]]
    V_ed_d <- V_ed[[d]]
    pos_y_obs_d_d <- pos_y_obs_d[[d]]
    Phi_d_a_d <- Phi_d_a[[d]]
    resid_d <- resid[d, ]
    
    ###
    ### Gd1
    ###
    
    Gd1 <- Phi_d %*% V_ed_0_inv_d %*% (V_ud + V_ed_d) %*% 
      V_ed_0_inv_d %*% Phi_d +
      V_ud - 2 * Phi_d %*% V_ed_0_inv_d %*% V_ud
    
    
    ###
    ### Gd2
    ### 
    
    # derivatives of h_d w.r.t. beta
    h_d_beta <- lapply(1:m, function(k){
      
      X_0_d_k <- matrix(0, nrow = m, ncol = p_k[k])
      X_0_d_k[k, ] <- x_d[[k]]
      
      if(k %in% pos_y_obs_d_d){
        sapply(1:p_k[k], function(beta_k_t){
          X_0_d_k[, beta_k_t] - Phi_d %*% V_ed_0_inv_d %*% X_0_d_k[, beta_k_t]
        })
      } else {
        X_0_d_k
      }
    })
    #names(h_d_beta) <- 1:m
    # list elements: k = 1, ..., m
    # each list element: 
    # dimension m x p_k 
    # rows correspond to h_d_beta_k_l in paper
    
    
    # trace(matrices h_d_beta_beta %*% covariance of corresponding beta)
    mat_tr_h_d_beta_cov_beta <- mapply(function(k, l){
      
      h_d_beta_k <- h_d_beta[[k]]
      h_d_beta_l <- h_d_beta[[l]]
      Q_obs_k_l <- Q_obs[beta_assign_k_mat[, k], beta_assign_k_mat[, l]]
      
      tmp_v1 <- mapply(function(var_a, var_b){
        # dimension h_d_beta_l_k_a_b: p_l x p_k
        h_d_beta_l_k_a_b <- h_d_beta_l[var_a, ] %*% t(h_d_beta_k[var_b, ])
        c_tmp <- h_d_beta_l_k_a_b %*% Q_obs_k_l
        sum(diag(c_tmp))
      }, order_a, order_b)
      # rows correspond to var_b
      # cols correspond to var_a
      matrix(tmp_v1, nrow = m, byrow = FALSE)
      
    }, order_a, order_b, SIMPLIFY = FALSE)
    # list elements correspond to G_d2, kl
    # each list:
    # var_a in columns
    # var_b in rows
    
    # names(mat_tr_h_d_beta_cov_beta) <- paste0(order_a, order_b)
    # all.equal(mat_tr_h_d_beta_cov_beta[["12"]], t(mat_tr_h_d_beta_cov_beta[["21"]]))
    # all.equal(mat_tr_h_d_beta_cov_beta[["13"]], t(mat_tr_h_d_beta_cov_beta[["31"]]))
    # all.equal(mat_tr_h_d_beta_cov_beta[["23"]], t(mat_tr_h_d_beta_cov_beta[["32"]]))
    
    
    # Gd2
    Gd2 <- Reduce('+', mat_tr_h_d_beta_cov_beta)
    
    
    ###
    ### Gd3
    ### 
    
    
    # derivatives of h_d w.r.t. theta
    h_d_theta <- sapply(1:theta_l, function(theta_a){
      Phi_d_a_d[[theta_a]] %*% V_ed_0_inv_d %*% resid_d
    })
    # rows correspond to h_d_theta_l in paper
    # columns correspond to elements of theta
    
    
    # trace(matrices h_d.heta_theta %*% variance of corresponding theta)
    tr_h_d_theta_var_theta <- mapply(function(var_a, var_b){
      tmp_v1 <-(h_d_theta[var_a, ] %*% t(h_d_theta[var_b, ])) %*% var_theta
      sum(diag(tmp_v1))
    }, order_a, order_b)
    mat_tr_h_d_theta_var_theta <- matrix(tr_h_d_theta_var_theta, 
                                         nrow = m, byrow = FALSE)
    # var_a in columns
    # var_b in rows
    
    
    # Gd3
    Gd3 <- mat_tr_h_d_theta_var_theta
    
    
    ### 
    ### MSE
    ### 
    
    
    # Gd1 + Gd2 + Gd3
    mse_d <- as.vector(Gd1 + Gd2 + 2 * Gd3)
    #mse1Gd3 <- Gd1 + Gd2 + 1 * Gd3
    names(mse_d) <- paste0(order_a, order_b)
    mse_d <- mse_d[order_a <= order_b]
    
    #mse1Gd3 <- as.vector(mse1Gd3)
    #names(mse1Gd3) <- paste0(order_a, order_b)
    #mse1Gd3 <- mse1Gd3[order_a <= order_b]
    
    # list("mse" = mse#, 
    #      #"mse1Gd3" = mse1Gd3
    #      )
    mse_d
  }))
  rownames(mse) <- rownames(y)
  
  
  # Output ----
  # ____________________________________________________________________________
  
  #if(MMFH.verbose.steps){ print("function output") }
  if(iter == maxiter){ 
    warning(paste0("Algorithm reached maximum number of ", 
                   maxiter, 
                   " iterations and hence did not converge.\n"))
    converged <- FALSE
  }
  
  #mse <- t(sapply(mse, function(m){ as.vector(m$mse) }))
  #mse1Gd3 <- t(sapply(mse, function(m){ as.vector(m$mse1Gd3) }))
  tmp_v1 <- paste0(order_a, order_b)
  tmp_v1[order_a == order_b] <- paste0("v", 1:m)
  tmp_v1[!order_a == order_b] <- paste0("cov", tmp_v1[!order_a == order_b])
  tmp_v1 <- tmp_v1[order_a <= order_b]
  colnames(mse) <- tmp_v1
  
  list(est = list(ebp = ebp, 
                  ref = ref, 
                  Xbeta = Xbeta, 
                  fit = list(method = method, 
                             covergence = converged, 
                             iterations = iter, 
                             estcoef = coef, 
                             refvar = theta, 
                             goodness = goodness)), 
       mse = mse)
}
