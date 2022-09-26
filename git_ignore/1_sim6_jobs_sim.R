# ______________________________________________________________________________
# ______________________________________________________________________________
#
# Multivariate Fay-Herriot estimator under partially missing values (MMFH) ----
#
# Code by Anna-Lena Wölwer in Mar 2022
# Dissertation
#
# ______________________________________________________________________________
# ______________________________________________________________________________

# ____________________________________________________________________________
# Set simulation ----
# ____________________________________________________________________________

# Aim: 
# EBLUPS, synthetic, MSE estimation
# Parameter estimation: REML
# D in 100, 200, 400 
# changing: rho_12



# ____________________________________________________________________________
# working directory ----
# ____________________________________________________________________________

setwd("C:/Users/alwoe/Seafile/me/Uni/Proj/Dis/github")

# ____________________________________________________________________________
# Load packages ----
# ____________________________________________________________________________


library(MASS)
library(Matrix)


# ____________________________________________________________________________
# simulation parameters ----
# ____________________________________________________________________________


  # settings for simulation
  verbose               <- FALSE
  eps                   <- 1e-5
  maxiter               <- 100


# ____________________________________________________________________________
# generate simulation data ----
# ____________________________________________________________________________

cat("\n", "generate simulation data", "\n")

  



###
### auxiliary data
###

### beta
# variable 1
beta.true.v1       <- c(2, 3)
p1                 <- length(beta.true.v1)
beta.true.v2       <- c(2, 3)
p2                 <- length(beta.true.v2)
beta.true.v3       <- c(2, 3)
p3                 <- length(beta.true.v3)
beta.true          <- c(beta.true.v1, beta.true.v2, beta.true.v3)
p                  <- p1 + p2 + p3

# beta positions for univariate FH
beta.pos.for.uv      <- list()
beta.pos.for.uv[[1]] <- 1 : (p1)
beta.pos.for.uv[[2]] <- (p1 + 1) : (p1 + p2)
beta.pos.for.uv[[3]] <- (p1 + p2 + 1) : p

### auxiliary data
range.x1           <- c(10, 100)
range.x2           <- c(10, 100)
range.x3           <- c(10, 100)

### X 
if (FALSE) {
  
  set.seed(1)
  d_tmp = 10
  
  for (d_tmp in c(10, 100, nr.areas.scen)) {
    
    # same X matrix for all scenarios and simulation runs
    
    # auxiliary variable 1
    x1 <- sapply(1:d_tmp, 
                 function(d) { runif(n   = 1, 
                                     min = range.x1[1], 
                                     max = range.x1[2]) }, simplify = "array")
    
    # auxiliary variable 2
    x2 <- sapply(1:d_tmp, 
                 function(d) { runif(n   = 1, 
                                     min = range.x2[1], 
                                     max = range.x2[2]) }, simplify = "array")
    
    # auxiliary variable 3
    x3 <- sapply(1:d_tmp, 
                 function(d) { runif(n   = 1, 
                                     min = range.x3[1], 
                                     max = range.x3[2]) }, simplify = "array")
    
    
    # same X matrix for all scenarios
    x      <- list()
    x[[1]] <- cbind("x11" = 1, "x12" = x1)
    x[[2]] <- cbind("x21" = 1, "x22" = x2)
    x[[3]] <- cbind("x31" = 1, "x32" = x3)
    
    saveRDS(x,
            file = paste0(path.user, "Dat/", name.sim, "/", "x123", "_D", d_tmp, ".RData"),
            version = 2)
  }
  rm(x1)
  rm(x2)
  rm(x3)
  rm(x)
}
x <- (readRDS(file = paste0(path.user, "Dat/", name.sim, "/", "x123", "_D", D, ".RData")))


###
### Xbeta
###

Xbeta.true           <- cbind( x[[1]] %*% beta.true[(1)           : (p1)   ], 
                               x[[2]] %*% beta.true[(1 + p1)      : (p1 + p2)], 
                               x[[3]] %*% beta.true[(1 + p1 + p2) : (p1 + p2 + p3)] ) 
dimnames(Xbeta.true) <- list("d" = 1:D, "var" = paste0("v", 1:m))
print(colMeans(Xbeta.true))
print(apply(Xbeta.true, 2, range))

p.vars               <- unlist(lapply(x, function(x_d) dim(x_d)[2]))
names(p.vars)        <- 1:m
p                    <- sum(p.vars)
a.tmp                <- unlist(lapply(1:m, function(var_k) { rep(var_k, times = p.vars[var_k]) }))
beta.names           <- list("d.v" = paste0("d", rep(1:D, each = m), "", "v", rep(1:m, times = D)),
                             "v.x" = paste0("v", a.tmp, "", "", 
                                            unlist(lapply(x, 
                                                          function(x_d) {
                                                            colnames(x_d)
                                                          } ))))[[2]]

X.bdiag               <- lapply (1:D, 
                                 function (d) { 
                                   as.matrix(do.call(Matrix::bdiag, 
                                                     lapply (x, 
                                                             function (x_d) { 
                                                               matrix(x_d[d,], nrow = 1) 
                                                             }) )) 
                                 })
X.bdiag               <- do.call("rbind", X.bdiag)
dimnames(X.bdiag)     <- list("d.v" = paste0("d", rep(1:D, each = m), ".", "v", rep(1:m, times = D)),
                              "v.x" = paste0(unlist(lapply (x, 
                                                            function (x_d) { 
                                                              colnames(x_d) 
                                                            } )) ))


# ____________________________________________________________________________
# example values for code debugging ----
# ____________________________________________________________________________

# data generation
s2u.start <- NULL
s2u       <- NULL

r = R.num_i[1] # simulation iteration r=1,...,R
d         <- 1 # domain d, d=1,...,D
method = "REML"

var_k     <- 1 # index variable k=1,...,m
var_l     <- 2 # index variable k=1,...,m
var_a     <- 2 # index variable k=1,...,m
var_b     <- 3 # index variable k=1,...,m
theta_a   <- 2 # index theta a=1,...,length(theta)
theta_b   <- 3 # index theta a=1,...,length(theta)
beta_k_t  <- 2 # index of p_k auxiliary variable for variable k 


# ____________________________________________________________________________
# Preparation simulation ----
# ____________________________________________________________________________

iter       <- 0
ll         <- 0
lambda     <- 0

dig.s2u    <- 3
dig.beta   <- 3
dig.lambda <- 3
dig.ll     <- 3
dig.iter   <- 0

wid.lambda <- 6
wid.iter   <- 5

s = 1

s2u.var.ij           <- cbind(rep(1:m, each  = m), rep(1:m, times = m))
s2u.var.ij           <- s2u.var.ij[s2u.var.ij[,1] <=  s2u.var.ij[,2],]
s2u.var.ij           <- rbind(s2u.var.ij[s2u.var.ij[,1] ==  s2u.var.ij[,2],],
                              s2u.var.ij[s2u.var.ij[,1] !=  s2u.var.ij[,2],])
colnames(s2u.var.ij) <- c("var_i","var_j")
s2u.var.ij           <- cbind(s2u.var.ij, "s2u" = 1:dim(s2u.var.ij)[1])
s2u.var.ij           <- data.frame(s2u.var.ij)
s2u.length           <- m * (m + 1) / 2
s2u.names            <- c(paste0("s2u", 1:m), 
                          paste0("rho", apply(s2u.var.ij[(m+1):s2u.length,], 1, 
                                              function (theta_a) { 
                                                paste0(theta_a[1], theta_a[2]) 
                                              }) ) )

r.not.conv <- vector()
r.not.calc <- vector()
r.iter.na  <- vector()
r.iter.1   <- vector()
r.MMFH.err <- vector()

var_k.order        <- var_a.order <- rep(1:m, each = m)
var_l.order        <- var_b.order <- rep(1:m, times = m)

mse_names                              <- paste0(var_a.order, var_b.order)
mse_names[ var_a.order == var_b.order] <- paste0("v", 1:m)
mse_names[!var_a.order == var_b.order] <- paste0("cov", mse_names[!var_a.order == var_b.order])
mse_names                              <- mse_names[var_a.order <= var_b.order]


# ____________________________________________________________________________
# Result arrays ----
# ____________________________________________________________________________

###
### results parameter estimates
###

FH_variant <- c("UFH", "MFH", "MMFH")
method.scen = c("REML")

# estimated parameters beta 
RES_param_beta <- array( dim = list( length(R.num.char_i),
                                     length(method.scen),
                                     length(FH_variant),
                                     (length(beta.true)) ),
                         dimnames = list( "R"     = R.num.char_i,
                                          "method" = method.scen,
                                          "FH"    = FH_variant,
                                          "param" = c(beta.names))   )

# estimated parameters theta 
RES_param_theta <- array( dim = list( length(R.num.char_i),
                                      length(method.scen),
                                      length(FH_variant),
                                      (length(s2u.true)) ),
                          dimnames = list( "R"     = R.num.char_i,
                                           "method" = method.scen,
                                           "FH"    = FH_variant,
                                           "param" = c(s2u.names))   )

# estimated parameters theta 
RES_param_theta <- array( dim = list( length(R.num.char_i),
                                      length(method.scen),
                                      length(FH_variant),
                                      (length(s2u.true)) ),
                          dimnames = list( "R"     = R.num.char_i,
                                           "method" = method.scen,
                                           "FH"    = FH_variant,
                                           "param" = c(s2u.names))   )

FH_variant_plus_true <- c("true", "UFH", "MFH", "MMFH")

# eblups
RES_EBLUPS <- array( dim = list( length(R.num.char_i),
                                 length(method.scen),
                                 length(FH_variant_plus_true),
                                 m,
                                 D ),
                     dimnames = list( "R"     = R.num.char_i,
                                      "method" = method.scen,
                                      "FH"    = FH_variant_plus_true,
                                      "Var" = 1:m,
                                      "d" = 1:D)   )

# synthetic 
RES_Syns <- array( dim = list( length(R.num.char_i),
                               length(method.scen),
                               length(FH_variant_plus_true),
                               m,
                               D ),
                   dimnames = list( "R"     = R.num.char_i,
                                    "method" = method.scen,
                                    "FH"    = FH_variant_plus_true,
                                    "Var" = 1:m,
                                    "d" = 1:D)   )

# mse eblup 1 
RES_EBLUP_MSE <- array( dim = list( length(R.num.char_i),
                                    length(method.scen),
                                    length(FH_variant),
                                    m*(m+1)/2,
                                    D ),
                        dimnames = list( "R"     = R.num.char_i,
                                         "method" = method.scen,
                                         "FH"    = FH_variant,
                                         "cov" = mse_names,
                                         "d" = 1:D)   )

# mse eblup 2 
RES_EBLUP_MSE_alt <- array( dim = list( length(R.num.char_i),
                                        length(method.scen),
                                        length(FH_variant),
                                        m*(m+1)/2,
                                        D ),
                            dimnames = list( "R"     = R.num.char_i,
                                             "method" = method.scen,
                                             "FH"    = FH_variant,
                                             "cov" = mse_names,
                                             "d" = 1:D)   )


# mse syn
RES_SYN_MSE <- array( dim = list( length(R.num.char_i),
                                  length(method.scen),
                                  length(FH_variant),
                                  m*(m+1)/2,
                                  D ),
                      dimnames = list( "R"     = R.num.char_i,
                                       "method" = method.scen,
                                       "FH"    = FH_variant,
                                       "cov" = mse_names,
                                       "d" = 1:D)   )


# ____________________________________________________________________________
#
# Run simulation ----
#
# ____________________________________________________________________________

cat("Run simulation, job =", job, ", ", name.data.gen, "\n")
start.time <- Sys.time()
set.seed(job)

print(range(R.num_i))
R_max <- max(R.num_i)
R_l <- length(R.num_i)

 # r = R.num_i[1]
  
  
  # ____________________________________________________________________________
  # Generate input data ----
  # ____________________________________________________________________________
  
  # random effects
  u                      <- sapply(1:D, 
                                   function(d) { 
                                     rmvnorm_copy(n = 1, mean = rep(0, m), sigma = V.ud.true)
                                   }, 
                                   simplify = "array")[1,,]
  u                      <- aperm(u, c(2,1))
  dimnames(u)            <- list("d" = as.character(1:D), "var" = paste0("v", 1:m))
  
  # sampling errors
  e                      <- sapply(1:D, 
                                   function(d) { 
                                     rmvnorm_copy(n = 1, mean = rep(0, m), sigma = V.ed[[d]])
                                   }, 
                                   simplify = "array")[1,,]
  e                      <- aperm(e, c(2,1))
  dimnames(e)            <- list("d" = as.character(1:D), "var" = paste0("v", 1:m))
  
  # yhat, y true and y observed
  y.true               <- Xbeta.true + u
  y.obs.r              <- y.true + e
  
  
  # ____________________________________________________________________________
  # missing values ----
  # ____________________________________________________________________________
  
  y              <- y.obs.r
  y[(1:(D/2)),1] <- NA
  
  for (var_k in 1:m) { # var_k = 1
    s2y[which(is.na(y[, paste0("v", var_k)]))] <- sapply(which(is.na(y[, paste0("v", var_k)])), function(i) { 
      s2y[[i]][var_k,     ] <- NA
      s2y[[i]][     ,var_k] <- NA
      s2y[[i]] 
    }, simplify = FALSE)
  }
  
  
  
  # ____________________________________________________________________________
  # verbose ----
  # ____________________________________________________________________________
  
  if (!test_code) {
    
    rep_tmp <- 20
    
    if (r %in% 1:rep_tmp | (r/rep_tmp) == round(r/rep_tmp)) {
      verbose = TRUE
    } else {
      verbose = FALSE
    }
    
    if (r %in% 1:rep_tmp | (r/rep_tmp) == round(r/rep_tmp)) {
      if (verbose) {
        time.between <- Sys.time()
        a.tmp        <- difftime(time.between, start.time, units = "hours") 
        name.unit    <- units(a.tmp)
        b.tmp        <- as.numeric(a.tmp)
        r_tmp <- which(R.num_i == r)
        cat(paste0("r=", r, "/", R_max, "; ",
                   "theta=",
                   paste0(as.vector(s2u.true), collapse = ","), ", ",
                   "beta=", 
                   paste0(beta.true, collapse = ","), ", ",
                   "job ",              job, ", ", 
                   "rho12", "_", cor.ref12.true, ", ",
                   "D", D, ", ",
                   "R_", R_i,
                   ", past=", round(b.tmp, digits = 4), " ", name.unit, 
                   ", remain=", round((R_l-r_tmp) * (b.tmp / r_tmp), digits = 4), " ", name.unit, "\n"))
      }
      if (!verbose) {
        cat(" r =", r, "/", R_max, "\n")
      }
    } else {
      cat(" r =", r, "/", R_max, "\n")
    }
  } else {
    # show all results
    verbose = TRUE
    cat("\n", " r =", r, "/", R_max, ";", 
        "true beta =", 
        paste0(beta.true, collapse = ", "), "|", 
        "true theta =",
        paste0(as.vector(s2u.true), collapse = ", "),
        "\n")
  }
  
  
  # ____________________________________________________________________________
  # ML and REML ----
  # ____________________________________________________________________________
  
  
  for (method in method.scen) { # method = "REML"
    
    RES_EBLUPS[as.character(r), method, "true", , ] <- t(y.true)
    RES_Syns[as.character(r), method, "true", , ]   <- t(y.true)
    
    
    # ____________________________________________________________________________
    # FH ----
    # ____________________________________________________________________________
    
    if (verbose) { cat(paste0("\n", "UFH ", method, ":", "\n")) }
    ls.bef.model <- ls()
    
    for (var_k in 1:m) { # var_k = 1
      
      obs.tmp  <- which(!is.na(y[, paste0("v",var_k)]))
      mis.tmp  <- which( is.na(y[, paste0("v",var_k)]))
      
      if (length(obs.tmp) > 1) {
        
        RES_UFH <- mseFH(formula = y[obs.tmp, var_k] ~ - 1 + x[[var_k]][obs.tmp,], 
                         vardir  = c(unlist(lapply(obs.tmp, function(d) { s2y[[d]][var_k, var_k] } ))), 
                         method  = method,
                         MAXITER = maxiter,
                         PRECISION = eps)
        
        if (!class(RES_UFH) == "try-error" ) {
          
          RES_param_theta[as.character(r), method, "UFH", var_k]                         <- RES_UFH$est$fit$refvar
          RES_param_beta[as.character(r), method, "UFH", beta.pos.for.uv[[var_k]]]       <- RES_UFH$est$fit$estcoef$beta
          
          RES_EBLUPS[as.character(r), method, "UFH", as.character(var_k), obs.tmp]       <- RES_UFH$est$eblup
          RES_EBLUP_MSE[as.character(r), method, "UFH", paste0("v", var_k), obs.tmp]     <- RES_UFH$mse
          RES_EBLUP_MSE_alt[as.character(r), method, "UFH", paste0("v", var_k), obs.tmp] <- RES_UFH$mse
          
          RES_Syns[as.character(r), method, "UFH", as.character(var_k), mis.tmp]         <- x[[var_k]][mis.tmp,] %*% RES_UFH$est$fit$estcoef$beta
          
          ref.v.tmp                 <- RES_UFH$est$fit$refvar
          V_hat.tmp                 <- diag(c(unlist(lapply(obs.tmp, function(d) { s2y[[d]][var_k, var_k] + ref.v.tmp } ))))
          V_hat_inv.tmp             <- tryCatch(expr  =              solve( V_hat.tmp ),
                                                error = function(e) { ginv( V_hat.tmp ) })
          X_V_hat_inv_X.tmp         <- t(x[[var_k]][obs.tmp,]) %*% V_hat_inv.tmp %*% x[[var_k]][obs.tmp,]
          X_V_hat_inv_X_inv.tmp     <- tryCatch(expr  =              solve( X_V_hat_inv_X.tmp ),
                                                error = function(e) { ginv( X_V_hat_inv_X.tmp ) })
          RES_SYN_MSE[as.character(r), method, "UFH", paste0("v", var_k), mis.tmp] <- c(unlist(lapply(mis.tmp, 
                                                                                                      function(d) { 
                                                                                                        tcrossprod( x[[var_k]][d,] %*% X_V_hat_inv_X_inv.tmp, x[[var_k]][d,] ) + 
                                                                                                          ref.v.tmp
                                                                                                      } )))
        } 
      }
      ls.tmp <- grep(".tmp", ls(), value = TRUE)
      rm(list = ls.tmp)
    }
    
    if (verbose) {
      cat(paste0("\n", "uv FH ", method, ":", "\n"))
      
      beta.new   <- RES_param_beta[as.character(r), method, "UFH", ]
      s2u.new    <- RES_param_theta[as.character(r), method, "UFH", ]     
      wid.s2u    <- range(nchar(round(s2u.new[!is.na(s2u.new)],  digits = dig.s2u)))[2] + 1;  if(wid.s2u < 5)  { wid.s2u <- 5 }
      wid.beta   <- range(nchar(round(beta.new[!is.na(beta.new)], digits = dig.beta)))[2] + 1; if(wid.beta < 5) { wid.beta <- 5 }
      wid.ll     <- range(nchar(round(ll,  digits = dig.ll)))[2]; if(wid.ll < 18)  { wid.ll <- 18 }
      
      verb1 <- c( formatC(c("iter"),                                           width = wid.iter),
                  formatC(c("loglike"),                                        width = wid.ll),
                  formatC(c("lambda"),                                         width = wid.lambda),
                  formatC(c(s2u.names),                                        width = wid.s2u),
                  formatC(c(paste0("b",beta.names)),                           width = wid.beta))
      verb2 <- c( formatC(c(iter),                                             width = wid.iter),
                  formatC(c(round(ll,                  digits = dig.ll)),      width = wid.ll),
                  formatC(c(as.character(round(lambda, digits = dig.lambda))), width = wid.lambda),
                  formatC(c(round(s2u.new,             digits = dig.s2u)),     width = wid.s2u),
                  formatC(c(round(beta.new,            digits = dig.beta)),    width = wid.beta))
      cat(verb1, "\n")
      cat(verb2, "\n")
      
    }          
    ls.aft.model <- ls()
    rm(list = ls.aft.model[!ls.aft.model %in% ls.bef.model])
    if (verbose) { cat(paste0("\n")) }
    
    
    # ____________________________________________________________________________
    # MFH ----
    # ____________________________________________________________________________
    
    if (verbose) { cat(paste0("\n", "MFH ", method, ":", "\n")) }
    ls.bef.model <- ls()
    
    s2u <- NULL
    s2u.start <- NULL
    
    obs.tmp  <- which(!is.na(rowMeans(y)))
    mis.tmp  <- which( is.na(rowMeans(y)))
    
    RES_MFH <- try(MMFH(y         = y[obs.tmp,],
                        x         = lapply(x, function (x_k) { x_k[obs.tmp,] }),
                        s2y       = s2y[obs.tmp],
                        s2u       = s2u,
                        s2u.start = s2u.start,
                        method    = method,
                        eps       = eps,
                        maxiter   = maxiter,
                        verbose   = verbose), silent = FALSE)
    
    
    if (!class(RES_MFH) == "try-error" ) {
      
      RES_param_theta[as.character(r), method, "MFH", ]          <- RES_MFH$est$fit$refvar
      RES_param_beta[as.character(r), method, "MFH", ]           <- RES_MFH$est$fit$estcoef$beta
      
      RES_EBLUPS[as.character(r), method, "MFH", , obs.tmp]        <- t(RES_MFH$est$ebp)
      RES_EBLUP_MSE[as.character(r), method, "MFH", , obs.tmp]     <- t(RES_MFH$mse)
      RES_EBLUP_MSE_alt[as.character(r), method, "MFH", , obs.tmp] <- t(RES_MFH$mse1Gd3)
      
      RES_Syns[as.character(r), method, "MFH", , mis.tmp]        <- t(sapply(1:m, function (var_k) {
        x[[var_k]][mis.tmp, ] %*% RES_param_beta[as.character(r), method, "MFH", beta.pos.for.uv[[var_k]]]
      }))
      
      ref.v.tmp                 <- RES_MFH$est$fit$refvar
      V_u_hat                   <- matrix(c(ref.v.tmp[1],
                                            rep(ref.v.tmp[4] * sqrt( ref.v.tmp[1] * ref.v.tmp[2] ), 1),
                                            rep(ref.v.tmp[5] * sqrt( ref.v.tmp[1] * ref.v.tmp[3] ), 1),
                                            #
                                            rep(ref.v.tmp[4] * sqrt( ref.v.tmp[1] * ref.v.tmp[2] ), 1),
                                            ref.v.tmp[2],
                                            rep(ref.v.tmp[6] * sqrt( ref.v.tmp[2] * ref.v.tmp[3] ), 1),
                                            #
                                            rep(ref.v.tmp[5] * sqrt( ref.v.tmp[1] * ref.v.tmp[3] ), 1),
                                            rep(ref.v.tmp[6] * sqrt( ref.v.tmp[2] * ref.v.tmp[3] ), 1),
                                            ref.v.tmp[3] ),
                                          byrow = TRUE,
                                          ncol = m)
      
      V_hat.tmp               <- lapply (obs.tmp, function (d) { s2y[[d]] + V_u_hat } )
      V_hat.tmp               <- Matrix::bdiag(V_hat.tmp)
      V_hat_inv.tmp             <- tryCatch(expr  =              solve( V_hat.tmp ),
                                            error = function(e) { ginv( V_hat.tmp ) })
      
      
      X.bdiag_obs               <- lapply (obs.tmp, 
                                           function (d) { 
                                             as.matrix(do.call(Matrix::bdiag, 
                                                               lapply (x, 
                                                                       function (x_d) { 
                                                                         matrix(x_d[d,], nrow = 1) 
                                                                       }) )) 
                                           })
      X.bdiag_obs               <- do.call("rbind", X.bdiag_obs)
      dimnames(X.bdiag_obs)     <- list("d.v" = paste0("d", rep(obs.tmp, each = m), ".", "v", rep(1:m, times = length(obs.tmp))),
                                        "v.x" = paste0(unlist(lapply (x, 
                                                                      function (x_d) { 
                                                                        colnames(x_d) 
                                                                      } )) ))
      X_V_hat_inv_X.tmp         <- t(X.bdiag_obs) %*% V_hat_inv.tmp %*% X.bdiag_obs
      rm(V_hat_inv.tmp)
      X_V_hat_inv_X_inv.tmp     <- tryCatch(expr  =              solve( X_V_hat_inv_X.tmp ),
                                            error = function(e) { ginv( X_V_hat_inv_X.tmp ) })
      RES_SYN_MSE[as.character(r), method, "MFH", , mis.tmp] <- t(do.call('rbind', lapply(mis.tmp, 
                                                                                          function(d) { 
                                                                                            x_d_tmp <- as.matrix(do.call(Matrix::bdiag, 
                                                                                                                         lapply (x, 
                                                                                                                                 function (x_d) { 
                                                                                                                                   matrix(x_d[d,], nrow = 1) 
                                                                                                                                 }) ))
                                                                                            res_tmp <- x_d_tmp %*% X_V_hat_inv_X_inv.tmp %*% t(x_d_tmp) + V_u_hat
                                                                                            res_tmp <- as.vector(res_tmp)
                                                                                            names(res_tmp) <- paste0(var_a.order, var_b.order)
                                                                                            res_tmp        <- res_tmp[var_a.order <= var_b.order]
                                                                                            res_tmp
                                                                                          } ) ))
    }
    
    
    # if (is.na(RES_MFH$est$fit$iterations) ) {
    #   r.iter.na <- c(r.iter.na, r)
    # } else {
    #   if (RES_MFH$est$fit$iterations == 1) {
    #     r.iter.1 <- c(r.iter.1, r)
    #   }
    # }
    # 
    # if (is.na(RES_MFH$est$fit$covergence)) {
    #   r.not.calc <- c(r.not.calc, r)
    # } else {
    #   if (RES_MFH$est$fit$covergence == FALSE) {
    #     r.not.conv <- c(r.not.conv, r)
    #   }
    # }
    
    ls.aft.model <- ls()
    rm(list = ls.aft.model[!ls.aft.model %in% ls.bef.model])
    if (verbose) { cat(paste0("\n")) }
    
    
    # ____________________________________________________________________________
    # MMFH ----
    # ____________________________________________________________________________
    
    if (verbose) { cat(paste0("\n", "MMFH ", method, ":", "\n")) }
    ls.bef.model <- ls()
    
    s2u <- NULL
    s2u.start <- NULL
    
    RES_MMFH <- try(MMFH(y         = y,
                         x         = x,
                         s2y       = s2y,
                         s2u       = s2u,
                         s2u.start = s2u.start,
                         method    = method,
                         eps       = eps,
                         maxiter   = maxiter,
                         verbose   = verbose), silent = FALSE)
    
    
    if (!class(RES_MMFH) == "try-error" ) {
      
      RES_param_theta[as.character(r), method, "MMFH", ]          <- RES_MMFH$est$fit$refvar
      RES_param_beta[as.character(r), method, "MMFH", ]           <- RES_MMFH$est$fit$estcoef$beta
      
      RES_EBLUPS[as.character(r), method, "MMFH", , ]        <- t(RES_MMFH$est$ebp)
      RES_EBLUP_MSE[as.character(r), method, "MMFH", , ]     <- t(RES_MMFH$mse)
      RES_EBLUP_MSE_alt[as.character(r), method, "MMFH", , ] <- t(RES_MMFH$mse1Gd3)
      
    }
    
    if (is.na(RES_MMFH$est$fit$iterations) ) {
      r.iter.na <- c(r.iter.na, r)
    } else {
      if (RES_MMFH$est$fit$iterations == 1) {
        r.iter.1 <- c(r.iter.1, r)
      }
    }
    
    if (is.na(RES_MMFH$est$fit$covergence)) {
      r.not.calc <- c(r.not.calc, r)
    } else {
      if (RES_MMFH$est$fit$covergence == FALSE) {
        r.not.conv <- c(r.not.conv, r)
      }
    }
    
    ls.aft.model <- ls()
    rm(list = ls.aft.model[!ls.aft.model %in% ls.bef.model])
    if (verbose) { cat(paste0("\n")) }
    
    
    # end simulation
  }


