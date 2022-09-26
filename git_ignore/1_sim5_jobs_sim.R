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
# Parameter estimation 
# ML and REML 
# D in 50, 100, 200, 300 
# everything else: fixed
warning("In diesem code sind die synthetischen Schätzer falsch, da müsste 2 mal x[[var_k]][obs.tmp,] durch x[[var_k]][mis.tmp,] ersetzt werden,
        aber diese Werte verwenden wir in sim5 nicht, nur in sim6, da wurde es korrigiert")


# ____________________________________________________________________________
# Set simulation ----
# ____________________________________________________________________________

start.code <- Sys.time()
name.sim   <- "sim5" # name of simulation
test_code  <- FALSE
job_tmp <- 67


# ____________________________________________________________________________
# working directory ----
# ____________________________________________________________________________

wd.tmp       <- getwd()
wd.tmp.start <- substr(wd.tmp, 1, 1)
if (wd.tmp.start == "C") { SERVER = FALSE } else { SERVER = TRUE }
path.end.sim <- "Uni/Proj/MultivarFHMis/Sim/"
path.begin.home <- "C:/Users/alwoe/Seafile/me/"
path.begin.server <- "/home/woelwer/"
if (SERVER == TRUE) {
  path.begin <- path.begin.server
  node.name  <- Sys.info()["nodename"]
  cat("\n\nnode.name =", node.name, "\n\n\n\n")
  test_code    <- FALSE
  #simulation   <- TRUE
  warning(node.name)
  print(sessionInfo())
  warning(sessionInfo())
}
if (SERVER == FALSE) {
  path.begin <- path.begin.home
}
path.user <- paste0(path.begin, path.end.sim)
cat("path.user =", path.user, "\n")
setwd(path.user)


# ____________________________________________________________________________
# Load packages ----
# ____________________________________________________________________________

list_packages <- c("MASS", "Matrix")
for (pack in list_packages) {
  cat(paste0("load  ", pack, "\n"))
  if (SERVER) {
    tryCatch(expr = {
      library(pack, character.only = TRUE, 
              lib.loc = paste0(path.user, "packages/to_load/"))
    },
    error = function(e) {
      message("Error loading from folder package ", pack)
      library(pack, character.only = TRUE)
    })
  } else {
    library(pack, character.only = TRUE)
  }
  if (SERVER) {
    cat(paste0("package version of ", pack, " = ", packageVersion(pack), "\n"))
    warning(paste0("package version of ", pack, " = ", packageVersion(pack), ""))
  }
}


# ____________________________________________________________________________
# simulation parameters ----
# ____________________________________________________________________________

# number of variables
m <- 3

# variance random effects
v.ref.true <- c(2, 2, 2)  

# correlation random effects
cor.ref13.true            <- c("p75") 
cor.ref13.true.num        <- c(.75) 
names(cor.ref13.true.num) <- cor.ref13.true

# correlation random effects
cor.ref23.true            <- c("p25") 
cor.ref23.true.num        <- c(.25) 
names(cor.ref23.true.num) <- cor.ref23.true

# variance sampling errors
v.rer.true <- c(3, 3, 3)  

# correlation sampling errors
cor.rer12.true            <- c("n5") 
cor.rer12.true.num        <- c(-0.5) 
names(cor.rer12.true.num) <- cor.rer12.true

# correlation sampling errors
cor.rer13.true            <- c("0") 
cor.rer13.true.num        <- c(0) 
names(cor.rer13.true.num) <- cor.rer13.true

# correlation sampling errors
cor.rer23.true            <- c("0") 
cor.rer23.true.num        <- c(0) 
names(cor.rer23.true.num) <- cor.rer23.true

if (!test_code) {
  # settings for simulation
  verbose               <- FALSE
  eps                   <- 1e-5
  maxiter               <- 100
} else {
  warning("test code")
  # settings for test runs
  verbose               <- TRUE
  eps                   <- 1e-1
  maxiter               <- 5
}

if (SERVER) { 
  print(sessionInfo())
  print(maxiter)
  print(eps)
}



# ____________________________________________________________________________
# simulation scenarios ----
# ____________________________________________________________________________

# universe correlation of random effects
# cor.ref12.true.scen            <- c("n75", "n5", "n25", "0", "p25", "p5", "p75") 
# cor.ref12.true.num.scen        <- c(-.75, -.5, -.25, 0, .25, .5, .75) 
# names(cor.ref12.true.num.scen) <- cor.ref12.true.scen

cor.ref12.true.scen            <- c("p5") 
cor.ref12.true.num.scen        <- c(.5) 
names(cor.ref12.true.num.scen) <- cor.ref12.true.scen


# number of areas
nr.areas.scen <- c(50, 100, 200, 300)

# R iterations
R                    <- 2050
R.subset.nr          <- 500
R.scen               <- 1:R.subset.nr
jobs.per.R.subset    <- round(R / R.subset.nr)
R.scen.num.list      <- list()
for (i in 1:R.subset.nr) {
  if (i == 1) { 
    beg_tmp <- 1
    end_tmp <- min(R, (jobs.per.R.subset+1))
  }
  R.scen.num.list[[i]] <- beg_tmp:end_tmp
  beg_tmp              <- end_tmp + 1
  end_tmp              <- min(beg_tmp + jobs.per.R.subset, R)
}
unlist(lapply(R.scen.num.list, length))


# possible combinations
sim_jobs <- data.frame(expand.grid(
  "rho12"       = cor.ref12.true.scen, 
  "D"           = nr.areas.scen,
  #"method"      = method.scen,
  "R"           = R.scen ))
sim_jobs$rho12  <- as.character(sim_jobs$rho12)
sim_jobs$D      <- as.character(sim_jobs$D)
#sim_jobs$method <- as.character(sim_jobs$method)
sim_jobs$job    <- 1:nrow(sim_jobs)

# delete some sim_jobs
# sim_jobs <- sim_jobs[sim_jobs$rho12 %in% "p5",]
# sim_jobs <- sim_jobs[!sim_jobs$rho12 %in% "n75",]
# sim_jobs <- sim_jobs[!sim_jobs$sig12 %in% "n75",]
# sim_jobs <- sim_jobs[!sim_jobs$rho12 %in% "n5",]
# sim_jobs <- sim_jobs[!sim_jobs$sig12 %in% "n5",]
dim(sim_jobs)
# 2000       

sim_jobs$R.num.range.char_i <- NA


# i = 1
for (i in unique(sim_jobs$R)) {
  
  a <- range(unlist(R.scen.num.list[[i]]))
  a2 <- paste0(a[1], "_", a[2])
  sim_jobs[sim_jobs$R == i, "R.num.range.char_i"] <- a2
  
}

a_tmp <- ""
sim_jobs$name.data.gen <- paste0("rho", "_", sim_jobs$rho12, "_", 
                                 "D", sim_jobs$D, "_",
                                 "R_", sim_jobs$R.num.range.char_i,
                                 a_tmp) 


# ____________________________________________________________________________
# generate submit file ----
# ____________________________________________________________________________

if (FALSE) {
  
  jobs_submit <- sim_jobs$job
  source(paste0(path.user, "Code/", name.sim, "/", "2_sim5_submit_gen_1", ".R"))
  
}

if (FALSE) {
  dir_to_create <- c( paste0(path.user, "Dat/", name.sim, "/"),
                      paste0(path.user, "Dat/", name.sim, "/sim_res_conv/"),
                      paste0(path.user, "Dat/", name.sim, "/sim_res/"),
                      paste0(path.user, "Dat/", name.sim, "/out/"),
                      paste0(path.user, "Dat/", name.sim, "/err/"),
                      paste0(path.user, "Dat/", name.sim, "/raw_sim_out/"),
                      paste0(path.user, "Dat/", name.sim, "/raw_sim_out/", "RES_param_beta", "/"),
                      paste0(path.user, "Dat/", name.sim, "/raw_sim_out/", "RES_param_theta", "/"),
                      paste0(path.user, "Dat/", name.sim, "/raw_sim_out/", "RES_EBLUPS", "/"),
                      paste0(path.user, "Dat/", name.sim, "/raw_sim_out/", "RES_EBLUP_MSE", "/"),
                      paste0(path.user, "Dat/", name.sim, "/raw_sim_out/", "RES_EBLUP_MSE_alt", "/"),
                      paste0(path.user, "Dat/", name.sim, "/raw_sim_out/", "RES_Syns", "/"),
                      paste0(path.user, "Dat/", name.sim, "/raw_sim_out/", "RES_SYN_MSE", "/"),
                      #
                      paste0(path.user, "Dat/", name.sim, "/raw_sim_out/", "r.iter.na", "/"),
                      paste0(path.user, "Dat/", name.sim, "/raw_sim_out/", "r.iter.1", "/"),
                      paste0(path.user, "Dat/", name.sim, "/raw_sim_out/", "r.not.conv", "/"),
                      paste0(path.user, "Dat/", name.sim, "/raw_sim_out/", "r.not.calc", "/")
  )
  for (d in dir_to_create) {
    if ( ! dir.exists(d) ) {
      dir.create(d, showWarnings = F)
    }
  }
}


# ____________________________________________________________________________
# current job information ----
# ____________________________________________________________________________

if (SERVER == TRUE) {
  job <- as.numeric(Sys.getenv("var")) # job number
  a <- round(runif(1, 1, 120))
  print(Sys.time())
  cat("\nrunif for", a, "seconds\n")
  Sys.sleep(a)
}
if (SERVER == FALSE) { 
  sim_jobs[sim_jobs$rho12 == "p5" & sim_jobs$sig12 == "p25" & sim_jobs$D == "100" & sim_jobs$R == 7,]
  job <- job_tmp
}
set.seed(job)

sim_job_i          <- sim_jobs[sim_jobs$job == job,]
#method             <- sim_job_i$method

cor.ref12.true     <- sim_job_i$rho12
cor.ref12.true.num <- cor.ref12.true.num.scen[cor.ref12.true]

D                  <- as.numeric(sim_job_i$D)
R_i                <- sim_job_i$R
R.num_i            <- R.scen.num.list[[R_i]]
if (test_code) {
  warning("test code")
  #R       <- 5
  R.num_i <- 9:11
  D       <- 100
}
D_char             <- as.character(D)
R.num.range_i      <- range(R.num_i)
R.num.range.char_i <- paste0(R.num.range_i[1], "_", R.num.range_i[2])
R.num.char_i       <- as.character(R.num_i)

cat(paste0("\n",
           "job ",               job,
           " | D = ",            D,
           #" | ",                method,
           " | rho = ",          cor.ref12.true,
           " | R = ",           R.num.range.char_i,
           "\n"))
if (SERVER) {
  warning(paste0("\n",
                 "job ",               job,
                 " | D = ",            D,
                # " | ",                method,
                 " | rho = ",          cor.ref12.true,
                 " | R = ",           R.num.range.char_i,
                 "\n"))
}
if (is.na(D))              { stop("is.na(D) = TRUE") }
#if (is.na(method))         { stop("is.na(method) = TRUE") }
if (is.na(cor.ref12.true)) { stop("is.na(cor.ref12.true) = TRUE") }
if (anyNA(R.num_i))        { stop("anyNA(R.num_i) = TRUE") }


# ____________________________________________________________________________
# names for output ----
# ____________________________________________________________________________

if (test_code) {
  warning("test code")
  a_tmp <- "_test"
} else {
  a_tmp <- ""
}
name.data.gen <- paste0("rho", "_", cor.ref12.true, "_", 
                        "D", D_char, "_",
                        "R_", R.num.range.char_i,
                        a_tmp) 

# file were results are later stored
paste0(path.user, "Dat/", name.sim, "/raw_sim_out/", "Res_param_beta", "/", name.data.gen, ".RData")


# ____________________________________________________________________________
# Load functions ----
# ____________________________________________________________________________

source(paste0(path.user, "/Code/", "MMFH_wt_packages", ".R"))

source(paste0(path.user, "/Code/", name.sim, "/", "0_eblupFH", ".R"))
eblupFH <- eblupFH_rebuilt; rm(eblupFH_rebuilt)
source(paste0(path.user, "/Code/", name.sim, "/", "0_mseFH", ".R"))
mseFH <- mseFH_rebuilt; rm(mseFH_rebuilt)

source(paste0(path.user, "/Code/", name.sim, "/", "0_rmvnorm_copy", ".R"))


# ____________________________________________________________________________
# generate simulation data ----
# ____________________________________________________________________________

cat("\n", "generate simulation data", "\n")

###
### variance-covariance matrices 
###

# variance-covariance matrix of sampling errors
s2y.true.v         <- v.rer.true
s2y.true.mat.cov   <- matrix(c(v.rer.true[1],
                               rep(cor.rer12.true.num * sqrt( v.rer.true[1] * v.rer.true[2] ), 1),
                               rep(cor.rer13.true.num * sqrt( v.rer.true[1] * v.rer.true[3] ), 1),
                               #
                               rep(cor.rer12.true.num * sqrt( v.rer.true[1] * v.rer.true[2] ), 1),
                               v.rer.true[2],
                               rep(cor.rer23.true.num * sqrt( v.rer.true[2] * v.rer.true[3] ), 1),
                               #
                               rep(cor.rer13.true.num * sqrt( v.rer.true[1] * v.rer.true[3] ), 1),
                               rep(cor.rer23.true.num * sqrt( v.rer.true[2] * v.rer.true[3] ), 1),
                               v.rer.true[3]
),
m, m, byrow = TRUE)
V.ed               <- lapply(1:D, function(a) { s2y.true.mat.cov })
s2y                <- lapply(1:D, function(d) { V.ed[[d]] })
V.ed_inv           <- lapply(1:D, 
                             function(d) {
                               tryCatch(expr  =               solve( s2y[[d]] ),
                                        error = function (e) { ginv( s2y[[d]] ) })
                             })


# variance-covariance matrix of random effects
s2u.true           <- c(v.ref.true, cor.ref12.true.num, cor.ref13.true.num, cor.ref23.true.num)
V.ud.true          <- matrix(c(v.ref.true[1],
                               rep(cor.ref12.true.num * sqrt( v.ref.true[1] * v.ref.true[2] ), 1),
                               rep(cor.ref13.true.num * sqrt( v.ref.true[1] * v.ref.true[3] ), 1),
                               #
                               rep(cor.ref12.true.num * sqrt( v.ref.true[1] * v.ref.true[2] ), 1),
                               v.ref.true[2],
                               rep(cor.ref23.true.num * sqrt( v.ref.true[2] * v.ref.true[3] ), 1),
                               #
                               rep(cor.ref13.true.num * sqrt( v.ref.true[1] * v.ref.true[3] ), 1),
                               rep(cor.ref23.true.num * sqrt( v.ref.true[2] * v.ref.true[3] ), 1),
                               v.ref.true[3]
),
m, m, byrow = TRUE)


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
method.scen = c("ML", "REML")

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


for (r in R.num_i) { # r = R.num_i[1]
  
  
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
  y.obs.r              <- y.true   + e
  
  
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
    
    for (var_k in 1:m) {
      
      obs.tmp  <- which(!is.na(y[, paste0("v",var_k)]))
      mis.tmp  <- which( is.na(y[, paste0("v",var_k)]))
      
      if (length(obs.tmp) > 1) {
        
        RES_UFH <- mseFH(formula = y[obs.tmp, var_k] ~ - 1 + x[[var_k]][obs.tmp,], 
                         vardir  = c(unlist(lapply(obs.tmp, function(d) { s2y[[d]][var_k, var_k] } ))), 
                         method  = method,
                         MAXITER = maxiter,
                         PRECISION = eps)
        
        if (!class(RES_UFH) == "try-error" ) {
          
          RES_param_theta[as.character(r), method, "UFH", var_k]                        <- RES_UFH$est$fit$refvar
          RES_param_beta[as.character(r), method, "UFH", beta.pos.for.uv[[var_k]]]      <- RES_UFH$est$fit$estcoef$beta
          
          RES_EBLUPS[as.character(r), method, "UFH", as.character(var_k), obs.tmp]      <- RES_UFH$est$eblup
          RES_EBLUP_MSE[as.character(r), method, "UFH", paste0("v", var_k), obs.tmp] <- RES_UFH$mse
          RES_EBLUP_MSE_alt[as.character(r), method, "UFH", paste0("v", var_k), obs.tmp] <- RES_UFH$mse
          
          RES_Syns[as.character(r), method, "UFH", as.character(var_k), mis.tmp]        <- x[[var_k]][obs.tmp,] %*% RES_UFH$est$fit$estcoef$beta
          
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
      
      RES_EBLUPS[as.character(r), method, "MFH", , obs.tmp]      <- t(RES_MFH$est$ebp)
      RES_EBLUP_MSE[as.character(r), method, "MFH", , obs.tmp] <- t(RES_MFH$mse)
      RES_EBLUP_MSE_alt[as.character(r), method, "MFH", , obs.tmp] <- t(RES_MFH$mse1Gd3)
      
      RES_Syns[as.character(r), method, "MFH", , mis.tmp]        <- t(sapply(1:m, function (var_k) {
        x[[var_k]][obs.tmp, ] %*% RES_param_beta[as.character(r), method, "MFH", beta.pos.for.uv[[var_k]]]
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
}


# ____________________________________________________________________________
# Save results ----
# ____________________________________________________________________________

paste0(path.user, "Dat/", name.sim, "/raw_sim_out/", "RES_param_beta", "/", name.data.gen, ".RData")

to_save <- grep("RES_", ls(), value = TRUE)
for ( i in to_save) {
  saveRDS(eval(parse(text = paste0(i))),
          file = paste0(path.user, "Dat/", name.sim, "/raw_sim_out/", i, "/", name.data.gen, ".RData"),
          version = 2)
}



# ____________________________________________________________________________
# warnings ----
# ____________________________________________________________________________

warnings()
cat("\n\n\n\n")

if (length(r.iter.na) > 0) {
  a.tmp <- length(r.iter.na); b.tmp <- ""
  if (a.tmp > 10) { a.tmp <- 10; b.tmp <- ", ..." }
  warning(paste0("\n \n \n MMFH: iter = NA in ",      name.data.gen, ", job ", job, ", for ", length(r.iter.na), "/", length(R.num_i), " cases", ", i.e. r = ", paste0(r.iter.na[1:a.tmp], collapse = ", "), b.tmp, "\n \n \n"))
  cat    (paste0("\n \n \n MMFH: iter = NA in ",      name.data.gen, ", job ", job, ", for ", length(r.iter.na), "/", length(R.num_i), " cases", ", i.e. r = ", paste0(r.iter.na[1:a.tmp], collapse = ", "), b.tmp, "\n \n \n"))
  saveRDS(r.iter.1,
          file = paste0(path.user, "Dat/", name.sim, "/raw_sim_out/", "r.iter.na", "/", name.data.gen, ".RData"),
          version = 2)
}
if (length(r.iter.1) > 0) {
  a.tmp <- length(r.iter.1); b.tmp <- ""
  if (a.tmp > 10) { a.tmp <- 10; b.tmp <- ", ..." }
  warning(paste0("\n \n \n MMFH: iter = 0 in ",       name.data.gen, ", job ", job, ", for ", length(r.iter.1), "/", length(R.num_i), " cases", ", i.e. r = ", paste0(r.iter.1[1:a.tmp], collapse = ", "), b.tmp, "\n \n \n"))
  cat    (paste0("\n \n \n MMFH: iter = 0 in ",       name.data.gen, ", job ", job, ", for ", length(r.iter.1), "/", length(R.num_i), " cases", ", i.e. r = ", paste0(r.iter.1[1:a.tmp], collapse = ", "), b.tmp, "\n \n \n"))
  saveRDS(r.iter.1,
          file = paste0(path.user, "Dat/", name.sim, "/raw_sim_out/", "r.iter.1", "/", name.data.gen, ".RData"),
          version = 2)
}
if (length(r.not.conv) > 0) {
  a.tmp <- length(r.not.conv); b.tmp <- ""
  if (a.tmp > 10) { a.tmp <- 10; b.tmp <- ", ..." }
  warning(paste0("\n \n \n MMFH not converged in ",   name.data.gen, ", job ", job, ", for ", length(r.not.conv), "/", length(R.num_i), " cases", ", i.e. r = ", paste0(r.not.conv[1:a.tmp], collapse = ", "), b.tmp, "\n \n \n"))
  cat    (paste0("\n \n \n MMFH not converged in ",   name.data.gen, ", job ", job, ", for ", length(r.not.conv), "/", length(R.num_i), " cases", ", i.e. r = ", paste0(r.not.conv[1:a.tmp], collapse = ", "), b.tmp, "\n \n \n"))
  saveRDS(r.not.conv,
          file = paste0(path.user, "Dat/", name.sim, "/raw_sim_out/", "r.not.conv", "/", name.data.gen, ".RData"),
          version = 2)
}
if (length(r.not.calc) > 0) {
  a.tmp <- length(r.not.calc); b.tmp <- ""
  if (a.tmp > 10) { a.tmp <- 10; b.tmp <- ", ..." }
  warning(paste0("\n \n \n MMFH not calculated in ",   name.data.gen, ", job ", job, ", for ", length(r.not.calc), "/", length(R.num_i), " cases", ", i.e. r = ", paste0(r.not.calc[1:a.tmp], collapse = ", "), b.tmp, "\n \n \n"))
  cat    (paste0("\n \n \n MMFH not calculated in ",   name.data.gen, ", job ", job, ", for ", length(r.not.calc), "/", length(R.num_i), " cases", ", i.e. r = ", paste0(r.not.calc[1:a.tmp], collapse = ", "), b.tmp, "\n \n \n"))
  saveRDS(r.not.calc,
          file = paste0(path.user, "Dat/", name.sim, "/raw_sim_out/", "r.not.calc", "/", name.data.gen, ".RData"),
          version = 2)
}


# ____________________________________________________________________________
# end code
# ____________________________________________________________________________

cat("End code \n")

end.code <- Sys.time()
print(difftime(end.code, start.code))


# ____________________________________________________________________________
# first results
# ____________________________________________________________________________

RES_param_beta[r,,,]
RES_param_theta[r,,,]
RES_EBLUPS[r,"REML",,,"1"]
RES_EBLUP_MSE[r,"REML",,,"1"]
RES_EBLUP_MSE_alt[r,"REML",,,"1"]
RES_Syns[r,"REML",,,"1"]
RES_SYN_MSE[r,"REML",,,"1"]

warnings()

print(name.data.gen)
