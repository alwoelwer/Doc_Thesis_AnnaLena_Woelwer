# The code presents an exemplary application of the MMFH model
# MMFH: Multivariate Fay-Herriot model with missing dependent variables
# In the code, we generate data according to a MMFH model and
# calculate the parameters and predictions of the MMFH model for the data. 
# This code can e.g. be used for Monte Carlo simulation studies.


# Packages ----
# ______________________________________________________________________________

library(Matrix)
library(mvtnorm)
library(reshape2)
library(ggplot2)
library(sae)


# Source functions ----
# ______________________________________________________________________________

setwd("C:/Users/alwoe/Seafile/me/Uni/Proj/Dis/github")
source(paste0(getwd(), "/gen_MMFH_data_m3", ".R"))
source(paste0(getwd(), "/MMFH", ".R"))


# Generate example data ----
# ______________________________________________________________________________

# Generate the model information which is fixed
# See the description of the function for more information on the input.
# The function generates auxiliary information from a uniform distribution.
D = 100
m = 3
v_ref = c(2,3,4)
cor_ref = c(.2,.3,.4)
v_rer = c(2.5,3.5,4.5)
cor_rer = c(.15,.25,.35)
beta = list(c(1.5, 2.5), 
            c(2.3,3.3,4.3),
            c(4.1,3.1,2.1, 2.2))
range_aux = c(10, 100)
perc_mis = c(5, 5, 2)


d_fix <- f_prep_MMFH_m3(seed = 56)

# Generate the model information which typically varies between Monte Carlo
# iterations, That is, the generation of the dependent variables.
# We furthermore set certain values of the dependent variables as missing
# (missing completely at random MCAR).
# In the example code, the domains with missing values are non-overlapping.
# That is, maximum one dependent variable is missing per domain.
# This can easily be changed by editing the code
d_var <- f_gen_MMFH_m3( x = d_fix$x,
                        beta = d_fix$beta,
                        V_ud = d_fix$V_ud,
                        V_ed = d_fix$V_ed,
                        seed = 67,
                        verbose = TRUE )

# Let us have a look at the generated dependent variables
d_plot <- reshape2::melt(d_var$y_mis)
d_plot$Var2 <- factor(d_plot$Var2)
colnames(d_plot) <- c("domain", "variable", "value")

ggplot2::ggplot(d_plot,
                aes(x = domain, y = value, group = variable, 
                    color = variable, fill = variable, 
                    shape = variable, lty = variable)) + 
  geom_line() +
  geom_point()

ggplot2::ggplot(d_plot,
                aes(x = variable, y = value,
                    fill = variable)) + 
  geom_boxplot(alpha = .6) +
  geom_jitter(pch = 21, alpha = .6)

# From the error message, you can see that there are (as we wanted) missing
# values in the dependent variables.
# The number of domains with missing values of variables 1, 2, and 3 is
colSums(is.na(d_var$y_mis))



# Apply MMFH fitting algorithm ----
# ______________________________________________________________________________

theta <- NULL
theta_start <- NULL
method <- "REML"
# method <- "ML"
verbose <- TRUE
eps <- 1e-5
maxiter <- 100
y    = d_var$y_mis
x    = d_fix$x
V_ed = d_fix$V_ed

res_MMFH <- f_MMFH(y       = d_var$y_mis,
                   x       = d_fix$x,
                   V_ed    = d_fix$V_ed,
                   method  = method,
                   eps     = eps,
                   maxiter = maxiter,
                   verbose = verbose)


res_MMFH_dt <- f_MMFH_dt(
  y = d_var$y_mis,
  x = d_fix$x,
  V_ed = d_fix$V_ed,
  method = method,
  eps = eps,
  maxiter = maxiter,
  verbose = verbose
)



# Compare MMFH fitting algorithm with sae::eblupFH for univariate FH models ----
# ______________________________________________________________________________

res_FH <- list()
for (k in 1:m) {
  a_tmp <- which(!is.na(d_var$y_mis[,k]))
  res_FH[[k]] <- sae::mseFH(d_var$y_mis[,k][a_tmp] ~ -1 + d_fix$x[[k]][a_tmp,],
                            vardir = sapply(V_ed[a_tmp], function (d){ d[k,k] }),
                            method = method,
                            PRECISION = eps,
                            MAXITER = maxiter)
}

# For variable 1: 
# Compare true values, EBLUPs of (univariate) FH model and
# EBPs of MMFH model
k = 1
a_tmp <- which(!is.na(d_var$y_mis[,k]))
eblup_tmp <- rep(NA, D)
mse_tmp <- rep(NA, D)
eblup_tmp[which(!is.na(d_var$y_mis[,k]))] <- as.vector(res_FH[[k]]$est$eblup) 
mse_tmp[which(!is.na(d_var$y_mis[,k]))] <- as.vector(res_FH[[k]]$mse) 

# EBP
cbind("true" = d_var$y_true[,k],
      "FH_EBLUP" = eblup_tmp,
      "FH_SYN" = as.vector(x[[k]] %*% res_FH[[k]]$est$fit$estcoef[,1]),
      "MMFH" = res_MMFH$est$ebp[,k])[1:10,]

# MSE
cbind("FH" = mse_tmp,
      "MMFH" = res_MMFH$mse[,"v1"])[1:10,]

