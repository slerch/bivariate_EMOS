rm(list=ls())

# load processed data
source("/home/sebastian/Projects/ATM_PP/code/postprocessing_tests/load_data.R")

## ------------------------------------------ ##

## Estimation of correlation function based on initial training period 

ind_trainCor <- which(dates <= 20160630)

source("/home/sebastian/Projects/ATM_PP/code/code_Nina/correlation_fun.R")

## following Nina's code:

# midpoint angles of sections of circle
theta <- c(22.5, 67.5, 112.5, 157.5, 202.5, 247.5, 292.5, 337.5)

# compute and plot empirical correlations
uwind_fc_means_train <- apply(uwind_fc_keep[,,,ind_trainCor], c(2,3,4), mean)
vwind_fc_means_train <- apply(vwind_fc_keep[,,,ind_trainCor], c(2,3,4), mean)
uwind_obs_train <- uwind_obs_keep[,,ind_trainCor]
vwind_obs_train <- vwind_obs_keep[,,ind_trainCor]

rad <- 8 # radius of inner cirle
cors <- correlation_fun(c(uwind_fc_means_train), 
                        c(uwind_obs_train),
                        c(vwind_fc_means_train), 
                        c(vwind_obs_train), 
                        radius = rad, main="", pch=20)

# set number of periods empirically
no.periods <- 2

# get predictor data and compute weights
rho.data <- cors[1,]
rho.weights <- cors[2,]/sum(cors[2,])

# estimate using weighted non-linear regression
rho.lm <- nls( rho.data ~ r * cos(2 * no.periods * pi/360 * theta + phi) + p, start=list(r=1, phi=0, p=0), weights=rho.weights)
rho.coef <- coef(rho.lm)

# function to compute correlation from angle
rho.fun <- function(x, no.periods){
  y <- rho.coef[1] * cos(2*no.periods*pi/360 * x + rho.coef[2]) + rho.coef[3]
  return(y)
}

# plot function
curve(rho.fun(x, no.periods), lty=2, lwd=2, add=TRUE)

# function to compute correlation from wind vector
rho.train.fun <- function(x.mean,y.mean, no.periods){
  phi.mean <- phi.help <- NULL
  for (i in 1:length(x.mean)){
    phi.help <- atan2(x.mean[i], y.mean[i])/(2*pi)*360
    if (phi.help < 0) phi.help <- 360 + phi.help
    phi.mean <- c(phi.mean,phi.help)
  }
  rho.mean <- rho.fun(phi.mean, no.periods)
  return(rho.mean)
}

# store for running EMOS
save(rho.fun, rho.train.fun, rho.coef, no.periods, file="/home/sebastian/Projects/ATM_PP/code/postprocessing_tests/localcorrelation.Rdata")
