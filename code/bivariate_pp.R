rm(list=ls())

library(MASS)
library(abind)
library(scoringRules)

# load processed data
source("/home/sebastian/Projects/ATM_PP/code/postprocessing_tests/load_data.R")
source("/home/sebastian/Projects/ATM_PP/code/helper_functions/multivariate_rank_histograms.R")

# load correlation data
load("/home/sebastian/Projects/ATM_PP/code/postprocessing_tests/localcorrelation.Rdata")

ind_eval <- which(dates >= 20160701) 

train_length <- 10

emos_coeffs_save <- NULL
uwind_pp_save <- array(NA, dim = c(50,201,81,length(ind_eval)))
vwind_pp_save <- array(NA, dim = c(50,201,81,length(ind_eval)))
es_raw_save <- NULL
es_pp_save <- NULL
vs_raw_save <- NULL
vs_pp_save <- NULL
avr_raw_save <- NULL
avr_pp_save <- NULL
bdr_raw_save <- NULL
bdr_pp_save <- NULL

for(ind_today in ind_eval){
  cat("\n", "Starting at", paste(Sys.time()), ":", ind_today, "of", 72,"\n"); flush(stdout())
  
  ind_training <- (ind_today-train_length):(ind_today-1)
  
  cat("... post-processing", "\n")
  # extract training u and v wind means and variances
  uwind_fc_means_train <- apply(uwind_fc_keep[,,,ind_training], c(2,3,4), mean)
  vwind_fc_means_train <- apply(vwind_fc_keep[,,,ind_training], c(2,3,4), mean)
  uwind_fc_vars_train <- apply(uwind_fc_keep[,,,ind_training], c(2,3,4), var)
  vwind_fc_vars_train <- apply(vwind_fc_keep[,,,ind_training], c(2,3,4), var)
  uwind_obs_train <- uwind_obs_keep[,,ind_training]
  vwind_obs_train <- vwind_obs_keep[,,ind_training]
  
  # linear regression to get mean vector
  U <- lm(c(uwind_obs_train) ~ c(uwind_fc_means_train))
  V <- lm(c(vwind_obs_train) ~ c(vwind_fc_means_train))
  U.coeffs <- U$coefficients
  V.coeffs <- V$coefficients
  
  # compute correlation of training data
  rho.training <- rho.train.fun(c(uwind_fc_means_train), c(vwind_fc_means_train), no.periods)
  
  # some terms to help with optimisation
  term.rho <- 1 - rho.training^2
  obs.mean.U <- c(uwind_obs_train) - (U.coeffs[1] + U.coeffs[2] * c(uwind_fc_means_train))
  obs.mean.V <- c(vwind_obs_train) - (V.coeffs[1] + V.coeffs[2] * c(vwind_fc_means_train)) 
  var.U <- c(uwind_fc_vars_train)
  var.V <- c(vwind_fc_vars_train)
  
  # log-likelihood function
  lik <- function(theta){
    gamma_U <- theta[1]; gamma_V <- theta[2]; delta_U <- theta[3]; delta_V <- theta[4]
    y <- sum( log(1/(2 * pi * sqrt(gamma_U^2 + delta_U^2 * var.U) *
                       sqrt(gamma_V^2 + delta_V^2 * var.V) * sqrt(term.rho)))
              - 0.5 * obs.mean.U^2 / term.rho / (gamma_U^2 + delta_U^2 * var.U)
              + rho.training * obs.mean.U * obs.mean.V
              / sqrt(gamma_U^2 + delta_U^2 * var.U)
              / sqrt(gamma_V^2 + delta_V^2 * var.V) / term.rho
              - 0.5 * obs.mean.V^2 / term.rho / (gamma_V^2 + delta_V^2 * var.V))
    return(y)
  }
  
  # optimise log-likelihood
  theta.start <- c(0.5,0.5,2,2)
  optim.out <- optim(theta.start, lik, method="BFGS", 
                     control=list(fnscale=-1, maxit=100000))
  
  # store output
  U.coeffs <- c(U.coeffs, optim.out$par[c(1,3)]^2)	
  V.coeffs <- c(V.coeffs, optim.out$par[c(2,4)]^2)
  emos.coeffs <- data.frame(t(U.coeffs), t(V.coeffs))
  names(emos.coeffs) <- c("a_U", "b_U", "c_U", "d_U", "a_V", "b_V", "c_V", "d_V")
  emos_coeffs_save <- rbind(emos_coeffs_save, emos.coeffs)
  
  cat("... generate post-processed ensemble", "\n")
  uwind_fc_means_test <- apply(uwind_fc_keep[,,,ind_today], c(2,3), mean)
  vwind_fc_means_test <- apply(vwind_fc_keep[,,,ind_today], c(2,3), mean)
  uwind_fc_vars_test <- apply(uwind_fc_keep[,,,ind_today], c(2,3), var)
  vwind_fc_vars_test <- apply(vwind_fc_keep[,,,ind_today], c(2,3), var)
  uwind_obs_test <- uwind_obs_keep[,,ind_today]
  vwind_obs_test <- vwind_obs_keep[,,ind_today]
  
  U.pred.mean <- U.coeffs[1] + U.coeffs[2] * c(uwind_fc_means_test)
  V.pred.mean <- V.coeffs[1] + V.coeffs[2] * c(vwind_fc_means_test)
  U.pred.sd <- sqrt(U.coeffs[3] + U.coeffs[4] * c(uwind_fc_vars_test))
  V.pred.sd <- sqrt(V.coeffs[3] + V.coeffs[4] * c(vwind_fc_vars_test))
  rho.pred <- rho.train.fun(U.pred.mean, V.pred.mean, no.periods)
  
  uwind_pp_test <- array(NA, dim = c(50,201,81))
  vwind_pp_test <- array(NA, dim = c(50,201,81))
  
  for(row in 1:201){
    # print(row)
    for(col in 1:81){
      pos_in_vector <- (col-1)*201 + row
      meanv <- c(U.pred.mean[pos_in_vector], V.pred.mean[pos_in_vector])
      covm <- matrix(NA, 2, 2)
      diag(covm) <- c(U.pred.sd[pos_in_vector]^2, V.pred.sd[pos_in_vector]^2)
      covm[1,2] <- covm[2,1] <- rho.pred[1]*U.pred.sd[pos_in_vector]*V.pred.sd[pos_in_vector]
      tmp <- mvrnorm(n = 50, mu = meanv, Sigma = covm)
      for(mem in 1:50){
        uwind_pp_test[mem, row, col] <- tmp[mem, 1] 
        vwind_pp_test[mem, row, col] <- tmp[mem, 2] 
      }
    }
  }
  
  # save
  uwind_pp_save[,,,which(ind_eval == ind_today)] <- uwind_pp_test
  vwind_pp_save[,,,which(ind_eval == ind_today)] <- vwind_pp_test
  
  cat("... compute scores", "\n")
  
  lonlat <- expand.grid(lon_values, lat_values)
  
  nIDs <- length(lon_values)*length(lat_values)
  es_out_raw <- numeric(length = nIDs)
  vs_out_raw <- numeric(length = nIDs)
  es_out_pp <- numeric(length = nIDs)
  vs_out_pp <- numeric(length = nIDs)
  
  uwind_obs_test <- uwind_obs_keep[,,ind_today]
  vwind_obs_test <- vwind_obs_keep[,,ind_today]
  
  for(id in 1:nIDs){
    # if(id %% 1000 == 0){print(id)}
    
    this_lon <- lonlat[id, 1]
    this_lat <- lonlat[id, 2]
    
    lonpos <- which(lon_values == this_lon)
    latpos <- which(lat_values == this_lat)
    
    this_obs <- c(uwind_obs_test[lonpos, latpos],
                  vwind_obs_test[lonpos, latpos])
    
    this_fc_raw <- rbind(uwind_fc_keep[,lonpos,latpos,ind_today], 
                         vwind_fc_keep[,lonpos,latpos,ind_today])
    this_fc_pp <- rbind(uwind_pp_test[,lonpos,latpos], 
                        vwind_pp_test[,lonpos,latpos])
    
    es_out_raw[id] <- es_sample(y = this_obs, dat = this_fc_raw)
    es_out_pp[id] <- es_sample(y = this_obs, dat = this_fc_pp)
    vs_out_raw[id] <- vs_sample(y = this_obs, dat = this_fc_raw)
    vs_out_pp[id] <- vs_sample(y = this_obs, dat = this_fc_pp)
  }
  
  es_raw_save <- c(es_raw_save, es_out_raw)
  es_pp_save <- c(es_pp_save, es_out_pp)
  vs_raw_save <- c(vs_raw_save, vs_out_raw)
  vs_pp_save <- c(vs_pp_save, vs_out_pp)
  
  cat("... compute ranks", "\n")
  
  fcraw_plot <- apply(abind(uwind_fc_keep[,,,ind_today], 
                            vwind_fc_keep[,,,ind_today], along = 4), 
                      c(1,4), rbind)
  fcpp_plot <- apply(abind(uwind_pp_test, vwind_pp_test, along = 4), c(1,4), rbind)
  
  obs_plot <- apply(abind(uwind_obs_test, vwind_obs_test, along = 3), c(3), rbind)
  
  avranks_raw <- avg_rh(fc = fcraw_plot, obs = obs_plot)
  avranks_pp <- avg_rh(fc = fcpp_plot, obs = obs_plot)
  
  avr_raw_save <- c(avr_raw_save, avranks_raw)
  avr_pp_save <- c(avr_pp_save, avranks_pp)
  
  bdranks_raw <- bd_rh(fc = fcraw_plot, obs = obs_plot)
  bdranks_pp <- bd_rh(fc = fcpp_plot, obs = obs_plot)
  
  bdr_raw_save <- c(bdr_raw_save, bdranks_raw)
  bdr_pp_save <- c(bdr_pp_save, bdranks_pp)
}
