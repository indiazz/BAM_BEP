#------------------set up--------------------
library(scoringRules)
library(MASS)
library(multDM)
library(dplyr)
library(tidyr)
# forecasts,observations and score sources
dir <- "C:/Users/20202943/Desktop/TUE/Year 6/BEP BAM/multiv_pp-master/multiv_pp-master/simulation code/source/"
source(paste0(dir, "generate_observations.R")) #niet nodig
source(paste0(dir, "generate_ensfc.R")) #niet nodig
source(paste0(dir, "postprocess_ensfc.R")) #neit nodig?
source(paste0(dir, "mvpp.R")) #niet nodig
source(paste0(dir, "evaluation_functions.R")) #nodig!
# MDM scources
dir <- "C:/Users/20202943/Documents/BAM BEP/"
source(paste0(dir, "MDM_scoretest.R"))
source(paste0(dir, "nieuwe MDM selection.R"))

set.seed(123)
MCrep <- 1000
output_all <- vector("list", MCrep)

#begin MC reps
for (mcrep in 1:MCrep){
  # make observations
  require(MASS)
  nout <- 1000
  ninit <- 500
  d <- 3
  rho0 <- c(0.1)
  obs_init <- array(NA, dim = c(ninit, d))
  obs <- array(NA, dim = c(nout, d))
  S <- matrix(NA, d, d)
  for(i in 1:d){
    for(j in 1:d){
      S[i,j] <- rho0^(abs(i-j))
    }
  }
  # mean vector
  mu <- rep(0, d)
  
  # observations
  obs_init <- mvrnorm(n = ninit, mu = mu, Sigma = S)
  obs <- mvrnorm(n = nout, mu = mu, Sigma = S)
  
  observations <- list("obs_init" = obs_init, "obs" = obs)
  
  # make raw ensemble forecasts
  
  samples <- c(1,2,3,4)
  for(s in samples){
    m <- 50
    
    # initialize output arrays
    ensfc_init <- array(NA, dim = c(ninit, m, d))
    ensfc <- array(NA, dim = c(nout, m, d))
    
    # Setting size -> no misspecification
    eps <- 0
    rho <- 0.1
    sigma <- 1
    # correlation matrix
    S <- matrix(NA, d, d)
    for(i in 1:d){
      for(j in 1:d){
        S[i,j] <- sigma*rho^(abs(i-j))
      }
    }
    
    # mean vector
    mu <- rep(eps, d)
    
    # generate forecasts
    for(nn in 1:ninit){
      tmp <- mvrnorm(n = m, mu = mu, Sigma = S)
      ensfc_init[nn,,] <- tmp
    }
    for(nn in 1:nout){
      tmp <- mvrnorm(n = m, mu = mu, Sigma = S)
      ensfc[nn,,] <- tmp
    }
    assign(paste0("sample_ensemble", s), list("ensfc_init" = ensfc_init, "ensfc" = ensfc))
  }
  #misspecified samples 5 and 6
  samples <- c(5,6)
  for(s in samples){
    m <- 50
    
    # initialize output arrays
    ensfc_init <- array(NA, dim = c(ninit, m, d))
    ensfc <- array(NA, dim = c(nout, m, d))
    
    # misspecification setting
    eps <- 0 #1
    rho <- 0.5 #0.1
    sigma <- 1
    # correlation matrix
    S <- matrix(NA, d, d)
    for(i in 1:d){
      for(j in 1:d){
        S[i,j] <- sigma*rho^(abs(i-j))
      }
    }
    
    # mean vector
    mu <- rep(eps, d)
    
    # generate forecasts
    for(nn in 1:ninit){
      tmp <- mvrnorm(n = m, mu = mu, Sigma = S)
      ensfc_init[nn,,] <- tmp
    }
    for(nn in 1:nout){
      tmp <- mvrnorm(n = m, mu = mu, Sigma = S)
      ensfc[nn,,] <- tmp
    }
    assign(paste0("sample_ensemble", s), list("ensfc_init" = ensfc_init, "ensfc" = ensfc))
  }
 
  #all samples
  samples <- list(
    "1_sample_ensemble" = sample_ensemble1, "2_sample_ensemble" = sample_ensemble2, "3_sample_ensemble" = sample_ensemble3,"4_sample_ensemble" = sample_ensemble4, "5_sample_ensemble" = sample_ensemble5, "6_sample_ensemble" = sample_ensemble6)
  samplenames <- names(samples)
  #score lists
  crps_list <- es_list <- vs1_list <- vs1w_list <- vs0_list <- vs0w_list <- obs_list <- vector("list", length(samplenames))
  names(crps_list) <- names(es_list) <- names(vs1_list) <- names(vs1w_list) <- names(vs0_list) <- names(vs0w_list) <- names(obs_list) <- samplenames
  
  for (nm in samplenames) {
    es_list[[nm]]  <- vs1_list[[nm]] <- vs1w_list[[nm]] <- vs0_list[[nm]] <- vs0w_list[[nm]] <- matrix(NA_real_, nrow = 1, ncol = nout)
    crps_list[[nm]] <- obs_list[[nm]] <- array(NA_real_, dim = c(1, nout, d))
  }
  
  eval_all_mult <- function(mvpp_out, obs){
    esout  <- es_wrapper(mvpp_out, obs, return_obs_list = TRUE)
    vs1out <- vs_wrapper(mvpp_out, obs, weight = FALSE, p = 1)
    vs1wout<- vs_wrapper(mvpp_out, obs, weight = TRUE,  p = 1)
    vs0out <- vs_wrapper(mvpp_out, obs, weight = FALSE, p = 0.5)
    vs0wout<- vs_wrapper(mvpp_out, obs, weight = TRUE,  p = 0.5)
    list(es = esout, vs1 = vs1out, vs1w = vs1wout, vs0 = vs0out, vs0w = vs0wout)
  }
  
  #score computations
  for (nm in samplenames) {
    ens <- samples[[nm]]$ensfc
    
    crps_list[[nm]][1,,] <- crps_wrapper(ens, observations$obs)
    
    tmp <- eval_all_mult(mvpp_out = ens, obs = observations$obs)
    
    obs_list[[nm]][1,,]  <- tmp$es$obs_list
    es_list[[nm]][1, ]   <- tmp$es$es
    vs1_list[[nm]][1, ]  <- tmp$vs1
    vs1w_list[[nm]][1, ] <- tmp$vs1w
    vs0_list[[nm]][1, ]  <- tmp$vs0
    vs0w_list[[nm]][1, ] <- tmp$vs0w
  }

  score_lists <- list(
    vs0  = vs0_list,
    vs0w = vs0w_list,
    vs1  = vs1_list,
    vs1w = vs1w_list,
    ES   = es_list
  )
  #choose which models are considered for size and power
  models_size  <- c("1_sample_ensemble","2_sample_ensemble","3_sample_ensemble","4_sample_ensemble")
  models_power <- samplenames
  model_sets <- list(size = models_size, power = models_power)
  #MDM over samples
  chunks <- vector("list", length(score_lists) * length(model_sets))
  k <- 1
  for (type in names(model_sets)){
    models_use <- model_sets[[type]]
  for (sc in names(score_lists)) {
    # filter models
    S_filt <- score_lists[[sc]][models_use]
    
    S_matrix <- do.call(rbind, S_filt)
    rownames(S_matrix) <- models_use
    stopifnot(nrow(S_matrix) == length(models_use))
    stopifnot(ncol(S_matrix) == nout)
    q_val <- floor(ncol(S_matrix)^(1/3))
    result <- MDM_on_scores_chain2( S_matrix = S_matrix, q= q_val,statistic = "Sc")
    result_selection <- selection(alpha = 0.05,q= q_val,statistic = "Sc",S_matrix  = S_matrix)
    
    outcomes <- result_selection$outcomes
    best_id  <- which.min(outcomes[, "Mean loss"])
    best_model <- rownames(outcomes)[best_id]
    
    # statistics
    T_value <- if (!is.null(result$parameter["T"])) {
      as.numeric(result$parameter["T"])
    } else {
      ncol(S_matrix)
    }
    
    df_value <- if (!is.null(result$parameter["df"])) {
      as.numeric(result$parameter["df"])
    } else {
      NA_real_
    }
    
    chunks[[k]] <- data.frame(
      type       = type,      
      score      = sc,
      n_models   = length(models_use),
      T          = T_value,
      q_value    = q_val,
      df         = df_value,
      Statistic  = as.numeric(result$statistic),
      p_value    = as.numeric(result$p.value),
      best_model = best_model,
      stringsAsFactors = FALSE
    )
    k <- k + 1
  }
}
  out_mcrep <- do.call(rbind, chunks)
  out_mcrep$mcrep <- mcrep
  output_all[[mcrep]] <- out_mcrep
  print(mcrep)
}
output_all_dfp <- do.call(rbind, output_all)
row.names(output_all_dfp) <- NULL
save(output_all_dfp, file = paste0("C:/Users/20202943/Documents/BAM BEP/","_power_MC_MDM_mock_powerdependence_SC_",".Rdata"))

#Analyses for RQ3
load("C:/Users/20202943/Documents/BAM BEP/_power_MC_MDM_mock_powerdependence_SC_.Rdata")
alpha <- 0.05
summary_rejections <- output_all_dfp %>%
  mutate(signif = p_value < alpha) %>%
  group_by(type, score) %>%
  summarise(rejection_rate = mean(signif), .groups = "drop")

print(summary_rates)

alpha <- 0.05
output_all_dfp$reject <- output_all_dfp$p_value < alpha

#power
powerss <- output_all_dfp %>%
  filter(type == "power")
empirical_power <- aggregate(reject ~ score, data = powerss, FUN = mean)
empirical_power$MCSE <- with(empirical_power, sqrt(reject * (1 - reject) / MCrep))

empirical_power

#size
size_sum <- summary_rates %>%
  filter(type == "size")

sizeee <- output_all_dfp %>%
  filter(type == "size")
empirical_size <- aggregate(reject ~ score, data = sizeee, FUN = mean)
empirical_size$MCSE <- with(empirical_size, sqrt(reject * (1 - reject) / MCrep))

empirical_size

