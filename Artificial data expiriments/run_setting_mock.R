rm(list=ls())
library(dplyr)
library(tidyr)
library(scoringRules)
library(MASS)
library(FactorCopula)
library(statmod)
library(abind)
# "source" directory 

dir <- "C:/Users/20202943/Documents/BAM BEP/Simulation code/source code/"

source(paste0(dir, "generate_observations.R"))
source(paste0(dir, "generate_ensfc.R"))
source(paste0(dir, "postprocess_ensfc_mock.R"))
source(paste0(dir, "mvpp_mock.R"))
source(paste0(dir, "evaluation_functions.R"))



eval_all_mult <- function(mvpp_out, obs){
  esout <- es_wrapper(mvpp_out, obs)
  vs1out <- vs_wrapper(mvpp_out, obs, weight = FALSE, p = 1)
  vs1wout <- vs_wrapper(mvpp_out, obs, weight = TRUE, p = 1)
  vs0out <- vs_wrapper(mvpp_out, obs, weight = FALSE, p = 0.5)
  vs0wout <- vs_wrapper(mvpp_out, obs, weight = TRUE, p = 0.5)
  return(list("es" = esout, "vs1" = vs1out, "vs1w" = vs1wout, "vs0" = vs0out, "vs0w" = vs0wout))
}

run_setting_artificial <- function(obsmodel, fcmodel, nout, ninit, nmembers,timeWindow, MCrep, rand_rep, progress_ind = FALSE, compute_crps, ...){
  ecc_m <- NULL
  d <- list(...)$d
  # generate objects to save scores to
  modelnames <- c("ens", "emos.q", "ecc.q", "ecc.s", "decc.q", "ssh", "gca","Clayton","Frank","Gumbel", "GOF", "Clayton_shuffle", "Frank_shuffle", "Gumbel_shuffle")
  crps_list <- es_list <- vs1_list <- vs1w_list <- vs0_list <- vs0w_list <- param_list  <- list()
  timing_list <- list()
  for(mm in 1:length(modelnames)){
    es_list[[mm]] <- vs1_list[[mm]] <- vs1w_list[[mm]] <- 
      vs0_list[[mm]] <- vs0w_list[[mm]] <- param_list[[mm]]  <- matrix(NA, nrow = MCrep, ncol = nout) 
    crps_list[[mm]] <- array(NA, dim = c(MCrep, nout, d))
    timing_list[[mm]] <- array(NA, dim = MCrep) 
  }
  names(crps_list) <- names(es_list) <- names(vs1_list) <- names(vs1w_list) <- 
    names(vs0_list) <- names(vs0w_list) <- names(param_list)  <- names(timing_list) <- modelnames
  
  for(rr in 1:MCrep){
    if(progress_ind){
      if(rr %% 1 == 0){
        cat("starting at", paste(Sys.time()), ": MC repetition", rr, "of", MCrep, "\n"); flush(stdout())
      }
    }
    #data laden 
    load_data <- function(file_name = "Mock_data") {
      # Load artificial data
      all_data <- read.csv(paste0("C:/Users/20202943/Documents/BAM BEP/artificial data experiments/data provided by paper COBASE/", file_name, ".csv"))
      
      # transformations to help R understand the data
      all_data$date <- as.Date(all_data$validTime)
      all_data$station <- factor(all_data$station)
      
      # remove unused columns
      all_data$X <- NULL
      all_data$validTime <- NULL
      all_data$runtime <- NULL
      
      #add dates
      all_data <- all_data %>%
        arrange(date) %>%
        mutate(td = dense_rank(date))
      
      all_data$ytime <- all_data$date
      
      return(all_data)
    }
    data_mock <- load_data("Mock_data")
    # set random seed
    set.seed(110+rr)
    
    # generate observations
    start_time <- Sys.time()
    obs <- generate_obs2( data = data_mock,
                          file_name               = "Mock_data",
                          observation_columns     = c("obs"),
                          ensemble_regex          = c("^M_"),
                          store_uvpp = FALSE,
                          output_dim_standard     = 10)
    end_time <- Sys.time()
    timing_list$obs[rr] <- end_time - start_time

    # generate ensemble forecasts
    start_time <- Sys.time()
    fc <- generate_ensfc2( data = data_mock,
                                       file_name               = "Mock_data",
                                       observation_columns     = c("obs"),
                                       ensemble_regex          = c("^M_"),
                                       store_uvpp = FALSE,
                                       output_dim_standard     = 10)

    end_time <- Sys.time()
    timing_list$fc[rr] <- end_time - start_time
    
    #postprocess ensemble forecasts
    start_time <- Sys.time()
    pp_out <- postproc(fcmodel = fcmodel, ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
                       obs = obs$obs, obs_init = obs$obs_init,
                       train = "init", trainlength = timeWindow, emos_plus = TRUE)
    end_time <- Sys.time()
    timing_list$uvpp[rr] <- end_time - start_time
    
    # if there are NaN's in pp output, generate new sets of obs and fc, and re-run pp code
    while(anyNA(pp_out)){
      set.seed(sample(1:100,1))
      obs <- generate_obs2( data = data_mock,
                            file_name               = "Mock_data",
                            observation_columns     = c("obs"),
                            ensemble_regex          = c("^M_"),
                            store_uvpp = FALSE,
                            output_dim_standard     = 10)
      fc <- generate_ensfc2( data = data_mock,
                             file_name               = "Mock_data",
                             observation_columns     = c("obs"),
                             ensemble_regex          = c("^M_"),
                             store_uvpp = FALSE,
                             output_dim_standard     = 10)
      pp_out <- postproc(fcmodel = fcmodel, ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
                         obs = obs$obs, obs_init = obs$obs_init,
                         train = "init", trainlength = TimeWindow, emos_plus = TRUE)
    }

    start_time <- Sys.time()
    if(compute_crps){
      crps_list$ens[rr, , ] <- crps_wrapper(fc$ensfc, obs$obs)
    }
    tmp <- eval_all_mult(mvpp_out = fc$ensfc, obs = obs$obs)
    es_list$ens[rr, ] <- tmp$es
    vs1_list$ens[rr, ] <- tmp$vs1
    vs1w_list$ens[rr, ] <- tmp$vs1w
    vs0_list$ens[rr, ] <- tmp$vs0
    vs0w_list$ens[rr, ] <- tmp$vs0w
    
    end_time <- Sys.time()
    
    timing_list$ens[rr] <- end_time - start_time
    # EMOS.Q
    start_time <- Sys.time()
    emos.q <- mvpp(method = "EMOS", variant = "Q", ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
                   obs = obs$obs, obs_init = obs$obs_init, postproc_out = pp_out, timeWindow = timeWindow)
 
    if(compute_crps){
      crps_list$emos.q[rr, , ] <- crps_wrapper(emos.q$mvppout, obs$obs)
    }
    
    
    tmp <- eval_all_mult(mvpp_out = emos.q$mvppout, obs = obs$obs)
    es_list$emos.q[rr, ] <- tmp$es
    vs1_list$emos.q[rr, ] <- tmp$vs1
    vs1w_list$emos.q[rr, ] <- tmp$vs1w
    vs0_list$emos.q[rr, ] <- tmp$vs0
    vs0w_list$emos.q[rr, ] <- tmp$vs0w
    
    end_time <- Sys.time()
    
    timing_list$emos.q[rr] <- end_time - start_time
    cat("prima tot hier ", "\n"); flush(stdout())
    # ECC.Q
    start_time <- Sys.time()
    ecc.q <- mvpp(method = "ECC", ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
                  obs = obs$obs, obs_init = obs$obs_init, postproc_out = pp_out,
                  EMOS_sample = emos.q$mvppout, timeWindow = timeWindow)
    if(compute_crps){
      crps_list$ecc.q[rr, , ] <- crps_wrapper(ecc.q$mvppout, obs$obs)
    }
    tmp <- eval_all_mult(mvpp_out = ecc.q$mvppout, obs = obs$obs)
    es_list$ecc.q[rr, ] <- tmp$es
    vs1_list$ecc.q[rr, ] <- tmp$vs1
    vs1w_list$ecc.q[rr, ] <- tmp$vs1w
    vs0_list$ecc.q[rr, ] <- tmp$vs0
    vs0w_list$ecc.q[rr, ] <- tmp$vs0w
    
    end_time <- Sys.time()
    
    timing_list$ecc.q[rr] <- end_time - start_time
    #
    # ECC.S 
    es_list_tmp <- vs1_list_tmp <- vs1w_list_tmp <-
      vs0_list_tmp <- vs0w_list_tmp <- matrix(NA, nrow = nout, ncol = rand_rep)
    crps_list_tmp <- array(NA, dim = c(nout, d, rand_rep))
    timing_list_tmp <- array(NA, dim = rand_rep)
    for(RR in 1:rand_rep){
      start_time <- Sys.time()
      emos.s <- mvpp(method = "EMOS", variant = "S", ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
                     obs = obs$obs, obs_init = obs$obs_init, postproc_out = pp_out, timeWindow = timeWindow)
      ecc.s <- mvpp(method = "ECC", ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
                    obs = obs$obs, obs_init = obs$obs_init, postproc_out = pp_out,
                    EMOS_sample = emos.s$mvppout, timeWindow = timeWindow)
      if(compute_crps){
        crps_list_tmp[,,RR] <- crps_wrapper(ecc.s$mvppout, obs$obs)
      }
      tmp <- eval_all_mult(mvpp_out = ecc.s$mvppout, obs = obs$obs)
      es_list_tmp[,RR] <- tmp$es
      vs1_list_tmp[,RR] <- tmp$vs1
      vs1w_list_tmp[,RR] <- tmp$vs1w
      vs0_list_tmp[,RR] <- tmp$vs0
      vs0w_list_tmp[,RR] <- tmp$vs0w
      
      end_time <- Sys.time()
      
      timing_list_tmp[RR] <- end_time - start_time
    }
    
    crps_list$ecc.s[rr,,] <- apply(crps_list_tmp, c(1,2), mean)
    es_list$ecc.s[rr, ] <- apply(es_list_tmp, 1, mean)
    vs1_list$ecc.s[rr, ] <- apply(vs1_list_tmp, 1, mean)
    vs1w_list$ecc.s[rr, ] <- apply(vs1w_list_tmp, 1, mean)
    vs0_list$ecc.s[rr, ] <- apply(vs0_list_tmp, 1, mean)
    vs0w_list$ecc.s[rr, ] <- apply(vs0w_list_tmp, 1, mean)
    timing_list$ecc.s[rr] <- apply(timing_list_tmp, 1, mean)
    
    # dECC.Q
    start_time <- Sys.time()
    decc.q <- mvpp(method = "dECC", ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
                   obs = obs$obs, obs_init = obs$obs_init, postproc_out = pp_out,
                   EMOS_sample = emos.q$mvppout, ECC_out = ecc.q$mvppout, timeWindow = timeWindow)
    if(compute_crps){
      crps_list$decc.q[rr, , ] <- crps_wrapper(decc.q$mvppout, obs$obs)
    }
    tmp <- eval_all_mult(mvpp_out = decc.q$mvppout, obs = obs$obs)
    es_list$decc.q[rr, ] <- tmp$es
    vs1_list$decc.q[rr, ] <- tmp$vs1
    vs1w_list$decc.q[rr, ] <- tmp$vs1w
    vs0_list$decc.q[rr, ] <- tmp$vs0
    vs0w_list$decc.q[rr, ] <- tmp$vs0w
    
    end_time <- Sys.time()
    
    timing_list$decc.q[rr] <- end_time - start_time
    
    # SSh 
    es_list_tmp <- vs1_list_tmp <- vs1w_list_tmp <-
      vs0_list_tmp <- vs0w_list_tmp <- matrix(NA, nrow = nout, ncol = rand_rep)
    crps_list_tmp <- array(NA, dim = c(nout, d, rand_rep))
    timing_list_tmp <- array(NA, dim = rand_rep)
    for(RR in 1:rand_rep){
      start_time <- Sys.time()
      ssh <- mvpp(method = "SSh", ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
                  obs = obs$obs, obs_init = obs$obs_init, postproc_out = pp_out,
                  EMOS_sample = emos.q$mvppout, timeWindow = timeWindow)
      if(compute_crps){
        crps_list_tmp[,,RR] <- crps_wrapper(ssh$mvppout, obs$obs)
      }
      tmp <- eval_all_mult(mvpp_out = ssh$mvppout, obs = obs$obs)
      es_list_tmp[,RR] <- tmp$es
      vs1_list_tmp[,RR] <- tmp$vs1
      vs1w_list_tmp[,RR] <- tmp$vs1w
      vs0_list_tmp[,RR] <- tmp$vs0
      vs0w_list_tmp[,RR] <- tmp$vs0w
      
      end_time <- Sys.time()
      
      timing_list_tmp[RR] <- end_time - start_time
    }
    crps_list$ssh[rr,,] <- apply(crps_list_tmp, c(1,2), mean)
    es_list$ssh[rr, ] <- apply(es_list_tmp, 1, mean)
    vs1_list$ssh[rr, ] <- apply(vs1_list_tmp, 1, mean)
    vs1w_list$ssh[rr, ] <- apply(vs1w_list_tmp, 1, mean)
    vs0_list$ssh[rr, ] <- apply(vs0_list_tmp, 1, mean)
    vs0w_list$ssh[rr, ] <- apply(vs0w_list_tmp, 1, mean)
    timing_list$ssh[rr] <- apply(timing_list_tmp, 1, mean)
    
   
    # Archimedean copulas 
    for (method in c("Clayton","Frank", "Gumbel")) { #GOF even weggehaalt
      es_list_tmp <- vs1_list_tmp <- vs1w_list_tmp <-
        vs0_list_tmp <- vs0w_list_tmp <- param_list_temp  <- matrix(NA, nrow = nout, ncol = rand_rep)
      crps_list_tmp <- array(NA, dim = c(nout, d, rand_rep))
      timing_list_tmp <- array(NA, dim = rand_rep)
      for(RR in 1:rand_rep){
        start_time <- Sys.time()
        mvd <- mvpp(method = method, ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
                    obs = obs$obs, obs_init = obs$obs_init, postproc_out = pp_out, timeWindow = timeWindow)
        if(compute_crps){
          crps_list_tmp[,,RR] <- crps_wrapper(mvd$mvppout, obs$obs)
        }
        
      
        
        tmp <- eval_all_mult(mvpp_out = mvd$mvppout, obs = obs$obs)
        es_list_tmp[,RR] <- tmp$es
        vs1_list_tmp[,RR] <- tmp$vs1
        vs1w_list_tmp[,RR] <- tmp$vs1w
        vs0_list_tmp[,RR] <- tmp$vs0
        vs0w_list_tmp[,RR] <- tmp$vs0w
        param_list_temp[,RR] <- mvd$params
        
        
        
        end_time <- Sys.time()
        
        timing_list_tmp[RR] <- end_time - start_time
      }
      
      crps_list[[method]][rr,,] <- apply(crps_list_tmp, c(1,2), mean)
      es_list[[method]][rr, ] <- apply(es_list_tmp, 1, mean)
      vs1_list[[method]][rr, ] <- apply(vs1_list_tmp, 1, mean)
      vs1w_list[[method]][rr, ] <- apply(vs1w_list_tmp, 1, mean)
      vs0_list[[method]][rr, ] <- apply(vs0_list_tmp, 1, mean)
      vs0w_list[[method]][rr, ] <- apply(vs0w_list_tmp, 1, mean)
      param_list[[method]][rr, ] <- apply(param_list_temp, 1, mean)
      timing_list[[method]][rr] <- apply(timing_list_tmp, 1, mean)
    }
    
   
    # Archimedean copulas shuffle
    mvppout_all <- list(Clayton = vector("list", rand_rep),
                        Frank   = vector("list", rand_rep),
                        Gumbel  = vector("list", rand_rep),
                        GOF     = vector("list", rand_rep))
    for (method in c("Clayton","Frank", "Gumbel")) { #GOF even weggehaalt
      es_list_tmp <- vs1_list_tmp <- vs1w_list_tmp <-
        vs0_list_tmp <- vs0w_list_tmp <- param_list_temp <- matrix(NA, nrow = nout, ncol = rand_rep)
      crps_list_tmp <- array(NA, dim = c(nout, d, rand_rep))
      timing_list_tmp <- array(NA, dim = rand_rep)
      for(RR in 1:rand_rep){
        start_time <- Sys.time()
        mvd <- mvpp(method = method, ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
                    obs = obs$obs, obs_init = obs$obs_init, postproc_out = pp_out, timeWindow = timeWindow)
        
        method_shuffle <- paste0(method, "_shuffle")
        COBASE <- mvpp(method = method_shuffle, ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
                       EMOS_sample = emos.q$mvppout, obs = obs$obs, obs_init = obs$obs_init, postproc_out = pp_out, MVPP_sample = mvd$mvpp, timeWindow = timeWindow)
        
        
        if(compute_crps){
          crps_list_tmp[,,RR] <- crps_wrapper(COBASE$mvppout, obs$obs)
        }
     
        
        tmp_COBASE <- eval_all_mult(mvpp_out = COBASE$mvppout, obs = obs$obs)
        es_list_tmp[,RR] <- tmp_COBASE$es
        vs1_list_tmp[,RR] <- tmp_COBASE$vs1
        vs1w_list_tmp[,RR] <- tmp_COBASE$vs1w
        vs0_list_tmp[,RR] <- tmp_COBASE$vs0
        vs0w_list_tmp[,RR] <- tmp_COBASE$vs0w
      
        
        end_time <- Sys.time()
        
        timing_list_tmp[RR] <- end_time - start_time
      }
      
      crps_list[[paste0(method, "_shuffle")]][rr,,] <- apply(crps_list_tmp, c(1,2), mean)
      es_list[[paste0(method, "_shuffle")]][rr, ] <- apply(es_list_tmp, 1, mean)
      vs1_list[[paste0(method, "_shuffle")]][rr, ] <- apply(vs1_list_tmp, 1, mean)
      vs1w_list[[paste0(method, "_shuffle")]][rr, ] <- apply(vs1w_list_tmp, 1, mean)
      vs0_list[[paste0(method, "_shuffle")]][rr, ] <- apply(vs0_list_tmp, 1, mean)
      vs0w_list[[paste0(method, "_shuffle")]][rr, ] <- apply(vs0w_list_tmp, 1, mean)
      param_list[[paste0(method, "_shuffle")]][rr, ] <- apply(param_list_temp, 1, mean)
      timing_list[[paste0(method, "_shuffle")]][rr] <- apply(timing_list_tmp, 1, mean)
    }
   
    
  }
  
  # return results
  out <- list("crps_list" = crps_list, "es_list" = es_list, "vs1_list" = vs1_list,
              "vs1w_list" = vs1w_list, "vs0_list" = vs0_list, "vs0w_list" = vs0w_list, "param_list" = param_list, "timing_list" = timing_list)
  return(out)
}


input_copula <- c("Frank","Gumbel", "Clayton", "Frank_shuffle", "Gumbel_shuffle", "Clayton_shuffle")  
input_d <- 3
repetitions <- 1
eps <- 1
sigma <- sqrt(5)
input_par <- expand.grid( input_copula, eps, sigma, input_d, 1:repetitions)
names(input_par) <- c("copula", "eps", "sigma", "d", "repetition")
MC_reps <- 1
evalDays <- 240
timeWindow <- 30
trainingDays <- 60
ensembleMembers <- 10
randomRepetitions <- 1
Rdata_dir <- "C:/Users/20202943/Documents/BAM BEP/artificial data experiments/data files/" # directory to save Rdata files to
Rout_dir <- "C:/Users/20202943/Documents/BAM BEP/artificial data experiments/data files/"  # directory to save Rout files to
run_wrapper <- function(runID){
  sink(file = paste0(Rout_dir, "setting_mock_with_timewindow_", runID, ".Rout"))
  tau <- runif(1, 0, 1)
  par_values <- as.numeric(input_par[ID, ])
  res <- run_setting_artificial(obsmodel = 7, fcmodel = 7, nout = evalDays, ninit = trainingDays, 
                      nmembers = ensembleMembers,
                      MCrep = MC_reps, rand_rep = randomRepetitions, 
                      progress_ind = TRUE, compute_crps = TRUE,
                      eps = input_par$eps[runID], 
                      sigma = input_par$sigma[runID], 
                      copula = input_par$copula[runID],
                      d = input_par$d[runID], timeWindow = timeWindow)
  savename <- paste0(Rdata_dir, "res_setting_mock_with_timewindow_", runID, ".Rdata")
  save(res, input_par, file = savename)
  sink()

}
ID <- 1
run_wrapper(runID = ID)
