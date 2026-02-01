rm(list=ls())


library(scoringRules)
library(MASS)
library(FactorCopula)
library(statmod)
library(abind)
# "source" directory 
dir <- "C:/Users/20202943/Documents/BAM BEP/CODE/Simulation code/source code/"

source(paste0(dir, "generate_observations.R")) 
source(paste0(dir, "generate_ensfc.R")) 
source(paste0(dir, "postprocess_ensfc_arch.R"))
source(paste0(dir, "mvpp_arch_shuffle.R")) 
source(paste0(dir, "evaluation_functions.R")) 




eval_all_mult <- function(mvpp_out, obs){
  esout <- es_wrapper(mvpp_out, obs)
  vs1out <- vs_wrapper(mvpp_out, obs, weight = FALSE, p = 1)
  vs1wout <- vs_wrapper(mvpp_out, obs, weight = TRUE, p = 1)
  vs0out <- vs_wrapper(mvpp_out, obs, weight = FALSE, p = 0.5)
  vs0wout <- vs_wrapper(mvpp_out, obs, weight = TRUE, p = 0.5)
  return(list("es" = esout, "vs1" = vs1out, "vs1w" = vs1wout, "vs0" = vs0out, "vs0w" = vs0wout))
}

most_occurring <- function(x) {
  return(names(sort(table(x),decreasing=TRUE))[1])
}


run_setting6 <- function(obsmodel, fcmodel, nout, ninit, nmembers, MCrep, rand_rep, progress_ind = FALSE, compute_crps, ...){
  ecc_m <- NULL
  d <- list(...)$d
  # generate objects to save scores to
  modelnames <- c("ens", "emos.q", "ecc.q", "ecc.s", "decc.q", "ssh", "gca","Clayton","Frank","Gumbel", "GOF", "Clayton_shuffle", "Frank_shuffle", "Gumbel_shuffle")
  crps_list <- es_list <- vs1_list <- vs1w_list <- vs0_list <- vs0w_list <- param_list <- list()
  # Timing_list stores the time needed to do all relevant computations
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
    
    
    # set random seed
    set.seed(110+rr)
    
    # generate observations
    start_time <- Sys.time()
    obs <- generate_obs(model = obsmodel, nout = nout, ninit = ninit, ...)
    end_time <- Sys.time()
    timing_list$obs[rr] <- end_time - start_time
    # generate ensemble forecasts
    start_time <- Sys.time()
    fc <- generate_ensfc(model = fcmodel, nout = nout, ninit = ninit, nmembers = nmembers, ...)
    end_time <- Sys.time()
    timing_list$fc[rr] <- end_time - start_time
    #cat("gelukt ", "\n"); flush(stdout())
    
    #postprocess ensemble forecasts
    start_time <- Sys.time()
    pp_out <- postproc(fcmodel = fcmodel, ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
                       obs = obs$obs, obs_init = obs$obs_init,
                       train = "init", trainlength = NULL, emos_plus = TRUE)
    end_time <- Sys.time()
    timing_list$uvpp[rr] <- end_time - start_time
    
    # if there are NaN's in pp output -> generate new sets of observations and ensemble forecasts -> re-run pp code
    while(anyNA(pp_out)){
      set.seed(sample(1:100,1))
      obs <- generate_obs(model = obsmodel, nout = nout, ninit = ninit, ...)
      fc <- generate_ensfc(model = fcmodel, nout = nout, ninit = ninit, nmembers = nmembers, ...)
      pp_out <- postproc(fcmodel = fcmodel, ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
                         obs = obs$obs, obs_init = obs$obs_init,
                         train = "init", trainlength = NULL, emos_plus = TRUE)
      print("Repeating observations")
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
    cat("prima tot hier ", "\n"); flush(stdout())
    # EMOS.Q
    start_time <- Sys.time()
    emos.q <- mvpp(method = "EMOS", variant = "Q", ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
                   obs = obs$obs, obs_init = obs$obs_init, postproc_out = pp_out)
    
    
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
    #cat("prima tot hier ", "\n"); flush(stdout())
    # ECC.Q
    start_time <- Sys.time()
    ecc.q <- mvpp(method = "ECC", ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
                  obs = obs$obs, obs_init = obs$obs_init, postproc_out = pp_out,
                  EMOS_sample = emos.q$mvppout)
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
    
    # ECC.S 
    es_list_tmp <- vs1_list_tmp <- vs1w_list_tmp <-
      vs0_list_tmp <- vs0w_list_tmp <- matrix(NA, nrow = nout, ncol = rand_rep)
    crps_list_tmp <- array(NA, dim = c(nout, d, rand_rep))
    timing_list_tmp <- array(NA, dim = rand_rep)
    for(RR in 1:rand_rep){
      start_time <- Sys.time()
      emos.s <- mvpp(method = "EMOS", variant = "S", ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
                     obs = obs$obs, obs_init = obs$obs_init, postproc_out = pp_out)
      ecc.s <- mvpp(method = "ECC", ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
                    obs = obs$obs, obs_init = obs$obs_init, postproc_out = pp_out,
                    EMOS_sample = emos.s$mvppout)
      print("  check!!")
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
                   EMOS_sample = emos.q$mvppout, ECC_out = ecc.q$mvppout)
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
                  EMOS_sample = emos.q$mvppout)
      if(compute_crps){
        crps_list_tmp[,,RR] <- crps_wrapper(ssh$mvppout, obs$obs)
      }
      #print("second checkpoint")
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
    
    # GCA
    es_list_tmp <- vs1_list_tmp <- vs1w_list_tmp <-
      vs0_list_tmp <- vs0w_list_tmp <- matrix(NA, nrow = nout, ncol = rand_rep)
    crps_list_tmp <- array(NA, dim = c(nout, d, rand_rep))
    timing_list_tmp <- array(NA, dim = rand_rep)
    for(RR in 1:rand_rep){
      start_time <- Sys.time()
      gca <- mvpp(method = "GCA", ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
                  obs = obs$obs, obs_init = obs$obs_init, postproc_out = pp_out)
      if(compute_crps){
        crps_list_tmp[,,RR] <- crps_wrapper(gca$mvppout, obs$obs)
      }
      tmp <- eval_all_mult(mvpp_out = gca$mvppout, obs = obs$obs)
      es_list_tmp[,RR] <- tmp$es
      vs1_list_tmp[,RR] <- tmp$vs1
      vs1w_list_tmp[,RR] <- tmp$vs1w
      vs0_list_tmp[,RR] <- tmp$vs0
      vs0w_list_tmp[,RR] <- tmp$vs0w
      
      end_time <- Sys.time()
      
      timing_list_tmp[RR] <- end_time - start_time
    }
    crps_list$gca[rr,,] <- apply(crps_list_tmp, c(1,2), mean)
    es_list$gca[rr, ] <- apply(es_list_tmp, 1, mean)
    vs1_list$gca[rr, ] <- apply(vs1_list_tmp, 1, mean)
    vs1w_list$gca[rr, ] <- apply(vs1w_list_tmp, 1, mean)
    vs0_list$gca[rr, ] <- apply(vs0_list_tmp, 1, mean)
    vs0w_list$gca[rr, ] <- apply(vs0w_list_tmp, 1, mean)
    timing_list$gca[rr] <- apply(timing_list_tmp, 1, mean)
    
    # Archimedean copulas 
    for (method in c("Clayton","Frank", "Gumbel")) { #Excluded GOF from FLOS copulas
      es_list_tmp <- vs1_list_tmp <- vs1w_list_tmp <-
        vs0_list_tmp <- vs0w_list_tmp <- param_list_temp  <- matrix(NA, nrow = nout, ncol = rand_rep)
      crps_list_tmp <- array(NA, dim = c(nout, d, rand_rep))
      timing_list_tmp <- array(NA, dim = rand_rep)
      for(RR in 1:rand_rep){
        start_time <- Sys.time()
        mvd <- mvpp(method = method, ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
                    obs = obs$obs, obs_init = obs$obs_init, postproc_out = pp_out)
        #print(mvd)
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
     # print(" aller laatste check!! :)")
    }
    

    # Archimedean copulas shuffle
    mvppout_all <- list(Clayton = vector("list", rand_rep),
                        Frank   = vector("list", rand_rep),
                        Gumbel  = vector("list", rand_rep))
    for (method in c("Clayton","Frank", "Gumbel")) { 
      es_list_tmp <- vs1_list_tmp <- vs1w_list_tmp <-
        vs0_list_tmp <- vs0w_list_tmp <- param_list_temp  <- matrix(NA, nrow = nout, ncol = rand_rep)
      crps_list_tmp <- array(NA, dim = c(nout, d, rand_rep))
      timing_list_tmp <- array(NA, dim = rand_rep)
      for(RR in 1:rand_rep){
        start_time <- Sys.time()
        mvd <- mvpp(method = method, ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
                    obs = obs$obs, obs_init = obs$obs_init, postproc_out = pp_out)
        
        method_shuffle <- paste0(method, "_shuffle")
        COBASE <- mvpp(method = method_shuffle, ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
                    EMOS_sample = emos.q$mvppout, obs = obs$obs, obs_init = obs$obs_init, postproc_out = pp_out, MVPP_sample = mvd$mvpp)
       
        
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
  
    # end loop over Monte Carlo repetitions
  }
  
  # return results, as a list
  out <- list("crps_list" = crps_list, "es_list" = es_list, "vs1_list" = vs1_list,
              "vs1w_list" = vs1w_list, "vs0_list" = vs0_list, "vs0w_list" = vs0w_list, "param_list" = param_list,
               "timing_list" = timing_list)
  return(out)
}

#------------------------------------------------------------------
# parameters to run

input_copula <- c("Frank","Gumbel", "Clayton", "Frank_shuffle", "Gumbel_shuffle", "Clayton_shuffle")
#reproduction data from paper of lerch et al.:
input_d <- 3
repetitions <- 1
eps <- 1
sigma <- sqrt(5)
rho0 <- c(0.25, 0.5, 0.75)
rho <- c(0.25, 0.5, 0.75)
input_par <- expand.grid(input_copula, rho0, eps, sigma, rho, input_d, 1:repetitions)
names(input_par) <- c("copula", "rho0", "eps", "sigma", "rho", "d", "repetition")
MC_reps <- 100

evalDays <- 150
timeWindow <- 30
trainingDays <- 50
ensembleMembers <- 50
randomRepetitions <- 1
Rdata_dir <- "C:/Users/20202943/Documents/BAM BEP/CODE/Simulation code/simulation data/Rdata_dir simulation/" # directory to save Rdata files to
Rout_dir <- "C:/Users/20202943/Documents/BAM BEP/CODE/Simulation code/simulation data/Rout_dir simulation/"  # directory to save Rout files to
run_wrapper <- function(runID){
  sink(file = paste0(Rout_dir, "setting_bep_", runID, ".Rout"))
  # Frank copula can only deal with positive tau values
  tau <- runif(1, 0, 1)
  par_values <- as.numeric(input_par[ID, ])
  res <- run_setting6(obsmodel = 7, fcmodel = 7, nout = evalDays, ninit = trainingDays, 
                      nmembers = ensembleMembers,
                      MCrep = MC_reps, rand_rep = randomRepetitions, 
                      progress_ind = TRUE, compute_crps = TRUE,
                      rho0 = input_par$rho0[runID], 
                      rho = input_par$rho[runID], 
                      eps = input_par$eps[runID], 
                      sigma = input_par$sigma[runID],
                      copula = input_par$copula[runID], 
                      d = input_par$d[runID])
  savename <- paste0(Rdata_dir, "res_setting_bep_", runID, ".Rdata")
  save(res, input_par, file = savename)
  sink()
}


for (ID in 1:dim(input_par)[1]) {
  closeAllConnections()
  print(ID)
  run_wrapper(runID = ID)
}


