# code to compute test statistics of DM tests from Flos et al.

rm(list=ls())
#libraries
library(forecast) 
#settings
input_copula <- c("Frank","Gumbel", "Clayton", "Frank_shuffle", "Gumbel_shuffle", "Clayton_shuffle")
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
fName <- paste0("res_setting6_")

df_raw <- data.frame(input_par)
df_raw$file <- 1:nrow(df_raw)

flist <- list.files("C:/Users/20202943/Documents/BAM BEP/setting 6/goed/data/")
existing <- as.numeric(sapply(flist, FUN = function(x) as.numeric(strsplit(strsplit(x, fName)[[1]][2], ".Rdata"))))
existing

df <- df_raw[which(is.element(df_raw$file, existing)),] 
print(df)
load(paste0("C:/Users/20202943/Documents/BAM BEP/setting 6/goed/data/",fName, "1.Rdata"))

input_models <- names(res$es_list)
input_models <-setdiff(names(res$es_list), "GOF")
input_scores <- names(res)
#dataframe creation
df_use <- as.data.frame(df[1,])
df_use$model <- as.character("a")
df_use$score <- as.character("a")

model_score_grid <- expand.grid(input_models, input_scores[!input_scores %in% c("param_list", "indep_list", "tau", "timing_list")])
model_score_grid <- expand.grid(input_models, input_scores)
print(model_score_grid)
for(i in 1:nrow(df)){
  df_use[((i-1)*nrow(model_score_grid)+1):(i*nrow(model_score_grid)),] <- df[i,]
  df_use[((i-1)*nrow(model_score_grid)+1):(i*nrow(model_score_grid)),]$model <- as.character(model_score_grid$Var1)
  df_use[((i-1)*nrow(model_score_grid)+1):(i*nrow(model_score_grid)),]$score <- as.character(model_score_grid$Var2)
}

dfmc <- data.frame(cbind(zoo::coredata(df_use)[rep(seq(nrow(df_use)),MC_reps),]))
dfmc$value <- NA
dfmc$p_value <- NA
print(dfmc)
#DM tests
for(ID in existing[!is.na(existing)]){
  load(paste0("C:/Users/20202943/Documents/BAM BEP/setting 6/goed/data/", fName, ID,".Rdata"))
  print(ID)
  for(this_model in input_models){
    for(this_score in input_scores){
      ind <- which(dfmc$file == ID & dfmc$model == this_model & dfmc$score == this_score)  
      
      if(this_score == "crps_list"){
        dm_teststat_vec <- rep(NA, MC_reps)
        dm_pvalue_vec <- rep(NA, MC_reps)
        for(MC_rep in 1:MC_reps){
          tmp <- NA
          tmp_p <- NA
          tryDM <- try(tmp_DM <- dm.test(e1 = res[[which(input_scores == this_score)]][[which(input_models == this_model)]][,,1][MC_rep,],
                                         e2 = res[[which(input_scores == this_score)]][[which(input_models == "ecc.q")]][,,1][MC_rep,], 
                                         #e2 = res[[which(input_scores == this_score)]][[which(input_models == "gca")]][,,1][MC_rep,], 
                                         #e2 = res[[which(input_scores == this_score)]][[which(input_models == "ssh")]][,,1][MC_rep,], 
                                         h = 1, power = 1), silent = TRUE)
          if(class(tryDM) != "try-error"){
            tmp <- tmp_DM$statistic
            tmp_p <- tmp_DM$p.value
          } else{
            tmp <- 0
            tmp_p <- 0
          }
          dm_teststat_vec[MC_rep] <- tmp
          dm_pvalue_vec[MC_rep] <- tmp_p
        }
      } else{
        dm_teststat_vec <- rep(NA, MC_reps)
        for(MC_rep in 1:MC_reps){
          tmp <- NA
          tmp_p <- NA
          tryDM <- try(tmp_DM <- dm.test(e1 = res[[which(input_scores == this_score)]][[which(input_models == this_model)]][MC_rep,],
                                         e2 = res[[which(input_scores == this_score)]][[which(input_models == "ecc.q")]][MC_rep,],
                                         #e2 = res[[which(input_scores == this_score)]][[which(input_models == "gca")]][,,1][MC_rep,], 
                                         #e2 = res[[which(input_scores == this_score)]][[which(input_models == "ssh")]][,,1][MC_rep,],
                                         h = 1, power = 1), silent = TRUE)
          if(class(tryDM) != "try-error"){
            tmp <- tmp_DM$statistic
            tmp_p <- tmp_DM$p.value
          } else{
            tmp <- 0
            tmp_p <-0
          }
          dm_teststat_vec[MC_rep] <- tmp
          dm_pvalue_vec[MC_rep] <- tmp_p
        }
      }
      
      
      dfmc$value[ind] <- dm_teststat_vec
      dfmc$p_value[ind] <- dm_pvalue_vec
    }
  }
  
}


save(dfmc, file = paste0("C:/Users/20202943/Documents/BAM BEP/setting 6/goed/data/","_settingpval2_ECC_Q_",6, "_obsmodel_",6,"_fcmodel_",6,".Rdata"))
#save(dfmc, file = paste0("C:/Users/20202943/Documents/BAM BEP/setting 6/goed/data/","_settingpval2GCA_",6, "_obsmodel_",6,"_fcmodel_",6,".Rdata"))
#save(dfmc, file = paste0("C:/Users/20202943/Documents/BAM BEP/setting 6/goed/data/","_settingpval2_SSH_",6, "_obsmodel_",6,"_fcmodel_",6,".Rdata"))


