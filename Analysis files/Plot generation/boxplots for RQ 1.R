
#libraries

library(dplyr)
library(tidyr)
library(ggplot2)

#set benchmark
benchmark_model <- "ecc.q"
# benchmark_model <- "ssh"
#benchmark_model <- "gca"

#input data for benchmark
#load(paste0("C:/Users/20202943/Documents/BAM BEP/setting 6/goed/data/_setting_6_obsmodel_6_fcmodel_6.Rdata"))
#load(paste0("C:/Users/20202943/Documents/BAM BEP/setting 6/goed/data/_withmeandiffwith0_6_obsmodel_6_fcmodel_6.Rdata")
load( paste0("C:/Users/20202943/Documents/BAM BEP/setting 6/goed/data/","_settingpval2_SSH_",6, "_obsmodel_",6,"_fcmodel_",6,".Rdata"))
load(paste0("C:/Users/20202943/Documents/BAM BEP/setting 6/goed/data/","_withmeandiff_ssh_",6, "_obsmodel_",6,"_fcmodel_",6,".Rdata"))


#figure settings
plotWidth <- 12
plotHeight <- 6 
plot_folder <- paste0("C:/Users/20202943/Documents/BAM BEP/setting 6/goed/data/")
dir.create(file.path(plot_folder), showWarnings = FALSE)
save_figure <- function(fileName, fig) {
  res <- 250
  ggsave(paste0(plot_folder, fileName),
    fig,
    width = plotWidth,
    height = plotHeight,
    dpi = res,
    limitsize = FALSE
  )
}

#adapt df
dfmc2 <- dfmc %>%
  filter(!score %in% c("param_list", "chosenCopula_list", "timing_list"))
dfmc2_copy <- dfmc2

dfmc2_copy$value <- (-1)*dfmc2_copy$value
input_scores <- unique(dfmc2_copy$score)

# drop model ensemble, emos and crps scores
dfmc2_copy_save <- dfmc2_copy 
dfmc2_copy <- subset(dfmc2_copy_save, model != "ens")
df2 <- subset(dfmc2_copy, model != "emos.q")
df2_rightscores <- subset(df2, score != "crps_list")

Plot_scores <- function(inputplot, thismodel){
  
  alpha <- 0.25
  
  #filter model
  dfm <- subset(inputplot, model == thismodel)
  
  quants <- unname(quantile(dfm$value, c(0.01, 0.99), na.rm = TRUE))
  ylimits <- c(
    1.5 * min(quants[1], qnorm(alpha)),
    1.5 * max(quants[2], qnorm(1 - alpha))
  )
  
  score_vec <- c(
    es_list   = "ES",
    vs0_list  = "VS (p=0.5)",
    vs1_list  = "VS (p=1)",
    vs0w_list = "VS0w",
    vs1w_list = "VS1w"
  )
  
  p1 <- ggplot(dfm, aes(x = score, y = value, colour = score)) +
    geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=qnorm(alpha), ymax=qnorm(1-alpha)),
              inherit.aes = FALSE,
              fill = "gray75", color="gray75", alpha=alpha) +
    geom_boxplot(outlier.shape = NA, na.rm = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray25") +
    facet_grid(rows = vars(rho), cols = vars(rho0),
               labeller = label_bquote(rows = rho==.(rho),
                                       cols = rho[0] == .(rho0))) +
    coord_cartesian(ylim = ylimits) +
    theme_bw() +
    theme(legend.position = "bottom") +
    xlab("Score") + ylab("DM test statistic") +
    ggtitle(paste0("Model: ", thismodel, " vs benchmark model ", benchmark_model)) +
    scale_x_discrete(labels = score_vec)
  
  return(p1)
}
#generate all plots for all models
input_models <- c("ecc.q","ecc.s","decc.q","gca", "ssh","Clayton","Frank","Gumbel","Clayton_shuffle", "Frank_shuffle", "Gumbel_shuffle" ) #gca and clayton were NA values
for(m in input_models){ 
  test <- Plot_scores(df2_rightscores, m)
  save_figure(paste0("boxplot_DM_per_scoretest_gca_",m,".png"),test)
}

#plots for artificial data
#load DM statistics
load("C:/Users/20202943/Documents/BAM BEP/setting 6/goed/data/artificial_univariate_DM.Rdata")
dm_to_long <- function(df_wide) {
  df_wide %>%
    pivot_longer(
      cols = matches("^bootstrap_\\d+$|^full$"),
      names_to  = "replicate",
      values_to = "value"
    )
}

inputplot_plain <- dm_to_long(dfmc_bootstrapped_artificial) 
#again drop ensemble, emos.q and crps
dfmc2_copy <- subset(inputplot_plain, model != "ens")
dfmc2_copy$value <- (-1)*dfmc2_copy$value
df2 <- subset(dfmc2_copy, model != "emos.q")
df2_rightscores <- subset(df2, score != "crps_list")

Plot_scores_artificial <- function(inputplot, thismodel, benchmark_model = NA) {
  
  alpha <- 0.25
  
  ## Filter to one model
  one_model <- subset(inputplot, model == thismodel)
  
  quants <- unname(quantile(one_model$value, c(0.01, 0.99), na.rm = TRUE))
  ylimits <- c(
    1.5 * min(quants[1], qnorm(alpha)),
    1.5 * max(quants[2], qnorm(1 - alpha))
  )
  
  score_vec <- c(
    es_list   = "ES",
    vs0_list  = "VS (p=0.5)",
    vs1_list  = "VS (p=1)",
    vs0w_list = "VS0w",
    vs1w_list = "VS1w"
  )
  p1 <- ggplot(one_model, aes(x = score, y = value, colour = score)) +
    geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=qnorm(alpha), ymax=qnorm(1-alpha)),
              inherit.aes = FALSE,
              fill = "gray75", color="gray75", alpha=alpha) +
    geom_boxplot(outlier.shape = NA, na.rm = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray25") +
    coord_cartesian(ylim = ylimits) +
    theme_bw() +
    theme(legend.position = "bottom") +
    xlab("Score") + ylab("DM test statistic") +
    ggtitle(paste0("Model: ", thismodel, " vs benchmark model ", benchmark_model)) +
    scale_x_discrete(labels = score_vec)
  
  return(p1)
}

#generate all plots for all models
input_models <- c("ecc.q","ecc.s","decc.q","ssh", "gca","Clayton","Frank","Gumbel","Clayton_shuffle", "Frank_shuffle", "Gumbel_shuffle" )
for(m in input_models){ 
  test <- Plot_scores_artificial(df2_rightscores, thismodel =m,benchmark_model = "ssh")
  save_figure(paste0("boxplot_DM_per_scoretest_testmock_ssh",m,".png"),test)
}
