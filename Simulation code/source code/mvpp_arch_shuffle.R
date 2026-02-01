# multivariate post-processing methods, applied to result of univariate post-processing

# the following methods are applied:
#   - EMOS-X: not accounting for any multivariate dependencies, independent samples from fc distributions in margins
#       EMOS-Q: equidistant quantiles at levels 1/m+1, ..., m/m+1
#       EMOS-OQ: (CRPS-) optimal quantiles, at levels (1-0.5)/m, ..., (m-0.5)/m
#       EMOS-R: random sample
#       EMOS-S: stratifed sampling approach of Hu et al (2016)
#       EMOS-T: fit parametric forecast distribution (eg normal) to ensemble,
#               and use quantiles corresponding to ens forecasts as quantile levels
#               for draws from post-processed distribution
#   - ECC-X: inheriting the multivariate dependency structure from the ensemble forecast
#             draw samples (as in EMOS-X), and re-order
#       ECC-Q: equidistant quantiles at levels 1/m+1, ..., m/m+1
#       ECC-OQ: (CRPS-) optimal quantiles, at levels (1-0.5)/m, ..., (m-0.5)/m
#       ECC-R: random sample
#       ECC-S: stratifed sampling approach of Hu et al (2016)
#       ECC-T: see above
#   - dECC: dynamic ECC, ECC variant proposed by Bouallegue et al. (2016, MWR)
#       ...
#       there are also different variants here!
#   - GCA: Gaussian copula approach from Pinson and Girard (2012)
#       GCA: with covariance matrix estimated from ensemble forecasts
#   - Schaake shuffle: reordering based on historical observations
#       SSh-Q: equidistant quantiles at levels 1/m+1, ..., m/m+1
#       SSh-OQ: (CRPS-) optimal quantiles, at levels (1-0.5)/m, ..., (m-0.5)/m
#       SSh-R: random sample
#       SSh-S: stratifed sampling approach of Hu et al (2016)
#       SSh-T: see above

# input:
#   method: indicate reordering method: "none", "ECC", "dECC", "GCA", "SSh",
#   variant: variant of re-ordering method to be used:
#       Q: ...
#       R: ...
#   ensfc, ensfc_init, obs, obs_init: output of generators of observations and ensemble fcsts
#     dimensions need to match
#   postproc_out: output of postproc code in postproc_ensfc.R (array of EMOS parameter values)
#   EMOS_sample: give specific EMOS sample to the function to use in re-ordering method
#     to ensure that for the methods with randomness, the same EMOS samples are used
#     otherwise, obtaining a new EMOS-R sample as basis for ECC-R means that different
#     random samples are used in each margin
#     This input is not required, and only used in the run_all function for the simulation


psurv_norm <- function(x, mean = 0, sd = 1) {
  return(1 - pnorm(x, mean = mean, sd = sd))
}

qsurv_norm <- function(x, mean = 0, sd = 1) {
  return(1 - qnorm(x, mean = mean, sd = sd))
}

dsurv_norm <- function(x, mean = 0, sd = 1) {
  return(-dnorm(x, mean = mean, sd = sd))
}

# getCentralDays <- function(value, minRange, maxRange) {
#   DAYS_PER_YEAR <- 365
# 
#   # Find smallest value
#   smallest_value <- as.integer(value - DAYS_PER_YEAR * floor((value - minRange) / DAYS_PER_YEAR))
# 
#   # Find largest value
#   largest_value <- as.integer(value + DAYS_PER_YEAR * floor((maxRange - value) / DAYS_PER_YEAR))
# 
#   return(seq(smallest_value, largest_value, by = DAYS_PER_YEAR))
# }
# 
# getIntervals <- function(day, minRange, maxRange, interval_length) {
#   central_days <- getCentralDays(day, minRange, maxRange)
# 
#   intervals <- unlist(lapply(central_days, function(x) {
#     seq(max(minRange, x - interval_length), min(maxRange, x + interval_length))
#   }))
# 
#   return(intervals)
# }


library(copula)
mvpp <- function(method, variant = NULL, ensfc, ensfc_init, obs, obs_init, postproc_out, sim_matrix = NULL, EMOS_sample = NULL, ECC_out = NULL, uvpp = NULL, MVPP_sample = NULL) {
  set.seed(2025)

  # generate array for output
  params <- array(NA, dim = dim(ensfc)[1])
  n <- dim(ensfc)[1]
  if (method %in% c("ECC", "dECC")) {
    m <- dim(ensfc)[2]
  } else {
    m <- 50 #ecc_m
  }
  d <- dim(ensfc)[3]
  mvppout <- array(NA, dim = c(n, m, d))

  # EMOS methods: no accounting for multivariate dependencies, simply draw from marginals
  if (method == "EMOS") {
    if (!is.null(EMOS_sample)) {
      stop("for method = 'EMOS', input 'EMOS_sample' must be 'NULL'")
    }
    if (is.null(variant)) {
      stop("for method = 'EMOS', a variant parameter must be specified")
    }

    if (variant == "R") {
      for (nn in 1:n) {
        for (dd in 1:d) {
          par <- postproc_out[nn, dd, ]
          mvppout[nn, , dd] <- rnorm(m, mean = par[1], sd = par[2])
        }
      }
    } else if (variant == "Q") {
      qlevels <- 1:m / (m + 1)
     
      for (nn in 1:n) {
        for (dd in 1:d) {
          par <- postproc_out[nn, dd, ]
        
          mvppout[nn, , dd] <- qnorm(qlevels, mean = par[1], sd = par[2])
        }
      }
    } else if (variant == "QO") {
      qlevels <- (1:m - 0.5) / m
      for (nn in 1:n) {
        for (dd in 1:d) {
          par <- postproc_out[nn, dd, ]
          mvppout[nn, , dd] <- qnorm(qlevels, mean = par[1], sd = par[2])
        }
      }
    } else if (variant == "S") {
      breakpoints <- 0:m / m
      qlevels <- runif(m, min = breakpoints[1:m], max = breakpoints[2:(m + 1)])
      for (nn in 1:n) {
        for (dd in 1:d) {
          par <- postproc_out[nn, dd, ]
          mvppout[nn, , dd] <- qnorm(qlevels, mean = par[1], sd = par[2])
        }
      }
    } else if (variant == "T") {
      for (nn in 1:n) {
        for (dd in 1:d) {
          ensfc_tmp <- ensfc[nn, , dd]
          ens_par <- c(mean(ensfc_tmp), sd(ensfc_tmp))
          qlevels <- pnorm(ensfc_tmp, mean = ens_par[1], sd = ens_par[2])
          postproc_par <- postproc_out[nn, dd, ]
          mvppout[nn, , dd] <- qnorm(qlevels, mean = postproc_par[1], sd = postproc_par[2])
        }
      }
    }
    # end of EMOS methods
  }

  # ECC code
  if (method == "ECC") {
    # if no EMOS_sample to base ECC on is given, recursively call 'mvpp' to generate such a sample
    if (is.null(EMOS_sample)) {
      EMOS_sample <- mvpp(
        method = "EMOS", variant = variant, postproc_out = postproc_out,
        ensfc = ensfc, ensfc_init = ensfc_init,
        obs = obs, obs_init = obs_init
      )
      message("no 'EMOS_sample' given for ECC, a new one is generated")
    } else if (!is.null(variant)) {
      message("'variant' parameter has no influence if EMOS_sample is supplied,
              make sure the EMOS_sample is produced with the desired variant")
    }
    
    # application of ECC is independent of 'variant' parameter, dependency is only through EMOS_sample
    for (nn in 1:n) {
      for (dd in 1:d) {
        ensfc_tmp <- ensfc[nn, , dd]
        mvppout[nn, , dd] <- EMOS_sample[nn, , dd][rank(ensfc_tmp, ties.method = "random")]
      }
    }
    # end of ECC code
  }

  # Schaake shuffle code
  if (method == "SSh") { #SSh-H
    # if no EMOS_sample to base SSh on is given, recursively call 'mvpp' to generate such a sample
    if (is.null(EMOS_sample)) {
      if (is.null(variant)) {
        stop("if no 'EMOS_sample' is given, 'variant' has to be specified to generate a new one")
      }
      EMOS_sample <- mvpp(
        method = "EMOS", variant = variant, postproc_out = postproc_out,
        ensfc = ensfc, ensfc_init = ensfc_init,
        obs = obs, obs_init = obs_init, ecc_m = ecc_m
      )
      message("no 'EMOS_sample' given for SSh, a new one is generated")
    } else if (!is.null(variant)) {
      message("'variant' parameter has no influence if EMOS_sample is supplied,
              make sure the EMOS_sample is produced with the desired variant")
    }

    obs_all <- rbind(obs_init, obs)

    # reorder post-processed forecast sample according to past observations
    for (nn in 1:n) {
      obs_IDs <- sample(x = 1:(dim(obs_init)[1] + nn - 1), size = m, replace = FALSE)
      for (dd in 1:d) {
        obs_tmp <- obs_all[obs_IDs, dd]
        mvppout[nn, , dd] <- EMOS_sample[nn, , dd][rank(obs_tmp, ties.method = "random")]
      }
    }
    # end of SSh code
  }
  
  
  # GCA code
  if(method == "GCA"){
    require(MASS)
    
    # if no EMOS_sample to base GCA on is given, recursively call 'mvpp' to generate such a sample
    if(any(!is.null(c(EMOS_sample, variant)))){
      message("'EMOS_sample' and 'variant' input have no effect for GCA")
    } 
    
    # concatenate obs_init and obs arrays to determine covariance matrix for Gaussian copulas
    obs_all <- rbind(obs_init, obs)
    
    for(nn in 1:n){
      obs_train <- obs_all[1:(dim(obs_init)[1]+nn-1), ]
      # estimate covariance matrix
      cov_obs <- cov(obs_train)
      # draw random sample from multivariate normal distribution with this covariance matrix
      mvsample <- mvrnorm(n = m, mu = rep(0,d), Sigma = cov_obs)
      for(dd in 1:d){
        par <- postproc_out[nn, dd, ]
        mvppout[nn, , dd] <- qnorm(pnorm(mvsample[, dd]), mean = par[1], sd = par[2])
      }
    } 
    
    # end of GCA code   
  }


  # Archimedean copula code
  if (any(match(c(
    "Clayton", "Frank", "Gumbel"
  ), method, nomatch = 0) != 0)) {
    require(MASS)

    # if no EMOS_sample to base Clayton on is given, recursively call 'mvpp' to generate such a sample
    if (any(!is.null(c(EMOS_sample, variant)))) {
      message("'EMOS_sample' and 'variant' input have no effect for Archimedean copulas")
    }

    obs_all <- rbind(obs_init, obs)

    for (nn in 1:n) {
      # Use a shifting time window of length obs_init to compute copula parameters
      obs_train <- obs_all[1:(dim(obs_init)[1]+nn-1), ]

      d <- dim(obs_train)[2]
      # The input for the copulas are the CDF values of the margins
      obs_train_CDF <- c()
      mean_vector <- c()
      sd_vector <- c()
      for (dd in 1:d) {
        # Use EMOS values for marginals
        if (!is.null(uvpp)) {
          dat <- subset(uvpp, stat == dd)
          averagedMean <- unlist(unname(dat[nn, ]["ens_mu"]))
          averagedSd <- unlist(unname(dat[nn, ]["ens_sd"]))
        } else {
          par <- postproc_out[nn, dd, ]
          averagedMean <- par[1]
          averagedSd <- par[2]
        }

        obs_train_CDF <- cbind(obs_train_CDF, pnorm(obs_train[, dd], mean = averagedMean, sd = averagedSd))

        mean_vector <- c(mean_vector, averagedMean)
        sd_vector <- c(sd_vector, averagedSd)
      }




      # Estimate the parameter and copula - itau method does not converge for some cases (without giving a warning/ error and runs indefinitely)
      # Copula parameters are bounded to prevent generating inf samples
      if (method == "Clayton") {
        fitcop <- tryCatch(
          {
            fitCopula(claytonCopula(dim = d), data = obs_train_CDF, method = "mpl", start = 1, optim.control = list(maxit = 1000), upper = 100)
          },
          error = function(e) {
            indepCopula(dim = d)
          }
        )
      } else if (method == "Frank") {
        fitcop <- tryCatch(
          {
            fitCopula(frankCopula(dim = d), data = obs_train_CDF, method = "mpl", start = 1, optim.control = list(maxit = 1000), upper = 100)
          },
          error = function(e) {
            
            indepCopula(dim = d)
          }
        )
      } else if (method == "Gumbel") {
        fitcop <- tryCatch(
          {
            fitCopula(gumbelCopula(dim = d), data = obs_train_CDF, method = "mpl", start = 1, optim.control = list(maxit = 1000), upper = 100)
          },
          error = function(e) {
            indepCopula(dim = d)
          }
        )
      } else {
        stop("Incorrect copula")
      }

      if (!(class(fitcop) == "indepCopula")) {
        cop <- fitcop@copula

        # save the fitted parameter
        params[nn] <- fitcop@estimate
      } else {
        cop <- fitcop
      }

      if (any(match(c(
        "Clayton", "Frank", "Gumbel"), method, nomatch = 0) != 0)) {
        # draw random sample from multivariate normal distribution with this estimated copula
        paramMargins <- list()
        surv_paramMargins <- list()

        for (i in 1:d) {
          paramMargins[[i]] <- list(mean = mean_vector[i], sd = sd_vector[i])
          surv_paramMargins[[i]] <- list(mean = -mean_vector[i], sd = sd_vector[i])
        }

      mvDistribution <- mvdc(
            copula = cop, margins = rep("norm", d),
            paramMargins = paramMargins
          )
    


        max_reps <- 2
        rep <- 0
        repeat {
          rep <- rep + 1
          mvsample <- rMvdc(m, mvDistribution)

          if (all(is.finite(mvsample)) | rep >= max_reps) {
            break
          }

          set.seed(sample(1:1000, 1))
        }

        for (dd in 1:d) {
          mvppout[nn, , dd] <- mvsample[, dd]
        }
      } else {
        U <- rCopula(m, copula = cop)

        for (dd in 1:d) {
          obs_tmp <- U[, dd]
          mvppout[nn, , dd] <- EMOS_sample[nn, , dd][rank(obs_tmp, ties.method = "random")]
        }
      }
    }

    # end of Archimedean Copula code
  }

    if (any(match(c(
      "Clayton_shuffle", "Gumbel_shuffle", "Frank_shuffle"), method, nomatch = 0) != 0)){
      if (is.null(EMOS_sample)) {
        EMOS_sample <- mvpp(
          method = "EMOS", variant = "Q", postproc_out = postproc_out,
          ensfc = ensfc, ensfc_init = ensfc_init,
          obs = obs, obs_init = obs_init
        )
        message("no 'EMOS_sample' given for shuffle, new one is generated, requires 'variant'")
      }
      EMOS_sample1 <-EMOS_sample
      mvppout <- array(NA, dim = dim(EMOS_sample1))
    
    for (nn in 1:n)
    {
      # impose dependence structure on post-processed forecasts
      for (dd in 1:d)
      {
        obs_tmp <- MVPP_sample[nn, , dd]
        mvppout[nn, , dd] <- EMOS_sample1[nn, , dd][rank(obs_tmp, ties.method = "random")]
      }
    }
    
    return(list("mvppout" = mvppout))
  }
  # dECC
  if (method == "dECC") {
    
    if (is.null(EMOS_sample)) {
      EMOS_sample <- mvpp(
        method = "EMOS", variant = variant, postproc_out = postproc_out, 
        ensfc = ensfc, ensfc_init = ensfc_init,
        obs = obs, obs_init = obs_init
      )
      message("no 'EMOS_sample' given for dECC, new one is generated, requires 'variant'")
    }

    if (is.null(ECC_out)) {
      message("no 'ECC_out' given, generating new one from input EMOS_sample")
      if (!is.null(variant)) {
        message("input 'variant' has no effect here,
                                    make sure variant in EMOS_sample is desired one")
      }
      ECC_out <- mvpp(
        method = "ECC", postproc_out = postproc_out,
        ensfc = ensfc, ensfc_init = ensfc_init,
        obs = obs, obs_init = obs_init, EMOS_sample = EMOS_sample
      )
    }

    
    # e in the dECC paper, eq (14)
    fcerror <- obs_init - apply(ensfc_init, c(1, 3), mean)
    # compute error correlation matrix and root
    R_e <- cor(fcerror)
    eig <- eigen(R_e)
    Re_root <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)

    # c in the dECC paper, eq (18)
    error_cor <- ECC_out - ensfc

    # breve c, eq (19)
    breve_c <- array(NA, dim = dim(ensfc))
    for (nn in 1:n) {
      for (mm in 1:m) {
        breve_c[nn, mm, ] <- Re_root %*% error_cor[nn, mm, ]
      }
    }

    # adjusted ensemble breve x, eq (22)
    breve_x <- ensfc + breve_c

    # apply ECC with breve_x as dependence template (Step 5 in paper)
    for (nn in 1:n) {
      for (dd in 1:d) {
        ensfc_tmp <- breve_x[nn, , dd]
        mvppout[nn, , dd] <- EMOS_sample[nn, , dd][rank(ensfc_tmp, ties.method = "random")]
      }
    }

    # end of dECC code
  }

  return(list("mvppout" = mvppout, "params" = params))

}