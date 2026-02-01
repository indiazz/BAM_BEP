# generate multivariate ensemble forecasts 
# Note that parts of the notation may differ from the notation in the paper
#   ... in particular regarding the numbering of the settings which was changed during the revision (2 -> S1; 3 -> 2; 4 -> 3A, new: 3B)

# input:
#   model: character string; indicating which model is used for generating the forecasts
#   nout: number of multivariate observations to be generated as evaluation period
#   ninit: additional initial training period for model estimation purposes
#   nmembers: number of ensemble members
#   d: dimension of the multivariate vectors
#   ... additional parameters, depending on the chosen model

# output:
#   list of two arrays "ensfc_init" and "ensfc" containing the ensemble forecasts
#   dimensions: ninit, nmembers, d; and nout, nmembers, d
#   content in array dimensions:
#     first dimension: forecast instance 
#     second dimension: ensemble member
#     third dimension: dimension in multivariate setting


generate_ensfc <- function(model, nout, ninit, nmembers, d, ...){
  require(MASS)
  
  # check input 
  if(any(!is.numeric(c(nout, ninit, nmembers, d)))){
    stop("Input 'nout', 'ninit', 'nmembers' and 'd' need to be numeric of length 1")
  }
  m <- nmembers
  
  # initialize output arrays
  ensfc_init <- array(NA, dim = c(ninit, m, d))
  ensfc <- array(NA, dim = c(nout, m, d))
  
  #Recreation of Lerch et al. ensemble forecasts
  if (model == 7){
    
    # check if appropriate additional parameters are given
    input <- list(...)
    required <- c("eps", "sigma", "rho")
    ind_check <- match(required, names(input), nomatch = 0)
    if(any(ind_check == 0)){
      stop(paste("Missing additional model-specific parameter",
                 paste("Given input:", paste(names(input), collapse=", ")),
                 paste("Required input:", paste(required, collapse=", ")),
                 sep="\n")
      )
    }
    
    # assign model-specific parameters from input
    eps <- input$eps
    sigma <- input$sigma
    rho <- input$rho
    
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
  }

  return(list("ensfc_init" = ensfc_init, "ensfc" = ensfc))
}

#generation of ensemble forecasts for the artificial data as done in FLOS COBASE paper
generate_ensfc2 <- function(data, file_name, observation_columns, ensemble_regex, store_uvpp = TRUE, ...){
  load_dimension_transform <- function(stations, observation_columns) {
    get_dimension <- function(stat, obs) {
      stat_idx <- match(stat, stations)
      obs_idx <- match(obs, observation_columns)
      
      return(stat_idx + length(stations) * (obs_idx - 1))
    }
    
    return(get_dimension)
  }
  print("running transform data...")
  # Stations and days from the data
  stations <- unique(data$station)
  days <- sort(unique(data$td))
  trainingDays <- 60
  nout <- length(days) - trainingDays
  d <- length(stations) * length(observation_columns)
  get_dimension <- load_dimension_transform(stations, observation_columns)
  
  # Ensemble members
  m <- sum(grepl(ensemble_regex[1], names(data))) # Each observation should have the same number of ens members
  
  
  cat("Formatting ensemble and observations...\n")
  
  # Create ensfc and obs data structure
  ensfc_init <- array(NA, dim = c(trainingDays, m, d))
  ensfc <- array(NA, dim = c(nout, m, d))
  
  
  # Dimensions consists of groups of (station, obs_type)
  for (var_idx in 1:length(observation_columns))
  {
    # Column names for the ensemble members
    ensemble_members <- grepl(ensemble_regex[var_idx], names(data))
    
    for (station_nr in stations)
    {
      # The dimension of the variable
      dd <- get_dimension(station_nr, obs_column)
      
      
      for (day in 1:length(days))
      {
        # Extract forecast for day and dim = index
        dat <- subset(data, td == days[day] & station == station_nr)
        if (day <= trainingDays) {
          ensfc_init[day, , dd] <- unlist(dat[ensemble_members], use.names = FALSE)
          
        } else {
          ensfc[day - trainingDays, , dd] <- unlist(dat[ensemble_members], use.names = FALSE)
          
        }
      }
    }
  }
  return(list(
    "ensfc"       = ensfc,
    "ensfc_init"  = ensfc_init,
    m))
}
# 
# fc <- generate_ensfc2( data = data_mock,
#                       file_name               = "Mock_data",
#                       observation_columns     = c("obs"),
#                       ensemble_regex          = c("^M_"),
#                       store_uvpp = FALSE,
#                       output_dim_standard     = 10)