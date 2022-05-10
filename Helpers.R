# Simulate the data
simulate_ZIDM = function(                              n_obs = 50,
                                                       n_vars = 50,
                                                       n_vars_theta = 50,
                                                       n_taxa = 100, # Changed to 5 for testing
                                                       n_relevant_vars = 4,
                                                       n_relevant_taxa = 4,
                                                       n_relevant_taxa_theta = 4, 
                                                       n_relevant_vars_theta = 4, 
                                                       beta_min = 1.25,
                                                       beta_max = 1.5,
                                                       beta_min_theta = 1.25, # Changed to zeros for testing
                                                       beta_max_theta = 1.5, # Changed to zeros for testing 
                                                       int_zero_min = 0,
                                                       int_zero_max = 1,
                                                       signoise = 1.0,
                                                       n_reads_min = 1000,
                                                       n_reads_max = 2000,
                                                       theta0 = 0.01,
                                                       rho = NULL, 
                                                       Sigma = NULL,
                                                       rho_theta = NULL, 
                                                       Sigma_theta = NULL, 
                                                       seed = 1 ){
  
  set.seed( seed )
  
  # check for required packages
  if(!require(dirmult)){
    stop("dirmult package required")
  }
  if(!require(MASS)){
    stop("MASS package required")
  }
  if(!require(matrixcalc)){
    stop("matrixcalc package required")
  }
  
  # Defense
  if( !is.null( rho ) & !is.null( Sigma ) ){
    stop("Bad input: Please provide either rho or Sigma for covariate correlation structure.")
  }
  
  if( is.null( rho_theta ) & is.null( Sigma_theta ) ){
    stop("Bad input: Please provide either rho or Sigma for covariate correlation structure.")
  }
  
  if( !is.null( rho ) ){
    if( rho > 1 | rho < 0 ){
      stop("Bad input: Please provide rho between 0 and 1.")
    }
  }
  
  if( !is.null( rho_theta ) ){
    if( rho_theta > 1 | rho_theta < 0 ){
      stop("Bad input: Please provide rho between 0 and 1.")
    }
  }
  
  if( !is.null( Sigma ) ){
    if( !is.positive.definite( Sigma ) ){
      stop("Bad input: Please provide positive definite covariance matrix.")
    }
  }
  
  if( !is.null( Sigma_theta ) ){
    if( !is.positive.definite( Sigma_theta ) ){
      stop("Bad input: Please provide positive definite covariance matrix.")
    }
  }
  
  if( !is.null( Sigma ) ){
    if( ncol(Sigma) != n_vars ){
      stop("Bad input: Please provide covariance matrix to match the number of covariates")
    }
  }  
  
  if( !is.null( Sigma_theta ) ){
    if( ncol(Sigma_theta) != n_vars_theta ){
      stop("Bad input: Please provide covariance matrix to match the number of covariates")
    }
  }  
  
  # covariance matrix for predictors
  if( !is.null( rho ) ){
    Sigma <- matrix( 0, n_vars, n_vars )
    Sigma = rho^abs(row(Sigma) - col(Sigma))
  }
  
  if( !is.null( rho_theta ) ){
    Sigma_theta <- matrix( 0, n_vars_theta, n_vars_theta )
    Sigma_theta = rho_theta^abs(row(Sigma_theta) - col(Sigma_theta))
  }
  
  # include the intercept
  XX = cbind(rep(1, n_obs),
             scale(MASS::mvrnorm(n = n_obs, mu = rep(0, n_vars), Sigma = Sigma)))
  
  # Do the same for the x_thetas 
  XX_theta = cbind(rep(1, n_obs),
                   scale(MASS::mvrnorm(n = n_obs, mu = rep(0, n_vars_theta), Sigma = Sigma_theta)))
  
  # empties
  ZZ = matrix(0, n_obs, n_taxa)
  betas = matrix(0, n_taxa, n_vars)
  beta_0 = matrix(0, n_obs, n_taxa)
  
  betas_theta = matrix(0, n_taxa, n_vars_theta) 
  
  # parameters with signs alternating
  if( n_relevant_taxa != 0 & n_relevant_vars != 0 ){
    st = 0
    low_side = beta_min
    high_side = beta_max
    if( n_relevant_taxa != 1){
      # warning if the lengths don't match
      coef = suppressWarnings(seq(low_side, high_side, len = n_relevant_taxa) * c(1, -1))
    }else{
      coef = (low_side + high_side) / 2
    }
    coef_g = rep(1.0, len = n_relevant_vars)
    for(ii in 1:n_relevant_vars){
      # overlap species
      betas[(st:(st + n_relevant_taxa - 1)) %% n_taxa + 1, 3 * ii - 2] = coef_g[ii] * sample(coef)[((ii - 1):(ii + n_relevant_taxa - 2)) %% n_relevant_taxa + 1]
      st = st + 1
    }
  }
  
  # -2.3 and 2.3 so that the intercept varies over three orders of magnitude
  intercept <-  runif(n_taxa, -2.3, 2.3)  # seq( -2.3, 2.3, length.out = n_taxa )  
  Beta <- cbind(intercept, signoise * betas)
  
  # Get probability of etas
  # parameters with signs alternating
  if( n_relevant_taxa_theta != 0 & n_relevant_vars_theta != 0 ){
    st = 0
    low_side = beta_min_theta
    high_side = beta_max_theta
    if( n_relevant_taxa_theta != 1){
      # warning if the lengths don't match
      coef = suppressWarnings(seq(low_side, high_side, len = n_relevant_taxa_theta) * c(1, -1))
    }else{
      coef = (low_side + high_side) / 2
    }
    coef_g = rep(1.0, len = n_relevant_vars)  
    for(ii in 1:n_relevant_vars_theta){
      # overlap species
      betas_theta[(st:(st + n_relevant_taxa_theta - 1)) %% n_taxa + 1, 3 * ii - 2] = coef_g[ii] * sample(coef)[((ii - 1):(ii + n_relevant_taxa_theta - 2)) %% n_relevant_taxa_theta + 1]
      st = st + 1
    }
  }
  
  # Coefficients for beta_theta
  intercept_theta <- runif(n_taxa, int_zero_min, int_zero_max) #  seq( int_zero_min, int_zero_max, length.out = n_taxa )# runif(n_taxa, -2.3, 2.3)    #
  Beta_theta <- cbind( intercept_theta, signoise * betas_theta )
  prob <- matrix(0, n_obs, n_taxa)
  
  # Generate eta_ij (indicator of at-risk element j for subject i )
  eta <- matrix( 0, n_obs, n_taxa )  
  
  for(ii in 1:n_obs){
    for(jj in 1: n_taxa ){
      thisrow = as.vector(exp(Beta_theta[ jj, ] %*% XX_theta[ii, ]))
      prob[ ii, jj ] = thisrow/(1 + thisrow)
      eta[ ii, jj ] <-  rbinom( 1, 1, prob[ ii, jj ] ) 
    }
  }
  
  # row totals 
  ct0 = sample(n_reads_min:n_reads_max, n_obs, rep = T)
  for(ii in 1:n_obs){
    thisrow = as.vector(exp(Beta %*% XX[ii, ]))*eta[ ii, ]
    beta_0[ii, ] = thisrow/sum(thisrow)
    ZZ[ii, ] = dirmult::simPop(J = 1, n = ct0[ii], pi = beta_0[ii, ], theta = theta0)$data[1, ]
  }
  
  return(list(X = XX, X_theta = XX_theta, Z = ZZ, alphas = intercept, betas = Beta, betas_theta = Beta_theta,
              n_reads_min = n_reads_min, n_reads_max = n_reads_max, eta = eta, prob = prob,
              theta0 = theta0, beta_0 = beta_0, Beta = Beta, rho = rho, signoise = signoise, Sigma = Sigma))
  
}

ZIDM_R = function( 
  Z = NULL, 
  iterations = 10000, 
  thin = 10, 
  sigma2_beta_gamma = sqrt( 5 ), 
  sigma2_beta_theta = sqrt( 5 ),
  seed = 1
){
  
  set.seed( seed )
  
  # Simulate data 
  if( is.null( Z ) ){
    stop("Bad input: Data must be supplied")
  }
  
  if( sigma2_beta_gamma < 0 | sigma2_beta_theta < 0 ){
    stop("Bad input: Variances must be positive")
  }
   
  # Get dimensions 
  n_obs <- nrow( Z )
  n_taxa <- ncol( Z )
  samples <- floor( iterations/ thin )    
  
  # Set intercept term
  X <- matrix( 1, nrow = n_obs, ncol = 1 )
  n_cov  <- ncol( X ) 
  
  # Initialize the model  
  varphi. <- array( 0, dim = c( n_taxa, ncol( X ) , samples ) )
  beta_gamma. <- array( 0, dim = c( n_taxa, ncol( X )  , samples ) )
  
  # Adjust initial values for  beta_gamma, varphi,  
  varphi.[ , 1, ]  <- 1
  beta_gamma.[ , , 1]  <-  rnorm( n = n_taxa*ncol( X )  )*varphi.[ , ,1]  
   
  # initialize latent variables
  cc. <- array( 0, dim = c( nrow( Z ), n_taxa, samples ) )  
  cc.[ , , 1 ] <-  Z 
  uu <- rep( 0 , nrow( Z ) ) 
  for( n in 1:nrow( Z ) ){
    sum_Z <- sum( Z[ n, ] )
    uu[ n ] <- rgamma( 1, sum_Z, sum_Z );
  } 
  
  # Memory for ZI variable eta
  eta. <- array( 0, dim = c( nrow( Z ), n_taxa, samples ) )  
  
  # Initialize eta at 1
  eta.[,,1] <- matrix( 1, nrow = nrow( Z ), ncol = ncol( Z ) )  
  
  # Get index list for eta 
  eta_accept. <- array( 0, dim = c( nrow( Z ), n_taxa, samples ) ) 
  
  eta_index <- list( ) 
  for( b in 1:n_taxa ){
    eta_index[[b]]  <- which( Z[ , b ]  == 0 ) - 1
    eta.[ which( Z[ ,b ]  == 0 ), b , 1] <- rbinom( 1, 1, 0.5 )
  }
 
  beta_theta. <-  array( 0, dim = c( n_taxa, ncol(X)  , samples ) )  
  zeta. <-  array( 0, dim = c( n_taxa, ncol(X)  , samples ) )  
  
  zeta.[ , 1, ] <- 1
  beta_theta.[ , , 1]  <- rnorm( n = n_taxa*ncol(X)  )*zeta.[ , , 1] 
  
  omega. <-  array( 0, dim = c( nrow( Z ), n_taxa, samples ) )
  omega.[ , , 1 ] <- 1  
  
  MCMC <-  zidm( iterations, thin, Z, X, X, beta_gamma., cc., uu, eta., eta_accept., eta_index, omega., beta_theta., zeta., varphi., sigma2_beta_gamma, sigma2_beta_theta )
  out <<- MCMC[1:8]
  names( out ) <- c( "varphi", "beta_gamma", "eta", "beta_theta", "zeta", "omega", "eta_accept", "cc")
  
return( out )
}
  
ZIDMbvs_R = function( Z = NULL, 
                      X = NULL, 
                      X_theta = NULL, 
                      iterations = 10000, 
                      thin = 10, 
                      sigma2_beta_gamma = sqrt( 5 ), 
                      sigma2_beta_theta = sqrt( 5 ),
                      a_varphi = 1,
                      b_varphi = 1,
                      a_zeta = 1,
                      b_zeta = 1,
                      seed = 1 ){
  
  set.seed( seed )
  
  # Simulate data 
  if( is.null( Z ) | is.null( X ) | is.null( X_theta ) ){
    stop("Bad input: Data must be supplied")
  }
  
  if( sigma2_beta_gamma < 0 | sigma2_beta_theta < 0 ){
    stop("Bad input: Variances must be positive")
  }
  
  # Get dimensions 
  n_obs <- nrow( Z )
  n_taxa <- ncol( Z )
  samples <- floor( iterations/ thin )    
  
  # Set intercept term
  X <- cbind( 1, X )
  X_theta <- cbind( 1, X_theta )
  n_cov  <- ncol( X ) 
  n_cov_theta <- ncol( X_theta )
  
  # Initialize the model  
  varphi. <- array( 0, dim = c( n_taxa, ncol( X ) , samples ) )
  beta_gamma. <- array( 0, dim = c( n_taxa, ncol( X )  , samples ) )
  
  # Adjust initial values for  beta_gamma, varphi,  
  varphi.[ , 1, ]  <- 1
  beta_gamma.[ , , 1]  <-  rnorm( n = n_taxa*ncol( X )  )*varphi.[ , ,1]  
  
  # initialize latent variables
  cc. <- array( 0, dim = c( nrow( Z ), n_taxa, samples ) )  
  cc.[ , , 1 ] <-  Z 
  uu <- rep( 0 , nrow( Z ) ) 
  for( n in 1:nrow( Z ) ){
    sum_Z <- sum( Z[ n, ] )
    uu[ n ] <- rgamma( 1, sum_Z, sum_Z );
  } 
  
  # Memory for ZI variable eta
  eta. <- array( 0, dim = c( nrow( Z ), n_taxa, samples ) )  
  
  # Initialize eta at 1
  eta.[,,1] <- matrix( 1, nrow = nrow( Z ), ncol = ncol( Z ) )  
  
  # Get index list for eta 
  eta_accept. <- array( 0, dim = c( nrow( Z ), n_taxa, samples ) ) 
  
  eta_index <- list( ) 
  for( b in 1:n_taxa ){
    eta_index[[b]]  <- which( Z[ , b ]  == 0 ) - 1
    eta.[ which( Z[ ,b ]  == 0 ), b , 1] <- rbinom( 1, 1, 0.5 )
  }
  
  beta_theta. <-  array( 0, dim = c( n_taxa, ncol(X_theta)  , samples ) )  
  zeta. <-  array( 0, dim = c( n_taxa, ncol(X_theta)  , samples ) )  
  
  zeta.[ , 1, ] <- 1
  beta_theta.[ , , 1]  <- rnorm( n = n_taxa*ncol(X_theta)  )*zeta.[ , , 1] 
  
  omega. <-  array( 0, dim = c( nrow( Z ), n_taxa, samples ) )
  omega.[ , , 1 ] <- 1  
  
  
  MCMC <- zidm_bvs( iterations, thin, Z, X, X_theta, beta_gamma., cc., uu, eta., eta_accept., eta_index, omega., beta_theta., zeta., varphi., sigma2_beta_gamma, a_varphi, b_varphi, a_zeta , b_zeta,  sigma2_beta_theta )
  out <<- MCMC[1:8]
  names( out ) <- c( "varphi", "beta_gamma", "eta", "beta_theta", "zeta", "omega", "eta_accept", "cc")
  
  return(out)
}
   




DM_R = function(  Z = NULL, 
                  iterations = 10000, 
                  thin = 10, 
                  sigma2_beta_gamma = sqrt( 5 ),  
                  seed = 1
){
  
  set.seed( seed )
  
  # Simulate data 
  if( is.null( Z ) ){
    stop("Bad input: Data must be supplied")
  }
  
  if( sigma2_beta_gamma < 0  ){
    stop("Bad input: Variance must be positive")
  }
  
  # Get dimensions 
  n_obs <- nrow( Z )
  n_taxa <- ncol( Z )
  samples <- floor( iterations/ thin )    
  
  # Set intercept term
  X <- matrix( 1, nrow = n_obs, ncol = 1 )
  n_cov  <- ncol( X ) 
  
  # Initialize the model  
  varphi. <- array( 0, dim = c( n_taxa, ncol( X ) , samples ) )
  beta_gamma. <- array( 0, dim = c( n_taxa, ncol( X )  , samples ) )
  
  # Adjust initial values for  beta_gamma, varphi,  
  varphi.[ , 1, ]  <- 1
  beta_gamma.[ , , 1]  <-  rnorm( n = n_taxa*ncol( X )  )*varphi.[ , ,1]  
  
  # initialize latent variables
  cc. <- array( 0, dim = c( nrow( Z ), n_taxa, samples ) )  
  cc.[ , , 1 ] <-  Z 
  uu <- rep( 0 , nrow( Z ) ) 
  for( n in 1:nrow( Z ) ){
    sum_Z <- sum( Z[ n, ] )
    uu[ n ] <- rgamma( 1, sum_Z, sum_Z );
  } 
  
  # Memory for ZI variable eta
  eta. <- array( 0, dim = c( nrow( Z ), n_taxa, samples ) )  
  
  # Initialize eta at 1
  eta.[,,1] <- matrix( 1, nrow = nrow( Z ), ncol = ncol( Z ) )  
  
  eta_index <- list( ) 
  for( b in 1:n_taxa ){
    eta_index[[b]]  <- which( Z[ , b ]  == 0 ) - 1
    eta.[ which( Z[ ,b ]  == 0 ), b , 1] <- rbinom( 1, 1, 0.5 )
  } 
  
  eta.[,,1] <- 1 
  
  MCMC_DM <-  dm( iterations, thin, Z, X, beta_gamma., cc., uu, eta., eta_index, varphi., sigma2_beta_gamma )
 
  
  out_DM <<- MCMC_DM[1:8]
  names( out_DM ) <- c( "varphi", "beta_gamma", "eta", "cc")
  
  return( out_DM )
}
  
DMbvs_R = function(   Z = NULL, 
                      X = NULL,  
                      iterations = 10000, 
                      thin = 10, 
                      sigma2_beta_gamma = sqrt( 5 ),  
                      a_varphi = 1,
                      b_varphi = 1, 
                      seed = 1 ){
  
  set.seed( seed )
  
  # Simulate data 
  if( is.null( Z ) | is.null( X ) ){
    stop("Bad input: Data must be supplied")
  }
  
  if( sigma2_beta_gamma < 0 ){
    stop("Bad input: Variances must be positive")
  }
  
  # Get dimensions 
  n_obs <- nrow( Z )
  n_taxa <- ncol( Z )
  samples <- floor( iterations/ thin )    
  
  # Set intercept term
  X <- cbind( 1, X ) 
  n_cov  <- ncol( X )  
  
  # Initialize the model  
  varphi. <- array( 0, dim = c( n_taxa, ncol( X ) , samples ) )
  beta_gamma. <- array( 0, dim = c( n_taxa, ncol( X )  , samples ) )
  
  # Adjust initial values for  beta_gamma, varphi,  
  varphi.[ , 1, ]  <- 1
  beta_gamma.[ , , 1]  <-  rnorm( n = n_taxa*ncol( X )  )*varphi.[ , ,1]  
  
  # initialize latent variables
  cc. <- array( 0, dim = c( nrow( Z ), n_taxa, samples ) )  
  cc.[ , , 1 ] <-  Z 
  uu <- rep( 0 , nrow( Z ) ) 
  for( n in 1:nrow( Z ) ){
    sum_Z <- sum( Z[ n, ] )
    uu[ n ] <- rgamma( 1, sum_Z, sum_Z );
  } 
  
  # Memory for ZI variable eta
  eta. <- array( 0, dim = c( nrow( Z ), n_taxa, samples ) )  
  
  # Initialize eta at 1
  eta.[,,1] <- matrix( 1, nrow = nrow( Z ), ncol = ncol( Z ) )  
  
  # Get index list for eta 
  eta_accept. <- array( 0, dim = c( nrow( Z ), n_taxa, samples ) ) 
  
  eta_index <- list( ) 
  for( b in 1:n_taxa ){
    eta_index[[b]]  <- which( Z[ , b ]  == 0 ) - 1
    eta.[ which( Z[ ,b ]  == 0 ), b , 1] <- rbinom( 1, 1, 0.5 )
  }
  
 
  eta.[,,1] <- 1 # 
  
  MCMC <- dm_bvs( iterations, thin, Z, X, beta_gamma., cc., uu, eta., eta_index, varphi., sigma2_beta_gamma, a_varphi, b_varphi )
  out <<- MCMC[1:4]
  names( out ) <- c( "varphi", "beta_gamma", "eta", "cc")
   
  return(out)
}

  
estimates_ZIDM = function( zidm_obj = NULL,  
                           burnin = 1, 
                           CI = 0.95
 ){
  # Set dimensions
  n_obs = dim(zidm_obj$cc)[1]
  n_taxa = dim(zidm_obj$cc)[2]
  samples = dim(zidm_obj$cc)[3]
  
  # Set CI range
  upper <- (1 - CI)/2 + CI
  lower <- (1 - CI)/2
  
  # Get estimates 
  post_beta_theta <- apply(1/(1+exp(zidm_obj$beta_theta[ , , (burnin + 1):(samples) ])), c(1), mean) 
  post_beta_theta_upper <- apply(1/(1+exp(zidm_obj$beta_theta[ , , (burnin + 1):(samples) ])), c(1), function(x){quantile(x, upper, na.rm = T)}) 
  post_beta_theta_lower <- apply(1/(1+exp(zidm_obj$beta_theta[ , , (burnin + 1):(samples) ])), c(1), function(x){quantile(x, lower, na.rm = T)}) 

  exp_beta_gamma <- exp(zidm_obj$beta_gamma[ , , (burnin + 1):(samples) ]) 
  prop <- exp_beta_gamma/colSums( exp_beta_gamma )
  post_prop <- rowMeans( prop )
  post_prop_upper <- apply( prop, c(1), function(x){quantile(x, upper, na.rm = T)}) 
  post_prop_lower <- apply( prop, c(1), function(x){quantile(x, lower, na.rm = T)}) 
   
  post_Tn_rep <- array(  0, dim = c( n_obs, n_taxa, length((burnin + 1):(samples) ) )  ) 
  post_Tn <- apply( zidm_obj$cc[ , , (burnin + 1):(samples) ], c(1,3), sum ) 
  for( i in 1:n_taxa ){ post_Tn_rep[,i,] <- post_Tn }
  post_cc_norm <- zidm_obj$cc[ , , (burnin + 1):(samples) ]/post_Tn_rep
  
  post_cc <- apply( post_cc_norm , c(1,2), mean ) 
  post_cc_upper <- apply( post_cc_norm , c(1,2), function(x){quantile(x, upper, na.rm = T)})  
  post_cc_lower <- apply( post_cc_norm , c(1,2), function(x){quantile(x, lower, na.rm = T)}) 
  
  return( list(post_theta = post_beta_theta, 
          post_theta_lower = post_beta_theta_lower,
          post_theta_upper = post_beta_theta_upper,
          post_gamma = post_prop,
          post_gamma_lower = post_prop_lower,
          post_gamma_upper =  post_prop_upper,
          post_psi = post_cc,
          post_psi_lower = post_cc_lower,
          post_psi_upper = post_cc_upper ) )
}

estimates_DM = function( dm_obj = NULL,  
                           burnin = 1, 
                           CI = 0.95
){
  # Set dimensions
  n_obs = dim(dm_obj$cc)[1]
  n_taxa = dim(dm_obj$cc)[2]
  samples = dim(dm_obj$cc)[3]
  
  # Set CI range
  upper <- (1 - CI)/2 + CI
  lower <- (1 - CI)/2
  
  # Get estimates 
  exp_beta_gamma <- exp(dm_obj$beta_gamma[ , , (burnin + 1):(samples) ]) 
  prop <- exp_beta_gamma/colSums( exp_beta_gamma )
  post_prop <- rowMeans( prop )
  post_prop_upper <- apply( prop, c(1), function(x){quantile(x, upper, na.rm = T)}) 
  post_prop_lower <- apply( prop, c(1), function(x){quantile(x, lower, na.rm = T)}) 
  
  post_Tn_rep <- array(  0, dim = c( n_obs, n_taxa, length((burnin + 1):(samples) ) )  ) 
  post_Tn <- apply( dm_obj$cc[ , , (burnin + 1):(samples) ], c(1,3), sum ) 
  for( i in 1:n_taxa ){ post_Tn_rep[,i,] <- post_Tn }
  post_cc_norm <- dm_obj$cc[ , , (burnin + 1):(samples) ]/post_Tn_rep
  
  post_cc <- apply( post_cc_norm , c(1,2), mean ) 
  post_cc_upper <- apply( post_cc_norm , c(1,2), function(x){quantile(x, upper, na.rm = T)})  
  post_cc_lower <- apply( post_cc_norm , c(1,2), function(x){quantile(x, lower, na.rm = T)}) 
  
  return( list( post_gamma = post_prop,
               post_gamma_lower = post_prop_lower,
               post_gamma_upper =  post_prop_upper,
               post_psi = post_cc,
               post_psi_lower = post_cc_lower,
               post_psi_upper = post_cc_upper ) )
}
  
select_perf <- function( selected, truth ){
  
  if( any( ! selected %in% c( 0, 1 ) ) ) {
    stop("Bad input: selected should be zero or one")
  }
  if( any( ! truth %in% c( 0, 1 ) ) ) {
    stop("Bad input: truth should be zero or one")
  }
  select <- which( selected == 1 )
  not_selected <- which( selected == 0 )
  included <- which( truth == 1 )
  excluded <- which( truth == 0 )
  
  TP <- sum( select %in% included )
  TN <- sum( not_selected %in% excluded )
  FP <- sum( select %in% excluded )
  FN <- sum( not_selected %in% included )
  sensitivity <- TP/( FN + TP )
  specificity <- TN/( FP + TN )
  mcc <- ( TP*TN - FP*FN )/(sqrt( TP + FP )*sqrt(TP + FN )*sqrt(TN + FP )*sqrt(TN + FN) )
  f1 <- 2*TP/( 2*TP + FN + FP )
  
  return( list( sens = sensitivity, spec = specificity, mcc = mcc, f1 = f1 ) )
}
  

  
  
  
  
  