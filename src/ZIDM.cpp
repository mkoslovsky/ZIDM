

/* Code for Zero-Inflated Dirichlet Distribution: Models and Computational Strategies for High-Dimensional Compositional Data
* Author: Fill in when accepted
*
* Performs Bayesian variable selection for zero-inflated multivariate count data 
*
*   References and Thanks to: 
*
*   Jesse Bennett Windle
*   Forecasting High-Dimensional, Time-Varying Variance-Covariance Matrices
*   with High-Frequency Data and Sampling Polya-Gamma Random Variates for
*   Posterior Distributions Derived from Logistic Likelihoods  
*   PhD Thesis, 2013   
*
*   Damien, P. & Walker, S. G. Sampling Truncated Normal, Beta, and Gamma Densities 
*   Journal of Computational and Graphical Statistics, 2001, 10, 206-215
*
*   Chung, Y.: Simulation of truncated gamma variables 
*   Korean Journal of Computational & Applied Mathematics, 1998, 5, 601-610
*
*   Makalic, E. & Schmidt, D. F. High-Dimensional Bayesian Regularised Regression with the BayesReg Package 
*   arXiv:1611.06649 [stat.CO], 2016 https://arxiv.org/pdf/1611.06649.pdf 
*/
  
#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma; 

// Mathematical constants
#define MATH_PI        3.141592653589793238462643383279502884197169399375105820974
#define MATH_PI_2      1.570796326794896619231321691639751442098584699687552910487
#define MATH_2_PI      0.636619772367581343075535053490057448137838582961825794990
#define MATH_PI2       9.869604401089358618834490999876151135313699407240790626413
#define MATH_PI2_2     4.934802200544679309417245499938075567656849703620395313206
#define MATH_SQRT1_2   0.707106781186547524400844362104849039284835937688474036588
#define MATH_SQRT_PI_2 1.253314137315500251207882642405522626503493370304969158314
#define MATH_LOG_PI    1.144729885849400174143427351353058711647294812915311571513
#define MATH_LOG_2_PI  -0.45158270528945486472619522989488214357179467855505631739
#define MATH_LOG_PI_2  0.451582705289454864726195229894882143571794678555056317392


// Define helper functions
namespace help{

// Simulate MVT normal data
arma::mat mvrnormArma( int n, arma::vec mu, arma::mat sigma ) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn( n, ncols );
  return arma::repmat( mu, 1, n ).t()  + Y * arma::chol( sigma );
}

// Generate exponential distribution random variates
double exprnd(double mu)
{
  return -mu * (double)std::log(1.0 - (double)R::runif(0.0,1.0));
}

// Function a_n(x) defined in equations (12) and (13) of
// Bayesian inference for logistic models using Polya-Gamma latent variables
// Nicholas G. Polson, James G. Scott, Jesse Windle
// arXiv:1205.0310
//
// Also found in the PhD thesis of Windle (2013) in equations
// (2.14) and (2.15), page 24
double aterm(int n, double x, double t)
{
  double f = 0;
  if(x <= t) {
    f = MATH_LOG_PI + (double)std::log(n + 0.5) + 1.5*(MATH_LOG_2_PI- (double)std::log(x)) - 2*(n + 0.5)*(n + 0.5)/x;
  }
  else {
    f = MATH_LOG_PI + (double)std::log(n + 0.5) - x * MATH_PI2_2 * (n + 0.5)*(n + 0.5);
  }    
  return (double)exp(f);
}

// Generate inverse gaussian random variates
double randinvg(double mu)
{
  // sampling
  double u = R::rnorm(0.0,1.0);
  double V = u*u;
  double out = mu + 0.5*mu * ( mu*V - (double)std::sqrt(4.0*mu*V + mu*mu * V*V) );
  
  if(R::runif(0.0,1.0) > mu /(mu+out)) {    
    out = mu*mu / out; 
  }    
  return out;
}

// Sample truncated gamma random variates
// Ref: Chung, Y.: Simulation of truncated gamma variables 
// Korean Journal of Computational & Applied Mathematics, 1998, 5, 601-610
double truncgamma()
{
  double c = MATH_PI_2;
  double X, gX;
  
  bool done = false;
  while(!done)
  {
    X = help::exprnd(1.0) * 2.0 + c;
    gX = MATH_SQRT_PI_2 / (double)std::sqrt(X);
    
    if(R::runif(0.0,1.0) <= gX) {
      done = true;
    }
  }
  
  return X;  
}

// Sample truncated inverse Gaussian random variates
// Algorithm 4 in the Windle (2013) PhD thesis, page 129
double tinvgauss(double z, double t)
{
  double X, u;
  double mu = 1.0/z;
  
  // Pick sampler
  if(mu > t) {
    // Sampler based on truncated gamma 
    // Algorithm 3 in the Windle (2013) PhD thesis, page 128
    while(1) {
      u = R::runif(0.0, 1.0);
      X = 1.0 / help::truncgamma();
      
      if ((double)std::log(u) < (-z*z*0.5*X)) {
        break;
      }
    }
  }  
  else {
    // Rejection sampler
    X = t + 1.0;
    while(X >= t) {
      X = help::randinvg(mu);
    }
  }    
  return X;
}


// Sample PG(1,z)
// Based on Algorithm 6 in PhD thesis of Jesse Bennett Windle, 2013
// URL: https://repositories.lib.utexas.edu/bitstream/handle/2152/21842/WINDLE-DISSERTATION-2013.pdf?sequence=1
double samplepg(double z)
{
  //  PG(b, z) = 0.25 * J*(b, z/2)
  z = (double)std::fabs((double)z) * 0.5;
  
  // Point on the intersection IL = [0, 4/ log 3] and IR = [(log 3)/pi^2, \infty)
  double t = MATH_2_PI;
  
  // Compute p, q and the ratio q / (q + p)
  // (derived from scratch; derivation is not in the original paper)
  double K = z*z/2.0 + MATH_PI2/8.0;
  double logA = (double)std::log(4.0) - MATH_LOG_PI - z;
  double logK = (double)std::log(K);
  double Kt = K * t;
  double w = (double)std::sqrt(MATH_PI_2);
  
  double logf1 = logA + R::pnorm(w*(t*z - 1),0.0,1.0,1,1) + logK + Kt;
  double logf2 = logA + 2*z + R::pnorm(-w*(t*z+1),0.0,1.0,1,1) + logK + Kt;
  double p_over_q = (double)std::exp(logf1) + (double)std::exp(logf2);
  double ratio = 1.0 / (1.0 + p_over_q); 
  
  double u, X;
  
  // Main sampling loop; page 130 of the Windle PhD thesis
  while(1) 
  {
    // Step 1: Sample X ? g(x|z)
    u = R::runif(0.0,1.0);
    if(u < ratio) {
      // truncated exponential
      X = t + help::exprnd(1.0)/K;
    }
    else {
      // truncated Inverse Gaussian
      X = help::tinvgauss(z, t);
    }
    
    // Step 2: Iteratively calculate Sn(X|z), starting at S1(X|z), until U ? Sn(X|z) for an odd n or U > Sn(X|z) for an even n
    int i = 1;
    double Sn = help::aterm(0, X, t);
    double U = R::runif(0.0,1.0) * Sn;
    int asgn = -1;
    bool even = false;
    
    while(1) 
    {
      Sn = Sn + asgn * help::aterm(i, X, t);
      
      // Accept if n is odd
      if(!even && (U <= Sn)) {
        X = X * 0.25;
        return X;
      }
      
      // Return to step 1 if n is even
      if(even && (U > Sn)) {
        break;
      }
      
      even = !even;
      asgn = -asgn;
      i++;
    }
  }
  return X;
}

// Call sample function from R
double sample_cpp( IntegerVector x ){
  Function f("sample");
  IntegerVector sampled = f(x, Named("size") = 1);
  return sampled[0];
}

double sample_prob_cpp( IntegerVector x, NumericVector prob ){
  Function f("sample");
  IntegerVector sampled = f(x, Named("size") = 1, Named("prob") = prob);
  return sampled[0];
}
 
// Get log derivative of matrix
double log_det_cpp( arma::mat x ){
  int p = x.n_cols;
  double sum = 0;
  
  arma::mat chole = chol(x);
  
  for( int j = 0; j < p; ++j ){
    sum += 2*log( chole( j, j ) );
  }
  
  return sum;
}
 
// Make a sequence
IntegerVector myseq( int first, int last) {
  IntegerVector y(abs(last - first) + 1);
  if (first < last) 
    std::iota(y.begin(), y.end(), first);
  else {
    std::iota(y.begin(), y.end(), last);
    std::reverse(y.begin(), y.end());
  }
  return y;
}
 

// Function :: Calculate log likelihood of X ( compositional ) across i
double log_multinomial_cpp( arma::mat z, arma::mat psi ){
  int n = psi.n_rows;
  arma::vec log_mult( n );
  
  for( int i = 0; i < n; ++i ){
    log_mult( i ) = lgamma( sum( z.row( i ) ) + 1 ) + sum( - lgamma( z.row( i ) + 1 ) + z.row( i )%log( psi.row( i )) );
  }
  
  double sum_log_mult = sum( log_mult );
  
  return sum_log_mult;
}
 

// Function :: Propose a new beta
double beta_prop_cpp( double beta ){
  // Simulate proposal with jitter around current beta_gamma
  // Make sure that beta_gamma is specified with iteration iterations and branch b
  double beta_prop = beta + sqrt(0.5)*rnorm( 1 )[ 0 ];
  
  // Return output
  return beta_prop ;
}

// Function :: Calculate log beta_gamma contribution
double log_beta_gamma_cpp( NumericMatrix beta_gamma, 
                           double sigma2_beta_gamma, 
                           NumericMatrix varphi ){
  
  // initiate log_beta_gamma and set size of vector
  double log_beta_gamma = 0;
  int B = beta_gamma.rows();
  int P = beta_gamma.cols();
  
  // Sum over all beta_gamma in B and P
  for( int b = 0; b < B; ++b ) {
    for( int p = 0; p < P; ++p ) {
      if( varphi( b, p ) == 1){
        log_beta_gamma += -0.50*log( 2*atan(1)*4*sigma2_beta_gamma ) - 1/( 2*sigma2_beta_gamma )*pow( beta_gamma( b, p ), 2 );
      }
    }
  }
  
  // Return output
  return log_beta_gamma;
}

// Function :: Calculate log beta contribution for individual p
double log_beta_p_cpp( double beta, 
                       double sigma2_beta){
  
  double log_beta = -0.50*log( 2*atan(1)*4*sigma2_beta ) - 1/( 2*sigma2_beta )*pow( beta, 2 );
  
  // Return output
  return log_beta;
}

// Function :: Calculate individual varphi
double log_varphi_pj_cpp( double t_pj, 
                          double a_varphi, 
                          double b_varphi
){
  double post_a_varphi = t_pj + a_varphi;
  double post_b_varphi = 1 - t_pj + b_varphi;
  double log_varphi_pj = lgamma( post_a_varphi ) + lgamma( post_b_varphi ) - lgamma( post_a_varphi + post_b_varphi ) - ( lgamma( a_varphi ) + lgamma( b_varphi ) - lgamma( a_varphi + b_varphi ) );
  
  // Return output
  return log_varphi_pj ;
}

// Function :: Calculate individual zeta
double log_zeta_p_cpp( double t_p, 
                          double a_zeta, 
                          double b_zeta
){
  double post_a_zeta = t_p + a_zeta;
  double post_b_zeta = 1 - t_p + b_zeta;
  double log_zeta_p = lgamma( post_a_zeta ) + lgamma( post_b_zeta ) - lgamma( post_a_zeta + post_b_zeta ) - ( lgamma( a_zeta ) + lgamma( b_zeta ) - lgamma( a_zeta + b_zeta ) );
  
  // Return output
  return log_zeta_p ;
}

 
// Function :: Between Step (jointly update beta_gamma and varphi) with add/delete using BB prior
List between_beta_gamma_varphi_augment_update_BB_cpp(
    arma::mat x,
    arma::mat varphi, 
    arma::mat beta_gamma, 
    arma::mat loggamma, 
    arma::mat cc,
    arma::mat eta, 
    double sigma2_beta_gamma,
    double a_varphi,
    double b_varphi
){
  
  // Between Model Step
  // Initiate memory space
  int J = varphi.n_rows;
  int P = varphi.n_cols; 
  int N = x.n_rows;
  
  // Get a covariate p to update in branch j 
  
  int j = help::sample_cpp( help::myseq( 0, J - 1) ); 
  for( int p = 1; p < P; ++p ){   // p = 1 skips intercept term
    
    // Calculate current log likelihood 
    double log_like = 0; 
    
    for( int n = 0; n < N; ++n ){
      if( eta( n, j ) == 1 ){
        log_like = log_like - lgamma( exp( loggamma( n , j ) ) ) + exp( loggamma( n , j ) )*log( cc( n, j ) );
      }
    }
    
    // Calculate proposed likelihood 
    double beta_gamma_prop;
    int varphi_prop; 
    
    if( varphi( j, p ) == 1 ){
      beta_gamma_prop = 0.0;
      varphi_prop = 0;
    }else{
      beta_gamma_prop = help::beta_prop_cpp( beta_gamma( j, p ) ); 
      varphi_prop = 1;
    }
    
    arma::vec loggamma_prop( N );
    loggamma_prop.zeros();
    
    for( int n = 0 ; n < N ; ++n ){
      loggamma_prop[ n ] = loggamma( n, j ) - beta_gamma( j, p ) * x( n, p ) + beta_gamma_prop * x( n, p );
    }
    
    double log_like_prop = 0; 
    for( int n = 0; n < N; ++n ){
      if( eta( n, j ) == 1 ){
        log_like_prop = log_like_prop - lgamma( exp( loggamma_prop[ n ] ) ) + exp( loggamma_prop[ n ] )*log( cc( n, j ) );
      }
    }
    
    // Calculate ratio for single varphi
    double r;
    
    // If BB
    if( varphi( j, p ) == 1 ){
      r = log_like_prop + help::log_varphi_pj_cpp( 0, a_varphi, b_varphi) - ( log_like + help::log_beta_p_cpp( beta_gamma( j, p ), sigma2_beta_gamma ) + help::log_varphi_pj_cpp( 1, a_varphi, b_varphi ) );
    }else{
      r = log_like_prop + help::log_beta_p_cpp( beta_gamma_prop, sigma2_beta_gamma ) + help::log_varphi_pj_cpp( 1, a_varphi, b_varphi) - ( log_like + help::log_varphi_pj_cpp( 0, a_varphi, b_varphi ) );
    }
    
    
    // Calculate acceptance probability
    double aa  = log( runif( 1 )[ 0 ] );
    
    // Determine acceptance
    if( aa < r ){
      varphi( j, p ) = varphi_prop;
      beta_gamma( j, p ) = beta_gamma_prop;
      loggamma.col( j ) = loggamma_prop;
    }
    
  }
  
  // Return output
  List between( 3 );
  between[ 0 ] = varphi;
  between[ 1 ] = beta_gamma;
  between[ 2 ] = loggamma;
  return between;
}
 

// Function :: Within Step ( update beta_gamma )
List within_beta_gamma_augment_cpp(
    arma::mat x,
    arma::mat varphi, 
    arma::mat beta_gamma, 
    arma::mat loggamma,
    arma::mat cc,
    arma::mat eta, 
    double sigma2_beta_gamma
){
  
  // Initiate memory space
  int J = varphi.n_rows;
  int P = varphi.n_cols;
  int N = x.n_rows;  
  
  // Get which covariates are included
  for( int j = 0; j < J; ++j ){
    for( int p = 0; p < P; ++p ){ // Update intercept too 
      
      if( varphi( j , p ) == 1 ){
        // Calculate current log likelihood 
        double log_like = 0; 
        
        for( int n = 0; n < N; ++n ){
          if( eta( n, j ) == 1 ){
            log_like = log_like - lgamma( exp( loggamma( n , j ) ) ) + exp( loggamma( n , j ) )*log( cc( n, j ) );
          }  
        }
        
        // Calculate proposed likelihood 
        double beta_gamma_prop = 0; 
        
        beta_gamma_prop = help::beta_prop_cpp( beta_gamma( j, p ) ); 
        
        arma::vec loggamma_prop( N );
        loggamma_prop.zeros();
        
        for( int n = 0 ; n < N ; ++n ){
          loggamma_prop[ n ] = loggamma( n, j ) - beta_gamma( j, p ) * x( n, p ) + beta_gamma_prop * x( n, p );
        }
        
        double log_like_prop = 0; 
        for( int n = 0; n < N; ++n ){
          if( eta( n, j ) == 1 ){
            log_like_prop = log_like_prop - lgamma( exp( loggamma_prop[ n ] ) ) + exp( loggamma_prop[ n ] )*log( cc( n, j ) );
          }
        }
        
        // Calculate ratio
        double r = 0;  
        r = log_like_prop + help::log_beta_p_cpp( beta_gamma_prop, sigma2_beta_gamma) - ( log_like + help::log_beta_p_cpp( beta_gamma( j, p ), sigma2_beta_gamma )  );
 
        // Calculate acceptance probability
        double a  = log( runif( 1 )[ 0 ] );
        
 
        // Determine acceptance
        if(a < r){ 
          beta_gamma( j, p ) = beta_gamma_prop;
          loggamma.col( j ) = loggamma_prop;
        }
      }
    }
  }

  // Return output
  List beta_gamma_out( 2 );
  beta_gamma_out[ 0 ] = beta_gamma;
  beta_gamma_out[ 1 ] = loggamma;
  return beta_gamma_out;
}

// Update latent variable cc
arma::mat cc_update_cpp( arma::mat z, arma::mat loggamma, arma::vec uu, arma::mat eta  ){
  int N = z.n_rows;
  int J = z.n_cols;
  
  arma::mat cc( N, J );
  cc.zeros();
  
  for( int n = 0; n < N; ++n ){
    for( int j = 0; j < J; ++j ){
      if( eta( n, j ) == 1 ){ 
        cc( n, j ) = rgamma( 1, z( n, j ) + exp(loggamma( n, j ) ), 1/( uu[ n ] + 1 )  )[ 0 ]; 
        if( cc( n, j ) < pow(10,-100)){
          cc( n, j ) = pow(10,-100);
        }
      }
    }
  } 
  
  return cc;
}

// Update latent variable uu 
arma::vec uu_update_cpp( arma::mat z, arma::mat cc ){
  int N = z.n_rows; 
  
  arma::vec uu( N );
  uu.zeros();
  
  for( int n = 0; n < N; ++n ){
    double z_dot_n = sum( z.row( n ) );
    double T_n = sum( cc.row( n ) );
    
    uu[ n ] = rgamma( 1, z_dot_n, 1/T_n )[ 0 ]; 
  }
  
  return uu; 
}

// Joint update for cc and eta
List between_cc_eta_update_cpp( arma::mat z, 
                                arma::mat x_theta, 
                                arma::mat beta_theta, 
                                arma::mat cc, 
                                arma::mat eta,  
                                List eta_index, 
                                arma::mat loggamma,
                                arma::mat varphi  ){
  // Get dimensions
    int B = eta.n_cols;   
    int N = eta.n_rows;
    arma::mat eta_accept( N, B );
    eta_accept.zeros();
  
  // Add/Delete for eta and c 
  for( int b = 0; b < B; ++b ){   // For updating all of the indicators
    // int b = help::sample_cpp( help::myseq( 0, B - 1) ); // For updating one at a time
    
    // Get current eta value and z_dot and tt
      arma::vec eta_index_b = as<arma::vec>( eta_index[ b ] );
      int num_zero_b = eta_index_b.size();
      for( int n = 0; n < num_zero_b; ++n ){
        int i_index = eta_index_b[ n ]; 
        int eta_c = eta( i_index, b );
        double cc_c = cc( i_index, b );
        double z_dot_n = sum( z.row( i_index ) );
        double T_n = sum( cc.row( i_index ) );
        arma::mat theta_i = exp( x_theta.row( i_index )*beta_theta.row( b ).t() )/( 1 + exp( x_theta.row( i_index )*beta_theta.row( b ).t() ) );
        double a_c_prop_scaled = exp( loggamma( i_index, b ) );  
        
        // Delete Step:
        if ( eta_c == 1 ){
          
          // Propose
          double cc_prop = 0;
          double eta_prop = 0;  
          double r = -z_dot_n*log( 1 - cc_c/T_n ) + lgamma( exp( loggamma( i_index, b ) ) ) - lgamma( a_c_prop_scaled ) + ( a_c_prop_scaled - exp( loggamma( i_index, b ) ) )*log( cc_c ) + log( 1 - theta_i[ 0 ] ) - log( theta_i[ 0 ] );
          
          // Calculate acceptance probability
          double a  = log( runif( 1 )[ 0 ] );
          
          // Determine acceptance
          if(a < r){
            cc( i_index, b ) = cc_prop;
            eta( i_index, b ) = eta_prop;
            eta_accept( i_index, b ) = 1; 
          }
          
        }else{ // Add step 
          
          // Propose
          double cc_prop = rgamma( 1, a_c_prop_scaled, 1 )[ 0 ];
          double eta_prop = 1;
          
          double r = -z_dot_n*log( 1 + cc_prop/T_n ) + lgamma( a_c_prop_scaled ) - lgamma( exp( loggamma( i_index, b ) ) ) + ( exp( loggamma( i_index, b ) ) - a_c_prop_scaled )*log( cc_prop ) - log( 1 - theta_i[ 0 ] ) + log( theta_i[ 0 ] );
          
          // Calculate acceptance probability
          double a  = log( runif( 1 )[ 0 ] );
          
          // Determine acceptance
          if(a < r){
            cc( i_index, b ) = cc_prop;
            eta( i_index, b ) = eta_prop;
            eta_accept( i_index, b ) = 1;
          }
        }
        
        }
  } 
  
  // Return output
  List cc_eta_out( 3 );
  cc_eta_out[ 0 ] = cc;
  cc_eta_out[ 1 ] = eta; 
  cc_eta_out[ 2 ] = eta_accept; 
  return cc_eta_out; 
}

// Update auxiliary parameters omega 
arma::mat update_omega( arma::mat x_theta,  
                        arma::mat beta_theta 
    ){ 
  int B = beta_theta.n_rows;
  int obs = x_theta.n_rows;
  
  // Make a home for W updates
  arma::mat updated_omega( obs, B ); 
  updated_omega.zeros();
   
  
  // Update each W individually
  for( int b = 0; b < B; ++b ){ 
    
    for( int n = 0; n < obs; ++n ){ 
      arma::mat x_sub = x_theta.row( n );
      arma::mat phi_i( 1,1 );
      phi_i = x_sub*beta_theta.row( b ).t(); 
      updated_omega( n, b ) = help::samplepg( phi_i[ 0 ] );
    }
  }
  return updated_omega;
}

 

// Function :: Update beta_theta Within fixed
arma::mat within_beta_theta( arma::mat x_theta,
                             arma::mat zeta,  
                             arma::mat omega, 
                             arma::mat eta,
                             double sigma2_beta_theta ){
  int P_theta = x_theta.n_cols;
  int B = zeta.n_rows;
  int obs = omega.n_rows;

  // Make place for update
  arma::mat beta_theta_update( B, P_theta );
  beta_theta_update.zeros();
  

  for( int b = 0; b < B; ++b ){ 
   
    // Make h_ij = k_ij/w_ij, k_ij = y_ij - 1/2
    arma::vec H = ( eta.col( b ) - 0.5)/omega.col( b );

    // Get number of active terms in zeta b
    int zeta_b_dim = sum( zeta.row( b ) ) ;

    // Adjust dimension of x and sigma2 by zeta and set to zero
    arma::mat x_theta_b_zeta( obs, zeta_b_dim );
    arma::mat sigma2_beta_theta_MAT( zeta_b_dim, zeta_b_dim );
    x_theta_b_zeta.zeros();
    sigma2_beta_theta_MAT.zeros();

    int count_b = 0;
    for( int p = 0; p < P_theta; ++p){
      if( zeta( b, p ) == 1 ){
        x_theta_b_zeta.col( count_b ) = x_theta.col( p );
        sigma2_beta_theta_MAT( count_b, count_b) = 1/sigma2_beta_theta;
        count_b += 1;
      }
    }

    // Make V_beta_theta and mu_beta_theta (and inside)
    arma::mat V_beta_theta( zeta_b_dim, zeta_b_dim );
    arma::mat mu_beta_theta( zeta_b_dim, 1 );
    arma::mat mu_beta_theta_inside( zeta_b_dim, 1 );
    V_beta_theta.zeros();
    mu_beta_theta.zeros();
    mu_beta_theta_inside.zeros();

    // Update for each individuals beta_temp
    for( int n = 0; n < obs; ++n ){
      V_beta_theta += omega( n , b )*x_theta_b_zeta.row( n ).t()*x_theta_b_zeta.row( n );
      mu_beta_theta_inside += omega( n , b )*x_theta_b_zeta.row( n ).t()*H[ n ];
    }

    V_beta_theta += sigma2_beta_theta_MAT;
    V_beta_theta = inv( V_beta_theta );
    mu_beta_theta = mu_beta_theta_inside;

    mu_beta_theta = V_beta_theta*mu_beta_theta;

    // Update bth element
    arma::mat beta_zeta_b_update( 1, zeta_b_dim );
    beta_zeta_b_update.zeros();

    beta_zeta_b_update =  help::mvrnormArma( 1, mu_beta_theta, V_beta_theta );

    // Re-index the beta updates appropriately
    count_b = 0;
    for( int p = 0; p < P_theta; ++p){
      if( zeta( b, p ) == 1 ){
        beta_theta_update( b, p ) = beta_zeta_b_update[ count_b ];
        count_b += 1;
      }
    }
  }

  return beta_theta_update;
}

// Function :: Between Step (jointly update beta_theta and zeta) with add/delete using BB prior
List between_beta_theta_zeta_update_BB_cpp(
    arma::mat x_theta,
    arma::mat zeta, 
    arma::mat beta_theta, 
    arma::mat eta,  
    double sigma2_beta_theta,
    double a_zeta,
    double b_zeta
){
  
  // Between Model Step
  // Initiate memory space
  int P_theta = zeta.n_cols;  
  int B = zeta.n_rows;
  int obs = x_theta.n_rows;
  
  // Choose a compositional element
   int b = help::sample_cpp( help::myseq( 0, B - 1) ); 
  //for( int b = 1; b < B; ++b ){

  
  // For each covariate, update
   for( int p = 1; p < P_theta; ++p ){   // p = 1 skips intercept term 
             //int p = help::sample_cpp( help::myseq( 1, P_theta - 1) );
    
    // Calculate current log likelihood of eta
    double log_like = 0;  
    
    for( int n = 0; n < obs; ++n ){ 
        arma::mat theta_i = exp( x_theta.row( n )*beta_theta.row( b ).t() )/( 1 + exp( x_theta.row( n )*beta_theta.row( b ).t() ) );
        log_like = log_like  + eta( n, b )*log( theta_i[ 0 ] ) + ( 1 - eta( n, b ) )*log( 1 - theta_i[ 0 ] ) ;
    }
    
    // Calculate proposed likelihood 
    arma::mat beta_theta_prop = beta_theta.row( b );
    arma::mat zeta_prop = zeta.row( b );
    
    if( zeta( b, p ) == 1 ){
      beta_theta_prop[ p ] = 0.0;
      zeta_prop[ p ] = 0;
    }else{
      beta_theta_prop[ p ] = help::beta_prop_cpp( beta_theta( b, p ) ); 
      zeta_prop[ p ] = 1;
    } 

    double log_like_prop = 0; 
    
    for( int n = 0; n < obs; ++n ){ 
      arma::mat theta_i_prop = exp( x_theta.row( n )*beta_theta_prop.t() )/( 1 + exp( x_theta.row( n )*beta_theta_prop.t() ) );
      log_like_prop = log_like_prop  + eta( n, b )*log( theta_i_prop[ 0 ] ) + ( 1 - eta( n, b ) )*log( 1 - theta_i_prop[ 0 ] ) ;
    }
 
    // Calculate ratio for single varphi
    double r = 0; 
    
    // If BB
    if( zeta( b, p ) == 1 ){
      r = log_like_prop + help::log_zeta_p_cpp( 0, a_zeta, b_zeta) - ( log_like + help::log_beta_p_cpp( beta_theta( b, p ), sigma2_beta_theta ) + help::log_zeta_p_cpp( 1, a_zeta, b_zeta ) );
    }else{
      r = log_like_prop + help::log_beta_p_cpp( beta_theta_prop[ p ], sigma2_beta_theta ) + help::log_zeta_p_cpp( 1, a_zeta, b_zeta ) - ( log_like + help::log_zeta_p_cpp( 0, a_zeta, b_zeta ) );
    }
    
    // Calculate acceptance probability
    double aa  = log( runif( 1 )[ 0 ] ); 
     
    // Determine acceptance
    if( aa < r ){
      zeta( b, p ) = zeta_prop[ p ];
      beta_theta( b, p ) = beta_theta_prop[ p ];
    }
  }
  //}
  
  // Return output
  List between( 2 );
  between[ 0 ] = zeta;
  between[ 1 ] = beta_theta; 
  return between;
}


IntegerVector sample_mult_cpp( IntegerVector x, int num ){
  // Calling sample()
  Function f( "sample" );
  IntegerVector sampled = f( x, Named( "size" ) = num, Named( "replace" ) = false );
  return sampled;
}

IntegerVector oneMultinomC(NumericVector probs) {
  int k = probs.size();
  SEXP ans;
  PROTECT(ans = Rf_allocVector(INTSXP, k));
  probs = Rf_coerceVector(probs, REALSXP);
  rmultinom(1, REAL(probs), k, &INTEGER(ans)[0]);
  UNPROTECT(1);
  return(ans);
}

void set_seed(double seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(std::floor(std::fabs(seed)));
}

double x_w( int B, int w, int zeta, int n ){
  double x_w_val = Rf_choose( B - w, B - zeta)/Rf_choose( n + B - w - 1, n );
  return x_w_val;
}

double w_x( int B, int w, int zeta, int n ){
   arma::vec w_x_values( zeta + 1);
   w_x_values.zeros();
   double w_x_val = 0;
   w_x_val = x_w( B, w, zeta, n ); 
   
   for( int i = 0; i < zeta + 1; ++i){
     w_x_values[ i ] = x_w( B, i, zeta, n );
    }
   
   w_x_val = w_x_val/sum( w_x_values );
   return w_x_val;
} 
 
double Enum( double x_i, double B, double zeta, double n, double w ){
  double val = ( ( x_i + 1 )/( n + B - w ) )*w_x( B, w, zeta, n ); 
  
  return val;
}

double Ezero( double B, double zeta, double n, double w ){
  double val = ( 1 - w/zeta )/( n + B - w )*w_x( B, w, zeta, n ); 
  return val;
  }

arma::mat rescaled( arma::mat beta_gamma){
  double sum_beta_gamma_0 = 0;
  int B = beta_gamma.n_rows; 
  int P = beta_gamma.n_cols; 
  arma::mat beta_gamma_hold( B, P );
  beta_gamma_hold = beta_gamma;
  
  for( int b = 0; b < B; ++b ){
    sum_beta_gamma_0 += abs( beta_gamma( b, 0 ) );
  }
  
  for( int b = 0; b < B; ++b ){
    beta_gamma_hold( b, 0 ) = beta_gamma( b, 0 )/sum_beta_gamma_0;
  }
 
  
  return beta_gamma_hold;
}

arma::mat rescaled2( arma::mat cc){

  int B = cc.n_cols; 
  int N = cc.n_rows;  
  arma::mat cc_hold( N, B );
  cc_hold.zeros(); 
  
  for( int i = 0; i < N; ++i ){
    double T_n = sum( cc.row( i ) );
     for( int b = 0; b < B; ++b ){
      cc_hold( i, b ) = cc( i, b )/T_n;
    }
  }
  
  return cc_hold;
}
 
 
} // For namespace 'help'

 
// Function :: MCMC algorithm
// Make sure arguments are indexed by iteration
// [[Rcpp::export]]
List zidm_bvs(
    int iterations,
    int thin, 
    arma::mat z,
    arma::mat x,
    arma::mat x_theta,
    arma::cube beta_gamma,  
    arma::cube cc,
    arma::vec temp_uu, 
    arma::cube eta, 
    arma::cube eta_accept,
    List eta_index,
    arma::cube omega,
    arma::cube beta_theta,
    arma::cube zeta, 
    arma::cube varphi, 
    double sigma2_beta_gamma, 
    double a_varphi, 
    double b_varphi,     
    double a_zeta,  
    double b_zeta,   
    double sigma2_beta_theta
){
  
  // Initiate memory
  List between_beta_gamma_varphi( 3 ); 
  List beta_gamma_out( 2 );
  List between_cc_eta( 3 );
  List between_beta_theta_zeta( 2 );
  int B = z.n_cols;
  int P = varphi.n_cols;
  int P_theta = x_theta.n_cols;
  int N = x.n_rows;   
  arma::mat temp_varphi( B, P );
  arma::mat temp_beta_gamma( B, P ); 
  arma::mat temp_eta( N, B ); 
  arma::mat temp_omega( N, B );
  arma::mat temp_beta_theta( B, P_theta );
  arma::mat temp_zeta( B, P_theta );
  arma::mat temp_eta_accept( B, P );
  arma::mat temp_cc( N, B );
   
  temp_varphi = varphi.slice( 0 );
  temp_beta_gamma = beta_gamma.slice( 0 );  
  temp_eta = eta.slice( 0 );
  temp_omega = omega.slice( 0 );
  temp_beta_theta = beta_theta.slice( 0 );
  temp_zeta = zeta.slice( 0 );
  temp_cc = cc.slice( 0 );  
  temp_eta_accept = eta_accept.slice( 0 );
 
  arma::mat temp_loggamma( N, B );
  temp_loggamma.zeros();
  for( int n = 0; n < N; ++n ){
    for( int j = 0; j < B; ++j ){
      temp_loggamma( n, j ) = ( x.row( n )*temp_beta_gamma.row( j ).t() ).eval()( 0, 0 ); 
    }
  }

  // Looping over the number of iterations specified by user
  for( int iter = 0; iter < iterations; ++iter ){
    
    // Add/Delete  
       between_beta_gamma_varphi = help::between_beta_gamma_varphi_augment_update_BB_cpp( x, temp_varphi, temp_beta_gamma, temp_loggamma, temp_cc, temp_eta, sigma2_beta_gamma, a_varphi, b_varphi );
       temp_varphi = as<arma::mat>( between_beta_gamma_varphi[ 0 ] );
       temp_beta_gamma = as<arma::mat>( between_beta_gamma_varphi[ 1 ] );
       temp_loggamma = as<arma::mat>( between_beta_gamma_varphi[ 2 ] ); 
    // Update within
       beta_gamma_out = help::within_beta_gamma_augment_cpp( x, temp_varphi, temp_beta_gamma, temp_loggamma, temp_cc, temp_eta, sigma2_beta_gamma );
       temp_beta_gamma =  as<arma::mat>( beta_gamma_out[ 0 ] );
       temp_loggamma = as<arma::mat>( beta_gamma_out[ 1 ] );
    
        // Update cc within
       temp_cc = help::cc_update_cpp( z, temp_loggamma, temp_uu, temp_eta );
       

    // Update uu
       temp_uu = help::uu_update_cpp( z, temp_cc ); 
 
    // Update cc and eta between
       between_cc_eta = help::between_cc_eta_update_cpp( z, x_theta, temp_beta_theta, temp_cc, temp_eta, eta_index, temp_loggamma, temp_varphi);
       temp_cc = as<arma::mat>( between_cc_eta[ 0 ] );
       temp_eta = as<arma::mat>( between_cc_eta[ 1 ] );
       temp_eta_accept = as<arma::mat>( between_cc_eta[ 2 ] );
       
    // Update omega
       temp_omega = help::update_omega( x_theta, temp_beta_theta );

    // Within beta_theta
       temp_beta_theta = help::within_beta_theta( x_theta, temp_zeta, temp_omega, temp_eta, sigma2_beta_theta );

    // Between beta_theta and zeta
       between_beta_theta_zeta = help::between_beta_theta_zeta_update_BB_cpp( x_theta, temp_zeta, temp_beta_theta, temp_eta, sigma2_beta_theta, a_zeta, b_zeta );
       temp_zeta = as<arma::mat>( between_beta_theta_zeta[ 0 ] );
       temp_beta_theta = as<arma::mat>( between_beta_theta_zeta[ 1 ] );
       
       // Rescale beta gamma
       //temp_beta_gamma = help::rescaled( temp_beta_gamma );
       
       // Rescalecc
       //temp_cc = help::rescaled2( temp_cc );
       

    // Set the starting values for the next iteration
    if( (iter + 1) % thin == 0 ){ 
      varphi.slice( (iter + 1)/thin - 1 ) = temp_varphi;
      beta_gamma.slice( (iter + 1)/thin - 1 ) = temp_beta_gamma; 
      eta.slice( (iter + 1)/thin - 1 ) = temp_eta; 
      omega.slice( (iter + 1)/thin - 1 ) = temp_omega;
      beta_theta.slice( (iter + 1)/thin - 1 ) = temp_beta_theta;
      zeta.slice( (iter + 1)/thin - 1 ) = temp_zeta; 
      cc.slice( (iter + 1)/thin - 1 ) = temp_cc; 
      eta_accept.slice( (iter + 1)/thin - 1 ) = temp_eta_accept; ; // Incorporates thinning
    }
     
    
    // Print out progress
    double printer = iter % 250;
    
    if( printer == 0 ){
      Rcpp::Rcout << "Iteration = " << iter << std::endl;
    }
  }
  
  // Return output
  List output( 8 ); 
  output[ 0 ] = varphi;
  output[ 1 ] = beta_gamma; 
  output[ 2 ] = eta;
  output[ 3 ] = beta_theta;
  output[ 4 ] = zeta;
  output[ 5 ] = omega;
  output[ 6 ] = eta_accept;
  output[ 7 ] = cc;
  return output ;
}


// Function :: MCMC algorithm
// Make sure arguments are indexed by iteration
// [[Rcpp::export]]
List zidm(
    int iterations,
    int thin, 
    arma::mat z,
    arma::mat x,
    arma::mat x_theta,
    arma::cube beta_gamma,  
    arma::cube cc,
    arma::vec temp_uu, 
    arma::cube eta, 
    arma::cube eta_accept,
    List eta_index,
    arma::cube omega,
    arma::cube beta_theta,
    arma::cube zeta, 
    arma::cube varphi, 
    double sigma2_beta_gamma,  
    double sigma2_beta_theta
){
  
  // Initiate memory
  List between_beta_gamma_varphi( 3 ); 
  List beta_gamma_out( 2 );
  List between_cc_eta( 3 );
  List between_beta_theta_zeta( 2 );
  int B = z.n_cols;
  int P = varphi.n_cols;
  int P_theta = x_theta.n_cols;
  int N = x.n_rows;   
  arma::mat temp_varphi( B, P );
  arma::mat temp_beta_gamma( B, P ); 
  arma::mat temp_eta( N, B ); 
  arma::mat temp_omega( N, B );
  arma::mat temp_beta_theta( B, P_theta );
  arma::mat temp_zeta( B, P_theta );
  arma::mat temp_eta_accept( B, P );
  arma::mat temp_cc( N, B );
  
  temp_varphi = varphi.slice( 0 );
  temp_beta_gamma = beta_gamma.slice( 0 );  
  temp_eta = eta.slice( 0 );
  temp_omega = omega.slice( 0 );
  temp_beta_theta = beta_theta.slice( 0 );
  temp_zeta = zeta.slice( 0 );
  temp_cc = cc.slice( 0 ); 
  temp_eta_accept = eta_accept.slice( 0 );
  
  arma::mat temp_loggamma( N, B );
  temp_loggamma.zeros();
  for( int n = 0; n < N; ++n ){
    for( int j = 0; j < B; ++j ){
      temp_loggamma( n, j ) = ( x.row( n )*temp_beta_gamma.row( j ).t() ).eval()( 0, 0 ); 
    }
  }
  
  // Looping over the number of iterations specified by user
  for( int iter = 0; iter < iterations; ++iter ){
    
    // Update within
    beta_gamma_out = help::within_beta_gamma_augment_cpp( x, temp_varphi, temp_beta_gamma, temp_loggamma, temp_cc, temp_eta, sigma2_beta_gamma );
    temp_beta_gamma =  as<arma::mat>( beta_gamma_out[ 0 ] );
    temp_loggamma = as<arma::mat>( beta_gamma_out[ 1 ] );

    // Update cc within
    temp_cc = help::cc_update_cpp( z, temp_loggamma, temp_uu, temp_eta );
    
    // Update uu
    temp_uu = help::uu_update_cpp( z, temp_cc );  
    
    // Update cc and eta between
    between_cc_eta = help::between_cc_eta_update_cpp( z, x_theta, temp_beta_theta, temp_cc, temp_eta, eta_index, temp_loggamma, temp_varphi);
    temp_cc = as<arma::mat>( between_cc_eta[ 0 ] );
    temp_eta = as<arma::mat>( between_cc_eta[ 1 ] );
    temp_eta_accept = as<arma::mat>( between_cc_eta[ 2 ] );
    
    // Update omega
    temp_omega = help::update_omega( x_theta, temp_beta_theta );
    
    // Within beta_theta
    temp_beta_theta = help::within_beta_theta( x_theta, temp_zeta, temp_omega, temp_eta, sigma2_beta_theta );
    
    // Set the starting values for the next iteration
    if( (iter + 1) % thin == 0 ){ 
      varphi.slice( (iter + 1)/thin - 1 ) = temp_varphi;
      beta_gamma.slice( (iter + 1)/thin - 1 ) = temp_beta_gamma; 
      eta.slice( (iter + 1)/thin - 1 ) = temp_eta; 
      omega.slice( (iter + 1)/thin - 1 ) = temp_omega;
      beta_theta.slice( (iter + 1)/thin - 1 ) = temp_beta_theta;
      zeta.slice( (iter + 1)/thin - 1 ) = temp_zeta; 
      cc.slice( (iter + 1)/thin - 1 ) = temp_cc; 
      eta_accept.slice( (iter + 1)/thin - 1 ) = temp_eta_accept;   
    }
    
    // Print out progress
    double printer = iter % 250;
    
    if( printer == 0 ){
      Rcpp::Rcout << "Iteration = " << iter << std::endl;
    }
  }
  
  // Return output
  List output( 8 ); 
  output[ 0 ] = varphi;
  output[ 1 ] = beta_gamma; 
  output[ 2 ] = eta;
  output[ 3 ] = beta_theta;
  output[ 4 ] = zeta;
  output[ 5 ] = omega;
  output[ 6 ] = eta_accept;
  output[ 7 ] = cc;
  return output ;
}


// Function :: MCMC algorithm
// Make sure arguments are indexed by iteration
// [[Rcpp::export]]
List dm_bvs(
    int iterations,
    int thin, 
    arma::mat z,
    arma::mat x, 
    arma::cube beta_gamma,  
    arma::cube cc,
    arma::vec temp_uu, 
    arma::cube eta,  
    List eta_index,   
    arma::cube varphi, 
    double sigma2_beta_gamma, 
    double a_varphi, // Beta-Bin
    double b_varphi   // Beta-Bin    
){
  
  // Initiate memory
  List between_beta_gamma_varphi( 3 ); 
  List beta_gamma_out( 2 );
  List between_cc_eta( 3 ); 
  int B = z.n_cols;
  int P = varphi.n_cols; 
  int N = x.n_rows;   
  arma::mat temp_varphi( B, P );
  arma::mat temp_beta_gamma( B, P ); 
  arma::mat temp_eta( N, B ); 
  arma::mat temp_omega( N, B ); 
  arma::mat temp_eta_accept( B, P );
  arma::mat temp_cc( N, B );
  
  temp_varphi = varphi.slice( 0 );
  temp_beta_gamma = beta_gamma.slice( 0 );  
  temp_eta = eta.slice( 0 );   
  temp_cc = cc.slice(0); 
  
  arma::mat temp_loggamma( N, B );
  temp_loggamma.zeros();
  for( int n = 0; n < N; ++n ){
    for( int j = 0; j < B; ++j ){
      temp_loggamma( n, j ) = ( x.row( n )*temp_beta_gamma.row( j ).t() ).eval()( 0, 0 ); 
    }
  }
  
  // Looping over the number of iterations specified by user
  for( int iter = 0; iter < iterations; ++iter ){
    
    // Add/Delete 
    between_beta_gamma_varphi = help::between_beta_gamma_varphi_augment_update_BB_cpp( x, temp_varphi, temp_beta_gamma, temp_loggamma, temp_cc, temp_eta, sigma2_beta_gamma, a_varphi, b_varphi );

    temp_varphi = as<arma::mat>( between_beta_gamma_varphi[ 0 ] );
    temp_beta_gamma = as<arma::mat>( between_beta_gamma_varphi[ 1 ] );
    temp_loggamma = as<arma::mat>( between_beta_gamma_varphi[ 2 ] );
    
    // Update within
    beta_gamma_out = help::within_beta_gamma_augment_cpp( x, temp_varphi, temp_beta_gamma, temp_loggamma, temp_cc, temp_eta, sigma2_beta_gamma );
    temp_beta_gamma =  as<arma::mat>( beta_gamma_out[ 0 ] );
    temp_loggamma = as<arma::mat>( beta_gamma_out[ 1 ] );
    
    // Update cc within
    temp_cc = help::cc_update_cpp( z, temp_loggamma, temp_uu, temp_eta );
    
    // Update uu
    temp_uu = help::uu_update_cpp( z, temp_cc ); 

    
    // Set the starting values for the next iteration
    if( (iter + 1) % thin == 0 ){ 
      varphi.slice( (iter + 1)/thin - 1 ) = temp_varphi;
      beta_gamma.slice( (iter + 1)/thin - 1 ) = temp_beta_gamma; 
      eta.slice( (iter + 1)/thin - 1 ) = temp_eta;  
      cc.slice( (iter + 1)/thin - 1 ) = temp_cc; 
    }
    
    // Print out progress
    double printer = iter % 250;
    
    if( printer == 0 ){
      Rcpp::Rcout << "Iteration = " << iter << std::endl;
    }
  }
  
  // Return output
  List output( 4 ); 
  output[ 0 ] = varphi;
  output[ 1 ] = beta_gamma; 
  output[ 2 ] = eta; 
  output[ 3 ] = cc;
  return output ;
}

// Function :: MCMC algorithm
// Make sure arguments are indexed by iteration
// [[Rcpp::export]]
List dm(
    int iterations,
    int thin, 
    arma::mat z,
    arma::mat x, 
    arma::cube beta_gamma,  
    arma::cube cc,
    arma::vec temp_uu, 
    arma::cube eta,  
    List eta_index, 
    arma::cube varphi, 
    double sigma2_beta_gamma 
){
  
  // Initiate memory
  List between_beta_gamma_varphi( 3 ); 
  List beta_gamma_out( 2 );
  List between_cc_eta( 3 ); 
  int B = z.n_cols;
  int P = varphi.n_cols; 
  int N = x.n_rows;   
  arma::mat temp_varphi( B, P );
  arma::mat temp_beta_gamma( B, P ); 
  arma::mat temp_eta( N, B );  
  arma::mat temp_cc( N, B );
  
  temp_varphi = varphi.slice( 0 );
  temp_beta_gamma = beta_gamma.slice( 0 );  
  temp_eta = eta.slice( 0 ); 
  temp_cc = cc.slice(0);
  
  arma::mat temp_loggamma( N, B );
  temp_loggamma.zeros();
  for( int n = 0; n < N; ++n ){
    for( int j = 0; j < B; ++j ){
      temp_loggamma( n, j ) = ( x.row( n )*temp_beta_gamma.row( j ).t() ).eval()( 0, 0 ); 
    }
  }
  
  // Looping over the number of iterations specified by user
  for( int iter = 0; iter < iterations; ++iter ){
 
    // Update within
    beta_gamma_out = help::within_beta_gamma_augment_cpp( x, temp_varphi, temp_beta_gamma, temp_loggamma, temp_cc, temp_eta, sigma2_beta_gamma );
    temp_beta_gamma =  as<arma::mat>( beta_gamma_out[ 0 ] );
    temp_loggamma = as<arma::mat>( beta_gamma_out[ 1 ] );
    
    // Update cc within
    temp_cc = help::cc_update_cpp( z, temp_loggamma, temp_uu, temp_eta );
    
    // Update uu
    temp_uu = help::uu_update_cpp( z, temp_cc ); 
     
    // Set the starting values for the next iteration
    if( (iter + 1) % thin == 0 ){ 
      varphi.slice( (iter + 1)/thin - 1 ) = temp_varphi;
      beta_gamma.slice( (iter + 1)/thin - 1 ) = temp_beta_gamma; 
      eta.slice( (iter + 1)/thin - 1 ) = temp_eta;  
      cc.slice( (iter + 1)/thin - 1 ) = temp_cc; 
    }
    
    // Print out progress
    double printer = iter % 250;
    
    if( printer == 0 ){
      Rcpp::Rcout << "Iteration = " << iter << std::endl;
    }
  }
  
  // Return output
  List output( 4 ); 
  output[ 0 ] = varphi;
  output[ 1 ] = beta_gamma; 
  output[ 2 ] = eta; 
  output[ 3 ] = cc;
  return output ;
}


// Function :: MCMC algorithm
// Make sure arguments are indexed by iteration
// [[Rcpp::export]]
arma::cube tuyl(
    int iterations,
    arma::mat z
){
  
  int N = z.n_rows; 
  int B = z.n_cols; 
  arma::cube samples = arma::zeros<arma::cube>( N, B, iterations );
 
  for( int i = 0; i < N; ++i ){
    int n = sum( z.row( i ) );
    int zeta = sum( z.row( i ) == 0 ); 
      
    for( int s = 0; s < iterations; ++s ){
         
       arma::vec prob( zeta + 1 );            
       prob.zeros();   
         
        for( int w = 0; w < zeta + 1; ++w ){
          prob[ w ] = Rf_choose( B - w, B - zeta)/Rf_choose( n + B - w - 1, n ); 
        }
         
        prob = prob/sum( prob ); 
  
        double sampled = help::sample_prob_cpp( help::myseq( 0, zeta ), wrap( prob ) ); 
       
         // Vector of zero components
         IntegerVector zero_i( 0 );
         for( int k = 0; k < B; ++k ){
           // Available zeros
           if( z( i, k ) == 0 ){
             zero_i.push_back( k );
           }
         }

         arma::vec dirichlet( B );
         dirichlet.zeros();

         for( int j = 0; j < B; ++j ){
           dirichlet[ j ] = rgamma( 1, z( i, j ) + 1, 1 )[ 0 ];
         }
          
         IntegerVector selected =  help::sample_mult_cpp( zero_i, sampled );
         for( int r = 0; r < sampled; ++r ){
            dirichlet[ selected[ r ] ] = 0;
         }
  
        samples.slice(s).row( i ) = ( dirichlet/sum( dirichlet ) ).t();
        }
    } 
  
  // Return output 
  return samples ;
}

// Function :: MCMC algorithm
// Make sure arguments are indexed by iteration
// [[Rcpp::export]]
arma::mat tuyl_meaner( 
    arma::mat z
){
  
  int N = z.n_rows; 
  int B = z.n_cols; 
  arma::mat estimates( N, B );
  estimates.zeros();
  
  for( int i = 0; i < N; ++i ){
    int n = sum( z.row( i ) );
    int zeta = sum( z.row( i ) == 0 ); 
    
     for( int j = 0; j < B; ++j ){ 
 
           if( z( i, j ) != 0 ){
             double sum = 0;
             for( int k = 0; k < (zeta + 1); ++k ){
               sum = sum + help::Enum( z(i,j), B, zeta, n, k );
             }
             estimates(i,j) = sum;
           }else{
             double sum_zero = 0;
             for( int k = 0; k < (zeta + 1); ++k ){
               sum_zero = sum_zero + help::Ezero( B, zeta, n, k );
             }
             estimates(i,j) = sum_zero;
           }
         }
    }
  
  // Return output 
  return estimates;
}


