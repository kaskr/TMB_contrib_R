
#include <TMB.hpp>

template<class Type>
Type dpoislinkgamma( Type x, Type n, Type w, Type cv, int give_log=0){
  Type dens;
  Type enc_prob = 1 - exp(-n);
  Type posmean = n * w / enc_prob;
  if( x==0 ){
    dens = 1 - enc_prob;
  }else{
    dens = enc_prob * dgamma(x, pow(cv,-2), posmean*pow(cv,2));
  }
  if(give_log) return log(dens); else return dens;
}
template<class Type>
Type ppoislinkgamma( Type x, Type n, Type w, Type cv){
  Type enc_prob = 1 - exp(-n);
  Type posmean = n * w / enc_prob;
  Type dist = 1 - enc_prob;
  if( x>0 ){
    dist += enc_prob * pgamma(x, pow(cv,-2), posmean*pow(cv,2));
  }
  return dist;
}


// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;

  // Data
  DATA_VECTOR( c_i );  // counts for observation i
  DATA_FACTOR( j_i );  // Random effect index for observation i

  // SPDE objects
  DATA_SPARSE_MATRIX(M0);
  DATA_SPARSE_MATRIX(M1);
  DATA_SPARSE_MATRIX(M2);
  DATA_VECTOR_INDICATOR(keep, c_i);

  // Parameters
  PARAMETER( beta0 );
  PARAMETER( ln_tau );
  PARAMETER( ln_kappa );
  PARAMETER( ln_w );
  PARAMETER( ln_cv );

  // Random effects
  PARAMETER_VECTOR( epsilon_j );

  // Objective funcction
  Type jnll = 0;

  // Derived quantities
  Type Range = sqrt(8) / exp( ln_kappa );
  Type SigmaE = 1 / sqrt(4 * M_PI * exp(2*ln_tau) * exp(2*ln_kappa));

  // Probability of random effects
  Eigen::SparseMatrix<Type> Q = exp(4*ln_kappa)*M0 + Type(2.0)*exp(2*ln_kappa)*M1 + M2;
  jnll += SCALE( GMRF(Q), 1/exp(ln_tau) )( epsilon_j );

  // Probability of data conditional on random effects
  vector<Type> cdf_i(c_i.size());
  vector<Type> mu_i(c_i.size());
  for( int i=0; i<c_i.size(); i++){
    mu_i(i) = exp( beta0 + epsilon_j(j_i(i)) );
    // oneStePredict using CDF
    cdf_i(i) = squeeze( ppoislinkgamma( c_i(i), mu_i(i), exp(ln_w), exp(ln_cv) ) );
    jnll -= keep.cdf_lower(i) * log( cdf_i(i) );
    jnll -= keep.cdf_upper(i) * log( 1.0 - cdf_i(i) );
    // oneStepPredict using oneStepGeneric
    jnll -= keep(i) * dpoislinkgamma( c_i(i), mu_i(i), exp(ln_w), exp(ln_cv), true );
  }

  // Reporting
  REPORT( jnll );
  REPORT( cdf_i );
  REPORT( mu_i )
  REPORT( Range );
  REPORT( SigmaE );

  return jnll;
}
