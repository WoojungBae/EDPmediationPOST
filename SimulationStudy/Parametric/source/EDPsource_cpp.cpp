// #define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
// #include <vector>
// #include <cmath>
// #include <Rmath.h>
//#include "utilFunctions.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;

template <typename T>
inline bool approx_equal_cpp(const T& lhs, const T& rhs, double tol = 1e-10) {
  return approx_equal(lhs, rhs, "absdiff", tol);
}

umat unique_rows(const umat& M) {
  int n_temp = M.n_rows;
  uvec ind_temp = zeros<uvec>(n_temp);
  
  for (int i = 0; i < n_temp; i++) {
    for (int j = i + 1; j < n_temp; j++) {
      if (approx_equal_cpp(M.row(i), M.row(j))) { 
        ind_temp(j) = 1; 
        break; 
      }
    }
  }
  
  ind_temp = find(ind_temp == 0);
  umat unique_M = M.rows(ind_temp);
  
  uvec unique_M1 = unique_M.col(0);
  uvec unique_M2 = unique_M.col(1);
  
  uvec sorted_unique_M1 = sort(unique_M1);
  uvec sorted_unique_M2 = sort(unique_M2);
  
  uvec M1_temp = unique(sorted_unique_M1);
  int n_temp_C1 = M1_temp.size();
  int count = 0;
  for (int i = 0; i < n_temp_C1; i++) {
    ind_temp = find(unique_M1 == M1_temp(i));
    int n_temp_C2 = ind_temp.size();
    uvec unique_M2_temp = sort(unique_M2(ind_temp));
    
    for (int j = 0; j < n_temp_C2; j++) {
      sorted_unique_M2(count) = unique_M2_temp(j);
      count++;
    }
  }
  
  umat sorted_unique_M = join_rows(sorted_unique_M1, sorted_unique_M2);
  
  return sorted_unique_M;
}

int mod(int const& a, int const& b){
  return a - floor(a/b)*b;
}

mat reorder(mat MAT, int const& k_old, int const& k_new) {
  vec MAT_curr = MAT.col(k_old);
  MAT.insert_cols(k_new,MAT_curr);
  
  if (k_old > k_new) {
    MAT.shed_col(k_old);
  } else {
    MAT.shed_col((k_old+1));
  }
  
  return(MAT);
}

// [[Rcpp::export]]
mat inv_cpp(mat const& MAT) {
  int p_temp = MAT.n_cols;
  int rank_temp = rank(MAT);
  
  mat MATinv;
  if (p_temp > rank_temp){
    MATinv = pinv(MAT);
  } else if (MAT.is_sympd()){
    MATinv = inv_sympd(symmatu(MAT));
  } else {
    MATinv = inv(MAT);
  }
  
  return MATinv;
}

// [[Rcpp::export]]
int rank_cpp(mat const& MAT) {
  return rank(MAT);
}

// [[Rcpp::export]]
vec invlogit_cpp(vec const& x_temp) { 
  return((1/(1+exp(-x_temp))));
}

// [[Rcpp::export]]
vec count_cpp(vec const& VEC, int const& p) {
  vec count(p);
  for (int q = 0; q < p; q++) {
    uvec index = find(VEC == q);
    count(q) = index.n_elem; 
  }
  
  return(count);
}

vec IQRoutlier_cpp(vec const& x_temp) { 
  vec Q_temp = {0.25,0.5,0.75};
  vec Qx_temp = quantile(x_temp, Q_temp);
  
  double med_temp = Qx_temp[1];
  double IQR_temp = Qx_temp[2] - Qx_temp[0];
  vec IQRoutlier = {med_temp - 1.5*IQR_temp, med_temp + 1.5*IQR_temp};
  return(IQRoutlier);
}

// [[Rcpp::export]]
vec rbeta_cpp(int const& N, double const& s, double const& r) {
  vec rgamma1 = randg(N, distr_param(s,1.0));
  vec rgamma2 = randg(N, distr_param(r,1.0));
  vec rbeta = rgamma1/(rgamma1+rgamma2);
  
  return(rbeta);
}

// [[Rcpp::export]]
double mvbeta_cpp(vec const& VEC, bool const& logt) {
  double mvbeta = sum(lgamma(VEC))-lgamma(sum(VEC));
  if (logt == false) {mvbeta = exp(mvbeta);}
  
  return(mvbeta);
}

// [[Rcpp::export]]
vec factorial_cpp(vec const& N, bool const& logt) {
  int p = N.size();
  
  vec factorial = zeros(p);
  for (int q = 0; q < p; q++) {
    int N_temp = N(q);
    if (N_temp > 0){
      for (int i = 1; i <= N_temp; i++) {
        factorial(q) += log(i);
      }
    }
  }
  
  if (logt == false) {factorial = exp(factorial);}
  
  return(factorial);
}

int rmultinom_cpp(vec const& prob) {
  int p = prob.size();
  
  vec cumprob = cumsum(prob/sum(prob));
  
  int rmultinom = 1;
  double rnd = randu(1)[0];
  for (int q = 0; q < p; q++) {
    if (rnd>cumprob(q)) {
      rmultinom = rmultinom+1;
    }
  }
  
  return(rmultinom);
}

int dmultinom_cpp(vec const& X, vec const& prob, bool const& logt) {
  vec N(1);
  N(0) = X.size();
  vec log_prob = log(prob);
  log_prob.replace(-datum::inf,0);
  
  double dmultinom = factorial_cpp(N,true)[0] - sum(factorial_cpp(X,true) + X % log_prob);
  if (logt == false) {dmultinom = exp(dmultinom);}
  
  return(dmultinom);
}

// [[Rcpp::export]]
mat rdirichlet_cpp(int const& N, vec const& alpha) {
  int p = alpha.size();
  
  mat rdirichlet(N,p);
  for (int q = 0; q < p; q++) {
    rdirichlet.col(q) = randg(N,distr_param(alpha(q),1.0));
  }
  vec sum_rdirichlet = sum(rdirichlet,1);
  rdirichlet.each_col() /= sum_rdirichlet;
  
  return rdirichlet.t();
}

// [[Rcpp::export]]
double ddirichlet_cpp(vec const& X, vec const& alpha, bool const& logt) {
  int p = alpha.size();
  vec ones_p = ones(p);
  
  double ddirichlet = lgamma(sum(alpha)) - sum(lgamma(alpha)) + sum((alpha - ones_p) % log(X));
  if (logt == false) {ddirichlet = exp(ddirichlet);}
  
  return ddirichlet;
}

// [[Rcpp::export]]
vec rscainvchisq_cpp(int const& N, double const& a_tau2, double const& b_tau2) {
  vec ones_N = ones(N);
  vec ab_tau_N = a_tau2 * b_tau2 * ones_N;
  
  vec rtau = ab_tau_N/chi2rnd(a_tau2,N);
  return rtau;
}

// [[Rcpp::export]]
double dscainvchisq_cpp(double const& tau, double const& a_tau2, double const& b_tau2, bool const& logt) {
  double ab_tau = a_tau2 * b_tau2;
  double a_tau_half = a_tau2/2;
  double ab_tau_half = ab_tau/2;
  
  double dtau = a_tau_half * log(ab_tau_half) - ab_tau_half/tau - 
    lgamma(a_tau_half) - (1+a_tau_half) * log(tau);
  if (logt == false) {dtau = exp(dtau);}
  
  return dtau;
}

// [[Rcpp::export]]
mat rmvn_cpp(int const& N, vec const& MU, mat const& SIG) {
  return mvnrnd(MU, SIG, N);
}

// [[Rcpp::export]]
double dmvn_cpp(vec const& X_temp, vec const& MU, mat const& SIG, bool const& logt) {
  int p = SIG.n_cols;
  mat Chol_Sig = chol(SIG);
  vec std_temp = solve(Chol_Sig, X_temp - MU); //check
  double std2_temp = sum(std_temp % std_temp);
  
  double dmvn = - (sum(log(Chol_Sig.diag())) + p*log(2*M_PI)/2 + std2_temp/2);
  if (logt == false) {dmvn = exp(dmvn);}
  
  return dmvn;
}

// [[Rcpp::export]]
vec rt_cpp(int const& N, double const& MU, double const& SIG2, int const& nu) {
  vec U = sqrt(nu/chi2rnd(nu, N));
  vec rt = sqrt(SIG2) * randn(N) % U  + MU * ones(N);
  
  return rt;
}

// [[Rcpp::export]]
mat rmvt_cpp(int const& N, vec const& MU, mat const& SIG, int const& nu) {
  int p = MU.size();
  
  vec U = sqrt(nu/chi2rnd(nu, N));
  mat Y = mvnrnd(zeros(p), SIG, N);
  
  mat rmvt(p,N);
  for(int i = 0; i < N; i++) {
    rmvt.col(i) = U(i) * Y.col(i) + MU;
  }
  
  return rmvt;
}

// [[Rcpp::export]]
double dmvt_cpp(vec const& X_temp, vec const& MU, mat const& SIG, int const& nu, bool logt) {
  int p = MU.size();
  
  double dmvt = (lgamma((nu+p)/2)-(lgamma(nu/2)+(p/2)*log(nu*M_PI)+log_det_sympd(SIG)/2))-
    ((nu+p)/2)*log(1+as_scalar((X_temp-MU).t()*inv_cpp(SIG)*(X_temp-MU))/nu);
  if(logt == false) {dmvt = exp(dmvt);}
  
  return dmvt;
}

// Define function to update alpha_theta prior Gamma(a_theta,b_theta) with shape and rate parameters
// [[Rcpp::export]]
double update_alpha_theta_cpp(int const& N, int const& Ky, double const& alpha_theta_curr,
                              double const& a_theta, double const& b_theta) {
  double logxi = log(R::rbeta(alpha_theta_curr + 1, N));
  double a_theta_new = a_theta + Ky - 1;
  double b_theta_new = b_theta - logxi;
  double binv_theta_new = 1/b_theta_new;
  
  double varrho = a_theta_new/(N * b_theta_new + a_theta_new);
  varrho = R::rbinom(1, varrho);
  
  // randg: shape and scale at default in Rcpp & RcppArmadillo
  // rgamma: shape and rate at default in R
  double alpha_theta_prop;
  if (varrho > 0){
    alpha_theta_prop = randg(distr_param(a_theta_new + 1, binv_theta_new));
  } else {
    alpha_theta_prop = randg(distr_param(a_theta_new, binv_theta_new));
  }
  
  return alpha_theta_prop;
}

// Define function to update alpha_omega using Metropolis Hastings
// [[Rcpp::export]]
double update_alpha_omega_cpp(int const& Kyx, int const& Ky, uvec const& Sy, 
                              uvec const& unique_Sy, double const& alpha_omega_curr,
                              double const& a_omega, double const& b_omega) {
  double bivn_omega = 1/b_omega;
  double alpha_omega_prop = randg(distr_param(a_omega,bivn_omega));
  
  int n_temp = 0;
  double prod_prop = 1;
  double prod_curr = 1;
  for (int k = 0; k < Ky; k++) {
    uvec ind_Sy = find(Sy == unique_Sy(k));
    n_temp = ind_Sy.n_elem; 
    
    if (n_temp>0){
      prod_prop = prod_prop * (alpha_omega_prop + n_temp) * R::beta(alpha_omega_prop + 1, n_temp);
      prod_curr = prod_curr * (alpha_omega_curr + n_temp) * R::beta(alpha_omega_curr + 1, n_temp);
    }
  }
  
  double like_prop = R::dgamma(alpha_omega_prop, a_omega, bivn_omega, false) *
    (pow(alpha_omega_prop, (Kyx - Ky))) * prod_prop;
  double like_curr = R::dgamma(alpha_omega_curr, a_omega, bivn_omega, false) *
    (pow(alpha_omega_curr, (Kyx - Ky))) * prod_curr;
  
  double ratio = (like_prop)/(like_curr);
  if (ratio < randu(1)[0]){ alpha_omega_prop = alpha_omega_curr; }
  
  return alpha_omega_prop;
}