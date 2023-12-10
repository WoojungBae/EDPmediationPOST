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
  if(varrho > 0){
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
  for(int k = 0; k < Ky; k++) {
    uvec ind_Sy = find(Sy == unique_Sy(k));
    n_temp = ind_Sy.n_elem; 
    
    if(n_temp>0){
      prod_prop = prod_prop * (alpha_omega_prop + n_temp) * R::beta(alpha_omega_prop + 1, n_temp);
      prod_curr = prod_curr * (alpha_omega_curr + n_temp) * R::beta(alpha_omega_curr + 1, n_temp);
    }
  }
  
  double like_prop = R::dgamma(alpha_omega_prop, a_omega, bivn_omega, false) *
    (pow(alpha_omega_prop, (Kyx - Ky))) * prod_prop;
  double like_curr = R::dgamma(alpha_omega_curr, a_omega, bivn_omega, false) *
    (pow(alpha_omega_curr, (Kyx - Ky))) * prod_curr;
  
  double ratio = (like_prop)/(like_curr);
  if(ratio < randu(1)[0]){ alpha_omega_prop = alpha_omega_curr; }
  
  return alpha_omega_prop;
}

/* FOR BIVARIATE confounders */
double update_conf_bin_cpp(double const& X_temp, double const& a_bin, double const& b_bin) {
  
  double a_bin_new = a_bin + X_temp;
  double b_bin_new = b_bin + 1 - X_temp;
  
  double bin_prop = R::rbeta(a_bin_new, b_bin_new);
  
  return bin_prop;
}

/* FOR CONTINUOUS confounders */
List update_conf_con_cpp(double const& X_temp, double const& a_var, double const& b_var,
                          double const& a_mean, double const& b_mean) {
  double ab_var = a_var * b_var;
  double b_mean_temp = b_mean/(b_mean + 1);
  double sum2_temp = pow((X_temp - a_mean),2);
  
  // Update the continuous covariate variance
  double a_var_new = a_var + 1;
  double b_var_new = (ab_var + b_mean_temp * sum2_temp)/a_var_new; 
  double var_prop = rscainvchisq_cpp(1, a_var_new, b_var_new)[0];
  
  // Update the continuous covariate mean
  double b_mean_new = b_mean + 1;
  double a_mean_new = (a_mean * b_mean + X_temp)/b_mean_new; 
  double var_new = var_prop/b_mean_new;
  double mean_prop = sqrt(var_new) * randn(1)[0] + a_mean_new;
  
  return List::create(_["var_par"]  = var_prop,
                      _["mean_par"] = mean_prop);
}

/* FOR CONTINUOUS OUTCOME */
List update_reg_con_cpp(double const& Y_temp, mat const& matX_temp,
                        double const& a_sig2, double const& b_sig2,
                        vec const& a_beta, mat const& B_beta,
                        mat const& Binv_beta, mat const& aBinv_beta) {
  
  // Set initial value
  vec beta_curr = mvnrnd(a_beta, B_beta);
  
  // Update variance
  vec error_temp = Y_temp - (matX_temp * beta_curr);
  double error2_temp = dot(error_temp,error_temp);
  
  double a_sig2_new = (a_sig2 + 1);
  double b_sig2_new = (a_sig2 * b_sig2 + error2_temp)/a_sig2_new;
  double sig2_prop = rscainvchisq_cpp(1, a_sig2_new, b_sig2_new)[0];
  
  // Update coefficient
  mat Binv_beta_new = (matX_temp.t() * matX_temp)/sig2_prop + Binv_beta;
  mat B_beta_new = inv_cpp(Binv_beta_new);
  vec a_beta_new = B_beta_new * ((matX_temp.t() * Y_temp)/sig2_prop + aBinv_beta);
  vec beta_prop = mvnrnd(a_beta_new, B_beta_new);
  
  return List::create(_["sig2_par"] = sig2_prop,
                      _["beta_par"] = beta_prop);
}

/* FOR BINARY OUTCOME */
vec update_reg_bin_cpp(int const& Y_temp, vec const& matX_temp, vec const& beta_curr,
                       vec const& a_beta, mat const& B_beta) {
  vec beta_prop = mvnrnd(beta_curr, B_beta);
  
  //careful this is a sum in original
  double loglik_prop = invlogit_cpp(matX_temp*beta_prop)[0] + 
    dmvn_cpp(beta_prop,a_beta,B_beta,true);
  double loglik_curr = invlogit_cpp(matX_temp*beta_curr)[0] + 
    dmvn_cpp(beta_curr,a_beta,B_beta,true);
  
  double ratio = exp(loglik_prop-loglik_curr);
  if (ratio < randu(1)[0]){ beta_prop = beta_curr; }
  
  return beta_prop;
}


// [[Rcpp::export]]
List EDPcluster_concon(int const& N, int const& py, int const& pm, 
                       int const& p_Z, int const& p_C1, int const& p_C2,
                       vec const& Y, vec const& M, mat const& X, 
                       mat const& matM, mat const& matX,
                       uvec Sy, uvec Sx, umat unique_Syx, 
                       mat YbetaPars, mat Ysig2Pars, mat MbetaPars, mat Msig2Pars,
                       mat ZpiPars, mat XpiPars, mat XmuPars, mat Xtau2Pars,
                       double const& alpha_theta, double const& alpha_omega,
                       vec const& f0_y, vec const& f0_m, vec const& f0_x,
                       double const& a_Ysig2, double const& b_Ysig2,
                       vec const& a_Ybeta, mat const& B_Ybeta, mat const& Binv_Ybeta, 
                       mat const& aBinv_Ybeta, double const& aaBinv_Ybeta,
                       double const& a_Msig2, double const& b_Msig2, 
                       vec const& a_Mbeta, mat const& B_Mbeta, mat const& Binv_Mbeta, 
                       mat const& aBinv_Mbeta, double const& aaBinv_Mbeta,
                       double const& a_pi, double const& b_pi,
                       double const& a_mu, double const& b_mu, 
                       double const& a_tau2, double const& b_tau2) {
  /* 
   * Y is nx1 vector of outcomes
   * M is nx1 vector of mediations
   * X is nxp matrix of confounders
   * Sy is nx1 vector of Y cluster memberships
   * Sx is nx1 vector of X cluster memberships
   * YbetaPars is kxp matrix of coef (k is number of Y clusters)
   * XpiPars is k*xp_C1 matrix of coef (k* is total # of clusters, p_C1 is # binary)
   * XmuPars is k*xp_C2 matrix of coef (p_C2 is # continuous covs)
   * Xtau2Pars is k*xp_C2 matrix of coef
   * alpha_omega is the current value of alpha_{omega}
   * alpha_theta is the current value of alpha_{theta}
   * f0_y is nx1 vector
   * f0_x is nx1 vector
   * unique_Syx is k*x2 matrix of ordered unique observations of Sy and Sx
   * indS is nx1 indicator of which row of unique_Syx everyone is in
   */
  
  // how many unique clusters
  int Kyx; // unique total clusters
  int Ky;      // unique Y-clusters
  
  // initialize some dummy vars to store things
  int num_Sy, num_Sy0, num_Syx, num_unique_Syx;
  uvec ind_Sy, ind_Sy0, ind_Syx, ind_unique_Syx;
  double Ylikreg, Mlikreg, Xlik;
  
  // containers for output
  int k_new;
  
  vec ZpiPars_prop(p_Z);
  vec XpiPars_prop(p_C1);
  vec XmuPars_prop(p_C2);
  vec Xtau2Pars_prop(p_C2);
  
  List outcome_list, mediator_list, convar_con_list;
  
  double Zpi_prop, Xpi_prop, Xmu_prop, Xtau2_prop;
  vec Ybeta_prop, Mbeta_prop, Ysig_prop, Msig_prop;
  double Msig2_prop, Ysig2_prop;
  vec Ybeta_curr, Mbeta_curr, Ysig_curr, Msig_curr;
  
  double Y_temp, M_temp;
  rowvec matM_temp, matX_temp, X_temp;
  uword Sy_temp, Sx_temp;
  
  // LOOP THROUGH EVERY PERSON AND CHANGE CLUSTER MEMBERSHIPS
  for(int i = 0; i < N; i++) {
    matM_temp = matM.row(i);
    matX_temp = matX.row(i);
    X_temp = X.row(i);
    
    Y_temp = Y(i);
    M_temp = M(i);
    
    Sy_temp = Sy(i);
    Sx_temp = Sx(i);
    
    // ---------------------------------------------------------------------------------------
    // ---------------------- Delete i-th person / Adjust the membershop ---------------------
    // ---------------------------------------------------------------------------------------
    // Rprintf("Loop: %d\n",i);
    // check if ith person is the lone person in his cluster
    ind_Syx = find(Sy == Sy_temp && Sx == Sx_temp);
    num_Syx = ind_Syx.size(); 
    // Rprintf("Cluster Y,X,#:%d,%d,%d\n",Sy_temp,Sx_temp,num_Syx);
    
    if(num_Syx==1) { //if lone person in X-Y cluster
      // DELETE ASSOCIATED COEFFICIENTS IN Y AND X CLUSTER
      
      ind_Sy = find(Sy == Sy_temp); //check if only person in Y cluster too
      num_Sy = ind_Sy.size();
      
      // Rprintf("num_Sy: %d\n",num_Sy);
      // Rprintf("Sy_temp: %d\n", Sy_temp);
      // Rcout << Sy.t() << std::endl;
      // Rcout << ind_Sy.t() << std::endl;
      // Rcout << ind_Sy.size() << std::endl;
      // delete Y coef if only one in Y cluster
      if(num_Sy==1) {
        YbetaPars.shed_col(Sy_temp-1); // NEW FOR OUTCOME
        Ysig2Pars.shed_col(Sy_temp-1);
        MbetaPars.shed_col(Sy_temp-1); // NEW FOR MEDIATION
        Msig2Pars.shed_col(Sy_temp-1);
      }
      // delete X coef
      // should find row in unique_Syx that corresponds to person i
      ind_unique_Syx = find(unique_Syx.col(0)==Sy_temp && unique_Syx.col(1)==Sx_temp);
      int unique_Syx_n_rows = unique_Syx.n_rows;
      
      // CHECK --- NEEDED TO DO THIS BECAUSE ind_unique_Syx is uvec
      ZpiPars.shed_col(ind_unique_Syx(0));
      XpiPars.shed_col(ind_unique_Syx(0)); 
      XmuPars.shed_col(ind_unique_Syx(0));
      Xtau2Pars.shed_col(ind_unique_Syx(0)); 
      
      // relabel X cluster
      for(int ii = 0; ii < N; ii++) {
        if(Sy(ii) == Sy_temp && Sx(ii) > Sx_temp){
          Sx(ii) = Sx(ii) - 1;
        }
      }
      
      for(int k = 0; k < unique_Syx_n_rows; k++) {
        if(unique_Syx(k,0) == Sy_temp && unique_Syx(k,1) > Sx_temp){
          unique_Syx(k,1) = unique_Syx(k,1) - 1;
        }
      }
      
      // relabel Y cluster (if needed)
      if(num_Sy == 1) {
        for(int ii = 0; ii < N; ii++) {
          if(Sy(ii) > Sy_temp){
            Sy(ii) = Sy(ii) - 1;
          }
        }
        
        for(int k = 0; k < unique_Syx_n_rows; k++) {
          if(unique_Syx(k,0) > Sy_temp){
            unique_Syx(k,0) = unique_Syx(k,0) - 1;
          }
        }
      }
      
      unique_Syx.shed_row(ind_unique_Syx(0)); // get rid of row
    }
    
    // NEED TO DELETE ROW OF Sy and Sx
    Sy.shed_row(i); Sx.shed_row(i);
    
    // ---------------------------------------------------------------------------------------
    // ------------------- Calculate the draw probabilty for a new cluster -------------------
    // ---------------------------------------------------------------------------------------
    // recalculate number of unique clusters
    Kyx = ZpiPars.n_cols;
    Ky = YbetaPars.n_cols;
    int K_totalposs = Kyx + Ky + 1;
    // Rprintf("Y clusters, total clusters, Total possible clusters: %d,%d,%d\n",Ky,Kyx,K_totalposs);
    
    int Kx_Sy; // # of X clusters within each Y cluster
    int n_k_woi, n_rk_woi; // counts for # in appropriate Y and X cluster, excluding the i-th person
    
    vec probs(K_totalposs);
    int count = 0;    
    for(int k = 0; k < Ky; k++) {
      // ---------------------------------------------------------------------------------------
      // ------------------------------- Existing y&m-x clusters -------------------------------
      // ---------------------------------------------------------------------------------------
      
      // get count of number of X clusters within k-th Y cluster
      ind_Sy0 = find(unique_Syx.col(0) == (k+1));
      Kx_Sy = ind_Sy0.size();
      // Rprintf("NumX cluster within Y cluster: %d\n",Kx_Sy);
      
      // get number of subjects within k-th cluster
      ind_Sy = find(Sy == (k+1));
      n_k_woi = ind_Sy.size();
      
      // Rprintf("n_k_woi: %d\n",n_k_woi);
      // likelihood for each existing Y cluster
      // Ylikreg = R::dbinom(Y_temp, 1 , invlogit_cpp(matX.row(i),YbetaPars.col(k))[0], false);
      Ylikreg = normpdf(Y_temp,dot(matM_temp,YbetaPars.col(k)),sqrt(Ysig2Pars(k)));
      Mlikreg = normpdf(M_temp,dot(matX_temp,MbetaPars.col(k)),sqrt(Msig2Pars(k)));
      
      // Rprintf("likelihood for Y: %.2f\n",Ylikreg);
      for(int r = 0; r < Kx_Sy; r++) {
        ind_Syx = find(Sy == (k+1) && Sx == (r+1));
        n_rk_woi = ind_Syx.size();
        // Rprintf("n_rk_woi: %d\n",n_rk_woi);
        
        Xlik = 1;
        // likelihood for binary confounders
        for (int q = 0; q < p_Z; q++) {
          Xlik *= R::dbinom(X_temp(q), 1, ZpiPars(q,count), false);
        }
        // likelihood for binary confounders
        for(int q = 0; q < p_C1; q++) {
          Xlik = Xlik * R::dbinom(X_temp(p_Z+q), 1, XpiPars(q,count), false);
        }
        // likelihood for continuous confounders
        for(int q = 0; q < p_C2; q++) {
          Xlik = Xlik * normpdf(X_temp(p_Z+p_C1+q), XmuPars(q,count), sqrt(Xtau2Pars(q,count)));
        }
        // Rprintf("Prodx2: %.2f\n",Xlik);
        
        probs(count) = (n_k_woi * (n_rk_woi/(n_k_woi+alpha_omega))) * Ylikreg * Mlikreg * Xlik;
        // Rprintf("\t current count and prob: %d, %.2f\n",count,probs(count));
        count++;
      }
      
      // ---------------------------------------------------------------------------------------
      // --------------------------- Existing y&m, but New x cluster ---------------------------
      // ---------------------------------------------------------------------------------------
      probs(Kyx+k) = (n_k_woi * (alpha_omega/(n_k_woi+alpha_omega))) * Ylikreg * Mlikreg * f0_x(i);
      // probs(Kyx+k) = ((n_k_woi * alpha_omega)/(n_k_woi+alpha_omega)) * Ylikreg * Mlikreg * f0_x_i;
    }
    
    // ---------------------------------------------------------------------------------------
    // ---------------------------------- New y&m-x clusters ---------------------------------
    // ---------------------------------------------------------------------------------------
    probs(Kyx+Ky) = alpha_theta * f0_y(i) * f0_m(i) * f0_x(i);
    // probs(Kyx+Ky) = alpha_theta * f0_y_i * f0_m_i * f0_x_i; 
    
    // ---------------------------------------------------------------------------------------
    // ------------------ USE MULTINOMIAL DISTRIBUTION TO CHOOSE NEW CLUSTER -----------------
    // ---------------------------------------------------------------------------------------
    k_new = rmultinom_cpp(probs);
    probs.zeros();
    
    // Rprintf("The new cluster is: %d\n", k_new);
    // need to map this integer to one of the clusters (or a new cluster)
    // ---------------------------------------------------------------------------------------
    // ------------------------------- Existing y&m-x clusters -------------------------------
    // ---------------------------------------------------------------------------------------
    if(k_new<=Kyx) {
      Sy.insert_rows(i,1);
      Sy(i) = unique_Syx(k_new-1,0);
      Sx.insert_rows(i,1);
      Sx(i) = unique_Syx(k_new-1,1);
    } else {
      // find out whether this is a new Y cluster or X cluster
      // ---------------------------------------------------------------------------------------
      // ---------------------------------- New y&m-x clusters ---------------------------------
      // ---------------------------------------------------------------------------------------
      if(k_new == K_totalposs) {
        Sx.insert_rows(i,1);
        Sx(i) = 1;
        
        // ---------------------------------------------------------------------------------------
        // --------------------------------- Draw New parameters ---------------------------------
        // ---------------------------------------------------------------------------------------
        // need functions update_reg_con_cpp
        outcome_list = update_reg_con_cpp(Y_temp,matM_temp,a_Ysig2,b_Ysig2,a_Ybeta,B_Ybeta,Binv_Ybeta,aBinv_Ybeta);
        Ysig2_prop = outcome_list["sig2_par"];
        Ybeta_prop = as<arma::vec>(outcome_list["beta_par"]);
        
        // need functions update_reg_con_cpp
        mediator_list = update_reg_con_cpp(M_temp,matX_temp,a_Msig2,b_Msig2,a_Mbeta,B_Mbeta,Binv_Mbeta,aBinv_Mbeta);
        Msig2_prop = mediator_list["sig2_par"];
        Mbeta_prop = as<arma::vec>(mediator_list["beta_par"]);
        
        // need functions update_conf_bin_cpp, update_conf_con_cpp
        for (int q = 0; q<p_Z; q++) {
          Zpi_prop = update_conf_bin_cpp(X_temp(q),a_pi,b_pi);
          ZpiPars_prop(q) = Zpi_prop;
        }
        
        for(int q = 0; q<p_C1; q++) {
          Xpi_prop = update_conf_bin_cpp(X_temp(p_Z+q),a_pi,b_pi);
          XpiPars_prop(q) = Xpi_prop;
        }
        
        for(int q = 0; q < p_C2; q++) {
          convar_con_list = update_conf_con_cpp(X_temp(p_Z+p_C1+q),a_tau2,b_tau2,a_mu,b_mu);
          Xtau2_prop = convar_con_list["var_par"];
          Xtau2Pars_prop(q) = Xtau2_prop;
          
          Xmu_prop = convar_con_list["mean_par"];
          XmuPars_prop(q) = Xmu_prop;
        }
        
        // ---------------------------------------------------------------------------------------
        // ----------------------------------- Save parameters -----------------------------------
        // ---------------------------------------------------------------------------------------
        YbetaPars.insert_cols(Ky,Ybeta_prop);
        Ysig2Pars.insert_cols(Ky,1);
        Ysig2Pars(Ky) = Ysig2_prop;
        
        MbetaPars.insert_cols(Ky,Mbeta_prop);
        Msig2Pars.insert_cols(Ky,1);
        Msig2Pars(Ky) = Msig2_prop;
        
        ZpiPars.insert_cols(Kyx,ZpiPars_prop);
        XpiPars.insert_cols(Kyx,XpiPars_prop);
        XmuPars.insert_cols(Kyx,XmuPars_prop);
        Xtau2Pars.insert_cols(Kyx,Xtau2Pars_prop);
        
        Sy.insert_rows(i,1);
        Sy(i) = Sy.max()+1;
        unique_Syx.insert_rows(Kyx,1);
        unique_Syx(Kyx,0) = Sy(i);
        unique_Syx(Kyx,1) = Sx(i);
      } else {
        // ---------------------------------------------------------------------------------------
        // --------------------------- Existing y&m, but New x cluster ---------------------------
        // ---------------------------------------------------------------------------------------
        // Y cluster should be k_new - Kyx
        // X cluster should be max + 1 of X clusters within Y
        
        Sy.insert_rows(i,1);
        Sy(i) = k_new - Kyx;
        ind_Sy0 = find(unique_Syx.col(0) == Sy(i));
        Sx.insert_rows(i,1);
        Sx(i) = ind_Sy0.size() + 1;
        
        // ---------------------------------------------------------------------------------------
        // --------------------------------- Draw New parameters ---------------------------------
        // ---------------------------------------------------------------------------------------
        // need functions update_conf_bin_cpp, update_conf_con_cpp
        for (int q = 0; q<p_Z; q++) {
          Zpi_prop = update_conf_bin_cpp(X_temp(q),a_pi,b_pi);
          ZpiPars_prop(q) = Zpi_prop;
        }
        
        for (int q = 0; q<p_C1; q++) {
          Xpi_prop = update_conf_bin_cpp(X_temp(p_Z+q),a_pi,b_pi);
          XpiPars_prop(q) = Xpi_prop;
        }
        
        for(int q = 0; q < p_C2; q++) {
          convar_con_list = update_conf_con_cpp(X_temp(p_Z+p_C1+q),a_tau2,b_tau2,a_mu,b_mu);
          Xtau2_prop = convar_con_list["var_par"];
          Xtau2Pars_prop(q) = Xtau2_prop;
          
          Xmu_prop = convar_con_list["mean_par"];
          XmuPars_prop(q) = Xmu_prop;
        }
        
        // ---------------------------------------------------------------------------------------
        // ----------------------------------- Save parameters -----------------------------------
        // ---------------------------------------------------------------------------------------
        ind_Sy0 = find(unique_Syx.col(0) <= Sy(i));
        num_Sy0 = ind_Sy0.size();
        
        ZpiPars.insert_cols(num_Sy0,ZpiPars_prop);
        XpiPars.insert_cols(num_Sy0,XpiPars_prop);
        XmuPars.insert_cols(num_Sy0,XmuPars_prop);
        Xtau2Pars.insert_cols(num_Sy0,Xtau2Pars_prop);
        
        unique_Syx.insert_rows(num_Sy0,1);
        unique_Syx(num_Sy0,0) = Sy(i);
        unique_Syx(num_Sy0,1) = Sx(i);
      }
    }
    
    // as necessary add back into Sy, Sx, unique_Syx, etc...
    // Rprintf("New cluster for person i: %d, %d\n\n",Sy(i),Sx(i));
    // Rcout << unique_Syx << std::endl;
    // Rcout << YbetaPars << std::endl;
    // Rprintf("Max Sy, Max unique_Syx: %d,%d\n",Sy.max(),unique_Syx.col(0).max());
  }

  uvec unique_Sy = unique(Sy);  
  umat Syx = join_rows(Sy, Sx);
  unique_Syx = unique_rows(unique_Syx);
  
  // ------------------------------------------------------------------------------
  // -------------------- Improve Mixing (Metropolis-Hastings) --------------------
  // ------------------------------------------------------------------------------
  Ky = unique_Sy.size();
  Kyx = unique_Syx.n_rows;
  int max_Sx = Sx.max();
  
  vec n_k = zeros(Ky);
  mat n_rk = zeros(Ky,max_Sx);
  vec max_Kx_Sy = zeros(Ky);
  for (int k = 0; k < Ky; k++) {
    ind_Sy = find(Sy == (k+1));
    num_Sy = ind_Sy.size();
    n_k(k) = num_Sy;
    
    ind_Sy0 = find(unique_Syx.col(0) == (k+1));
    num_Sy0 = ind_Sy0.size();
    max_Kx_Sy(k) = num_Sy0;
    for(int r = 0; r < num_Sy0; r++) {
      ind_Syx = find(Sy == (k+1) && Sx == (r+1));
      num_Syx = ind_Syx.size();
      n_rk(k,r) = num_Syx;
    }
  }
  
  // initialize some dummy vars to store things
  double ratio;
  
  int num_max_Kx_Sy, num_max_Kx_Sy_temp;
  vec max_Kx_Sy_temp, cumsum_max_Kx_Sy;
  
  uvec ind_max_Kx_Sy, ind_max_Kx_Sy_temp, uvec_cumsum_max_Kx_Sy, unique_Sy_temp;
  umat unique_Syx_temp;
  
  int n_k_curr, n_k_prop, n_rk_curr;
  int Kx1_curr, Kx1_prop, Kx2_curr, Kx2_prop;
  int k_curr, r_curr, k_prop, r_prop;
  int k_curr_1, r_curr_1, k_prop_1;
  
  // ---------------------------------- 1st move ----------------------------------
  if(Kyx>Ky && Ky>1){
    ind_max_Kx_Sy = find(max_Kx_Sy == 1);
    num_max_Kx_Sy = ind_max_Kx_Sy.size();
    
    Kx2_curr = Kyx - num_max_Kx_Sy;
    
    // Choose cluster (k,r) and h
    int choose_kr_temp = rmultinom_cpp(ones(Kx2_curr));
    int choose_new_k_temp = rmultinom_cpp(ones(Ky-1));
    
    unique_Syx_temp = unique_Syx;
    if(num_max_Kx_Sy > 0) {
      ind_max_Kx_Sy = find(max_Kx_Sy != 1);
      cumsum_max_Kx_Sy = cumsum(max_Kx_Sy) - ones(Ky);
      cumsum_max_Kx_Sy = cumsum_max_Kx_Sy(ind_max_Kx_Sy);
      uvec_cumsum_max_Kx_Sy = conv_to<uvec>::from(cumsum_max_Kx_Sy); 
      unique_Syx_temp.rows(uvec_cumsum_max_Kx_Sy);
    }
    
    k_curr = unique_Syx_temp(choose_kr_temp-1,0);
    r_curr = unique_Syx_temp(choose_kr_temp-1,1);
    k_curr_1 = k_curr - 1;
    r_curr_1 = r_curr - 1;
    
    unique_Sy_temp = unique_Sy;
    unique_Sy_temp.shed_row(k_curr_1);
    k_prop = unique_Sy_temp(choose_new_k_temp-1);
    
    ind_unique_Syx = find(unique_Syx.col(0) == k_prop);
    num_unique_Syx = ind_unique_Syx.size();
    r_prop = num_unique_Syx + 1;
    k_prop_1 = k_prop - 1;
    
    // Calculate K_{x,2+}^{*}
    max_Kx_Sy_temp = zeros(Ky+1);
    max_Kx_Sy_temp.head(Ky) = max_Kx_Sy;
    max_Kx_Sy_temp(k_curr_1) = max_Kx_Sy_temp(k_curr_1) - 1;
    max_Kx_Sy_temp(k_prop_1) = max_Kx_Sy_temp(k_prop_1) + 1;
    ind_max_Kx_Sy_temp = find(max_Kx_Sy_temp == 1);
    num_max_Kx_Sy_temp = ind_max_Kx_Sy_temp.size();
    
    Kx2_prop = Kyx - num_max_Kx_Sy_temp;
    
    // Calculate the acceptance probability
    n_k_prop = n_k(k_prop_1);
    
    n_k_curr = n_k(k_curr_1);
    n_rk_curr = n_rk(k_curr_1,r_curr_1);
    
    ind_Syx = find(Sy == k_curr && Sx == r_curr);
    num_Syx = ind_Syx.size();
    
    vec Y_temp = Y(ind_Syx);
    vec M_temp = M(ind_Syx);
    mat matM_temp = matM.rows(ind_Syx);
    mat matX_temp = matX.rows(ind_Syx);
    
    Ybeta_curr = YbetaPars.col(k_curr_1);
    Mbeta_curr = MbetaPars.col(k_curr_1);
    Ysig_curr = sqrt(Ysig2Pars(k_curr_1)) * ones(num_Syx);
    Msig_curr = sqrt(Msig2Pars(k_curr_1)) * ones(num_Syx);
    
    Ybeta_prop = YbetaPars.col(k_prop_1);
    Mbeta_prop = MbetaPars.col(k_prop_1);
    Ysig_prop = sqrt(Ysig2Pars(k_prop_1)) * ones(num_Syx);
    Msig_prop = sqrt(Msig2Pars(k_prop_1)) * ones(num_Syx);
    
    ratio =
      ((lgamma(n_k_curr - n_rk_curr) + lgamma(n_k_prop + n_rk_curr)) - (lgamma(n_k_curr) + lgamma(n_k_prop))) +
      ((lgamma(alpha_omega + n_k_curr) + lgamma(alpha_omega + n_k_prop)) - (lgamma(alpha_omega + n_k_curr - n_rk_curr) + lgamma(alpha_omega + n_k_prop + n_rk_curr))) +
      as_scalar((sum(log_normpdf(Y_temp, matM_temp * Ybeta_prop, Ysig_prop)) + sum(log_normpdf(M_temp, matX_temp * Mbeta_prop, Msig_prop))) -
      (sum(log_normpdf(Y_temp, matM_temp * Ybeta_curr, Ysig_curr)) + sum(log_normpdf(M_temp, matX_temp * Mbeta_curr, Msig_curr)))) +
      log(Kx2_curr/Kx2_prop);
    ratio = exp(ratio);
    
    // If ratio > runif(1), retrun "proposed" and o.w.,do nothing, i.e. retrun "current".
    if (ratio > randu(1)[0]){
      // reorder omega parameters
      uvec ind_unique_Syx_curr = find(unique_Syx.col(0)==k_curr && unique_Syx.col(1)==r_curr);
      uvec ind_Sy0_prop = find(unique_Syx.col(0) <= k_prop);
      int num_Sy0_prop = ind_Sy0.size();
      ZpiPars = reorder(ZpiPars,ind_unique_Syx_curr(0),num_Sy0_prop);
      XpiPars = reorder(XpiPars,ind_unique_Syx_curr(0),num_Sy0_prop);
      XmuPars = reorder(XmuPars,ind_unique_Syx_curr(0),num_Sy0_prop);
      Xtau2Pars = reorder(Xtau2Pars,ind_unique_Syx_curr(0),num_Sy0_prop);
      
      // relabel clusters
      ind_Syx = find(Sy == k_curr && Sx == r_curr);
      num_Syx = ind_Syx.size();
      for (int ii = 0; ii < num_Syx; ii++) {
        Sy(ind_Syx(ii)) = k_prop;
        Sx(ind_Syx(ii)) = r_prop;
      }
      
      // relabel X cluster
      ind_Syx = find(Sy == k_curr && Sx >= r_curr);
      num_Syx = ind_Syx.size();
      for (int ii = 0; ii < num_Syx; ii++) {
        Sx(ind_Syx(ii)) -= 1;
      }
      
      // wrap-up
      Syx = join_rows(Sy, Sx);
      
      // Determine unique clusters from the cluster membership variable, Syx
      unique_Syx = unique_rows(Syx);
      
      // Make vector of y-clusters
      unique_Sy = unique(Sy);
      
      // Calculate the number of clusters
      Kyx = unique_Syx.n_rows;
      
      // Calculate the number of y-clusters
      Ky = unique_Sy.size();
      
      // Find the largest number of x clusters
      max_Sx = Sx.max();
      
      // Calculate the number of subjects in each y-cluster and store in vector n_k.
      // Calculate the number of subjects in each x-cluster and store in matrix n_rk.
      // Use k to index y-clusters and r to index x-clusters.
      n_k.resize(Ky);
      n_rk.resize(Ky,max_Sx);
      max_Kx_Sy.resize(Ky);
      for (int k = 0; k < Ky; k++) {
        ind_Sy = find(Sy == (k+1));
        num_Sy = ind_Sy.size();
        n_k(k) = num_Sy;
        
        ind_Sy0 = find(unique_Syx.col(0) == (k+1));
        num_Sy0 = ind_Sy0.size();
        max_Kx_Sy(k) = num_Sy0;
        for(int r = 0; r < num_Sy0; r++) {
          ind_Syx = find(Sy == (k+1) && Sx == (r+1));
          num_Syx = ind_Syx.size();
          n_rk(k,r) = num_Syx;
        }
      }
    }
  }
  
  // 2nd move or 3rd move
  if(randu(1)[0] < 0.5){
    // ---------------------------------- 2nd move ----------------------------------
    if(Kyx>Ky){
      // Calculate K_{x,1} and K_{x,2+}
      ind_max_Kx_Sy = find(max_Kx_Sy == 1);
      num_max_Kx_Sy = ind_max_Kx_Sy.size();
      
      Kx2_curr = Kyx - num_max_Kx_Sy;
      
      // Choose cluster (k,r) and h
      int choose_kr_temp = rmultinom_cpp(ones(Kx2_curr));
      
      unique_Syx_temp = unique_Syx;
      if(num_max_Kx_Sy > 0) {
        ind_max_Kx_Sy = find(max_Kx_Sy != 1);
        cumsum_max_Kx_Sy = cumsum(max_Kx_Sy) - ones(Ky);
        cumsum_max_Kx_Sy = cumsum_max_Kx_Sy(ind_max_Kx_Sy);
        uvec_cumsum_max_Kx_Sy = conv_to<uvec>::from(cumsum_max_Kx_Sy);
        unique_Syx_temp.rows(uvec_cumsum_max_Kx_Sy);
      }
      
      k_curr = unique_Syx_temp(choose_kr_temp-1,0);
      r_curr = unique_Syx_temp(choose_kr_temp-1,1);
      k_prop = Ky + 1;
      r_prop = 1;
      
      k_curr_1 = k_curr - 1;
      r_curr_1 = r_curr - 1;
      k_prop_1 = k_prop - 1;
      
      // Calculate K_{x,1}^{*}
      max_Kx_Sy_temp = zeros(Ky+1);
      max_Kx_Sy_temp.head(Ky) = max_Kx_Sy;
      max_Kx_Sy_temp(k_curr_1) = max_Kx_Sy_temp(k_curr_1) - 1;
      max_Kx_Sy_temp(k_prop_1) = max_Kx_Sy_temp(k_prop_1) + 1;
      ind_max_Kx_Sy_temp = find(max_Kx_Sy_temp == 1);
      num_max_Kx_Sy_temp = ind_max_Kx_Sy_temp.size();
      
      Kx1_prop = num_max_Kx_Sy_temp;
      
      // Calculate the acceptance probability
      n_k_curr = n_k(k_curr_1);
      n_rk_curr = n_rk(k_curr_1,r_curr_1);
      
      ind_Syx = find(Sy == k_curr && Sx == r_curr);
      num_Syx = ind_Syx.size();
      
      vec Y_temp = Y(ind_Syx);
      vec M_temp = M(ind_Syx);
      mat matM_temp = matM.rows(ind_Syx);
      mat matX_temp = matX.rows(ind_Syx);
      
      Ybeta_curr = YbetaPars.col(k_curr_1);
      Mbeta_curr = MbetaPars.col(k_curr_1);
      Ysig_curr = sqrt(Ysig2Pars(k_curr_1)) * ones(num_Syx);
      Msig_curr = sqrt(Msig2Pars(k_curr_1)) * ones(num_Syx);
      
      ratio =
        (lgamma(n_k_curr - n_rk_curr) + lgamma(n_rk_curr) - lgamma(n_k_curr)) +
        ((lgamma(alpha_omega + n_k_curr) + lgamma(alpha_omega)) - (lgamma(alpha_omega + n_k_curr - n_rk_curr) + lgamma(alpha_omega + n_rk_curr))) +
        as_scalar((sum(log(f0_y(ind_Syx) % f0_m(ind_Syx)))) - (sum(log_normpdf(Y_temp, matM_temp * Ybeta_curr, Ysig_curr)) + sum(log_normpdf(M_temp, matX_temp * Mbeta_curr, Msig_curr)))) +
        log(Kx2_curr/(Kx1_prop*Kyx)) + log(alpha_theta);
      ratio = exp(ratio);
      
      // If ratio > runif(1), retrun "proposed" and o.w.,do nothing, i.e. retrun "current".
      if (ratio > randu(1)[0]){
        // Add theta parameters
        Ysig2_prop = rscainvchisq_cpp(1, a_Ysig2, b_Ysig2)[0];
        Ybeta_prop = mvnrnd(a_Ybeta, B_Ybeta);
        YbetaPars.insert_cols(Ky,Ybeta_prop);
        Ysig2Pars.insert_cols(Ky,1);
        Ysig2Pars(Ky) = Ysig2_prop;
        
        Msig2_prop = rscainvchisq_cpp(1, a_Msig2, b_Msig2)[0];
        Mbeta_prop = mvnrnd(a_Mbeta, B_Mbeta);
        MbetaPars.insert_cols(Ky,Mbeta_prop);
        Msig2Pars.insert_cols(Ky,1);
        Msig2Pars(Ky) = Msig2_prop;
        
        // reorder omega parameters
        uvec ind_unique_Syx_curr = find(unique_Syx.col(0)==k_curr && unique_Syx.col(1)==r_curr);
        ZpiPars = reorder(ZpiPars,ind_unique_Syx_curr(0),Kyx-1);
        XpiPars = reorder(XpiPars,ind_unique_Syx_curr(0),Kyx-1);
        XmuPars = reorder(XmuPars,ind_unique_Syx_curr(0),Kyx-1);
        Xtau2Pars = reorder(Xtau2Pars,ind_unique_Syx_curr(0),Kyx-1);
        
        // relabel clusters
        ind_Syx = find(Sy == k_curr && Sx == r_curr);
        num_Syx = ind_Syx.size();
        for(int ii = 0; ii < num_Syx; ii++) {
          Sy(ind_Syx(ii)) = k_prop;
          Sx(ind_Syx(ii)) = r_prop;
        }
        
        // relabel X cluster
        ind_Syx = find(Sy == k_curr && Sx > r_curr);
        num_Syx = ind_Syx.size();
        for(int ii = 0; ii < num_Syx; ii++) {
          Sx(ind_Syx(ii)) -= 1;
        }
        
        // wrap-up
        Syx = join_rows(Sy, Sx);
        
        // Determine unique clusters from the cluster membership variable, Syx
        unique_Syx = unique_rows(Syx);
        
        // Make vector of y-clusters
        unique_Sy = unique(Sy);
        
        // Calculate the number of clusters
        Kyx = unique_Syx.n_rows;
        
        // Calculate the number of y-clusters
        Ky = unique_Sy.size();
        
        // Find the largest number of x clusters
        max_Sx = Sx.max();
        
        // Calculate the number of subjects in each y-cluster and store in vector n_k.
        // Calculate the number of subjects in each x-cluster and store in matrix n_rk.
        // Use k to index y-clusters and r to index x-clusters.
        n_k.resize(Ky);
        n_rk.resize(Ky,max_Sx);
        max_Kx_Sy.resize(Ky);
        for (int k = 0; k < Ky; k++) {
          ind_Sy = find(Sy == (k+1));
          num_Sy = ind_Sy.size();
          n_k(k) = num_Sy;
          
          ind_Sy0 = find(unique_Syx.col(0) == (k+1));
          num_Sy0 = ind_Sy0.size();
          max_Kx_Sy(k) = num_Sy0;
          for(int r = 0; r < num_Sy0; r++) {
            ind_Syx = find(Sy == (k+1) && Sx == (r+1));
            num_Syx = ind_Syx.size();
            n_rk(k,r) = num_Syx;
          }
        }
      }
    }
    
  } else {
    // ---------------------------------- 3rd move ----------------------------------
    // Calculate K_{x,1} and K_{x,2+}
    ind_max_Kx_Sy = find(max_Kx_Sy == 1);
    num_max_Kx_Sy = ind_max_Kx_Sy.size();
    
    if(num_max_Kx_Sy>0 && Ky>1){
      
      Kx1_curr = num_max_Kx_Sy;
      
      // Choose cluster (k,r) and h
      int choose_kr_temp = rmultinom_cpp(ones(Kx1_curr));
      int choose_new_k_temp = rmultinom_cpp(ones(Ky-1));
      
      cumsum_max_Kx_Sy = cumsum(max_Kx_Sy) - ones(Ky);
      cumsum_max_Kx_Sy = cumsum_max_Kx_Sy(ind_max_Kx_Sy);
      uvec_cumsum_max_Kx_Sy = conv_to<uvec>::from(cumsum_max_Kx_Sy);
      unique_Syx_temp = unique_Syx.rows(uvec_cumsum_max_Kx_Sy);
      
      k_curr = unique_Syx_temp(choose_kr_temp-1,0);
      r_curr = 1;
      k_curr_1 = k_curr - 1;
      r_curr_1 = r_curr - 1;
      
      unique_Sy_temp = unique_Sy;
      unique_Sy_temp.shed_row(k_curr_1);
      k_prop = unique_Sy_temp(choose_new_k_temp-1);
      
      ind_unique_Syx = find(unique_Syx.col(0) == k_prop);
      num_unique_Syx = ind_unique_Syx.size();
      r_prop = num_unique_Syx + 1;
      k_prop_1 = k_prop - 1;
      
      // Calculate K_{x,2+}^{*}
      max_Kx_Sy_temp = zeros(Ky+1);
      max_Kx_Sy_temp.head(Ky) = max_Kx_Sy;
      max_Kx_Sy_temp(k_curr_1) = max_Kx_Sy_temp(k_curr_1) - 1;
      max_Kx_Sy_temp(k_prop_1) = max_Kx_Sy_temp(k_prop_1) + 1;
      ind_max_Kx_Sy_temp = find(max_Kx_Sy_temp == 1);
      num_max_Kx_Sy_temp = ind_max_Kx_Sy_temp.size();
      
      Kx2_prop = Kyx - num_max_Kx_Sy_temp;
      
      // Calculate the acceptance probability
      n_k_prop = n_k(k_prop_1);
      
      n_k_curr = n_k(k_curr_1);
      n_rk_curr = n_rk(k_curr_1,r_curr_1);
      
      ind_Syx = find(Sy == k_curr && Sx == r_curr);
      num_Syx = ind_Syx.size();
      
      vec Y_temp = Y(ind_Syx);
      vec M_temp = M(ind_Syx);
      mat matM_temp = matM.rows(ind_Syx);
      mat matX_temp = matX.rows(ind_Syx);
      
      Ybeta_prop = YbetaPars.col(k_prop_1);
      Mbeta_prop = MbetaPars.col(k_prop_1);
      Ysig_prop = sqrt(Ysig2Pars(k_prop_1)) * ones(num_Syx);
      Msig_prop = sqrt(Msig2Pars(k_prop_1)) * ones(num_Syx);
      
      ratio =
        (lgamma(n_k_prop + n_rk_curr) - (lgamma(n_rk_curr) + lgamma(n_k_prop))) +
        ((lgamma(alpha_omega + n_rk_curr) + lgamma(alpha_omega + n_k_prop)) - (lgamma(alpha_omega) + lgamma(alpha_omega + n_k_prop + n_rk_curr))) +
        as_scalar((sum(log_normpdf(Y_temp, matM_temp * Ybeta_prop, Ysig_prop)) + sum(log_normpdf(M_temp, matX_temp * Mbeta_prop, Msig_prop))) - (sum(log(f0_y(ind_Syx) % f0_m(ind_Syx))))) -
        log((Kx1_curr * Kyx - 1)/Kx2_prop) - log(alpha_theta);
      ratio = exp(ratio);
      
      // If ratio > runif(1), retrun "proposed" and o.w.,do nothing, i.e. retrun "current".
      if (ratio > randu(1)[0]){
        // Delete theta parameters
        YbetaPars.shed_col(k_curr_1);
        Ysig2Pars.shed_col(k_curr_1);
        MbetaPars.shed_col(k_curr_1);
        Msig2Pars.shed_col(k_curr_1);
        
        // reorder omega parameters
        uvec ind_unique_Syx_curr = find(unique_Syx.col(0)==k_curr && unique_Syx.col(1)==r_curr);
        uvec ind_Sy0_prop = find(unique_Syx.col(0) <= k_prop);
        int num_Sy0_prop = ind_Sy0.size();
        ZpiPars = reorder(ZpiPars,ind_unique_Syx_curr(0),num_Sy0_prop);
        XpiPars = reorder(XpiPars,ind_unique_Syx_curr(0),num_Sy0_prop);
        XmuPars = reorder(XmuPars,ind_unique_Syx_curr(0),num_Sy0_prop);
        Xtau2Pars = reorder(Xtau2Pars,ind_unique_Syx_curr(0),num_Sy0_prop);
        
        // relabel clusters
        ind_Syx = find(Sy == k_curr && Sx == r_curr);
        num_Syx = ind_Syx.size();
        for(int ii = 0; ii < num_Syx; ii++) {
          Sy(ind_Syx(ii)) = k_prop;
          Sx(ind_Syx(ii)) = r_prop;
        }
        
        // relabel Y cluster
        ind_Sy = find(Sy >= k_curr);
        num_Sy = ind_Sy.size();
        for(int ii = 0; ii < num_Sy; ii++) {
          Sy(ind_Sy(ii)) -= 1;
        }
        
        // wrap-up
        Syx = join_rows(Sy, Sx);
        
        // Determine unique clusters from the cluster membership variable, Syx
        unique_Syx = unique_rows(Syx);
        
        // Make vector of y-clusters
        unique_Sy = unique(Sy);
        
        // Calculate the number of clusters
        Kyx = unique_Syx.n_rows;
        
        // Calculate the number of y-clusters
        Ky = unique_Sy.size();
        
        // Find the largest number of x clusters
        max_Sx = Sx.max();
        
        // Calculate the number of subjects in each y-cluster and store in vector n_k.
        // Calculate the number of subjects in each x-cluster and store in matrix n_rk.
        // Use k to index y-clusters and r to index x-clusters.
        n_k.resize(Ky);
        n_rk.resize(Ky,max_Sx);
        max_Kx_Sy.resize(Ky);
        for (int k = 0; k < Ky; k++) {
          ind_Sy = find(Sy == (k+1));
          num_Sy = ind_Sy.size();
          n_k(k) = num_Sy;
          
          ind_Sy0 = find(unique_Syx.col(0) == (k+1));
          num_Sy0 = ind_Sy0.size();
          max_Kx_Sy(k) = num_Sy0;
          for(int r = 0; r < num_Sy0; r++) {
            ind_Syx = find(Sy == (k+1) && Sx == (r+1));
            num_Syx = ind_Syx.size();
            n_rk(k,r) = num_Syx;
          }
        }
      }
    }
  }
  
  return List::create(_["Sy"]  = Sy,
                      _["Sx"]  = Sx,
                      _["Syx"] = Syx,
                      _["unique_Sy"]  = unique_Sy,
                      _["unique_Syx"] = unique_Syx,
                      _["YbetaPars"] = YbetaPars,
                      _["Ysig2Pars"] = Ysig2Pars,
                      _["MbetaPars"] = MbetaPars,
                      _["Msig2Pars"] = Msig2Pars,
                      _["ZpiPars"]   = ZpiPars,
                      _["XpiPars"]   = XpiPars,
                      _["XmuPars"]   = XmuPars,
                      _["Xtau2Pars"] = Xtau2Pars,
                      _["Ky"]   = Ky, 
                      _["Kyx"]  = Kyx, 
                      _["n_k"]  = n_k, 
                      _["n_rk"] = n_rk,
                      _["max_Sx"]    = max_Sx,
                      _["max_Kx_Sy"] = max_Kx_Sy);
}

// [[Rcpp::export]]
List EDPpostmc_concon(int const& p_C1, int const& p_C2,
                      vec const& n_k, mat const& n_rk, 
                      int const& Ky, vec const& max_Kx_Sy, 
                      int const& D, int const& num_MC, 
                      mat const& YbetaPars, mat const& Ysig2Pars,
                      mat const& MbetaPars, mat const& Msig2Pars,
                      cube const& ZpiPars_array, cube const& XpiPars_array, 
                      cube const& XmuPars_array, cube const& Xtau2Pars_array,
                      double const& alpha_theta, double const& alpha_omega,
                      vec const& f0_z0_mc, vec const& f0_z1_mc, 
                      mat const& Ybeta_prior, 
                      mat const& Mbeta_prior, vec const& Msig2_prior,
                      double const& a_Ysig2, double const& b_Ysig2,
                      vec const& a_Ybeta, mat const& B_Ybeta,
                      double const& a_Msig2, double const& b_Msig2,
                      vec const& a_Mbeta, mat const& B_Mbeta,
                      double const& a_pi, double const& b_pi,
                      double const& a_mu, double const& b_mu,
                      double const& a_tau2, double const& b_tau2){
  // initialize some dummy vars to store things
  double c_mc;
  double Zpi_prop, Xpi_prop, Xmu_prop, Xtau_prop, Xtau2_prop;
  double Ysig2_prop, Msig2_prop;
  vec Ybeta_prop, Mbeta_prop;
  
  vec zeros_D = zeros(D);  
  vec ones_D = ones(D);
  vec z0_D = zeros_D; 
  vec z1_D = ones_D;
  
  uvec uones_D = conv_to<uvec>::from(ones_D);
  
  // Initialize matrices and vectors
  int p_Z = 1;
  
  mat Ysig2_mc = Ysig2Pars;
  mat Ybeta_mc = YbetaPars;
  
  mat Msig2_mc = Msig2Pars;
  mat Mbeta_mc = MbetaPars;
  
  cube Zpi_mc = ZpiPars_array;
  cube Xpi_mc = XpiPars_array;
  cube Xmu_mc = XmuPars_array;
  cube Xtau2_mc = Xtau2Pars_array;
  
  vec M0_mc(D);
  vec M1_mc(D);
  
  mat C_mc(D, (p_C1+p_C2));
  
  uvec Sy_mc = uones_D;
  uvec Sx_mc = uones_D;
  
  // Assign D observations one of the current y-clusters or to a new
  // y-cluster using draws from a multinomial distribution with probabilities 'prob_ycluster'
  vec prob_ycluster(Ky + 1);
  for(int k = 0; k < Ky; k++) {
    prob_ycluster(k) = n_k(k);
  }
  prob_ycluster(Ky) = alpha_theta;
  
  for(int d = 0; d < D; d++) {
    Sy_mc(d) = rmultinom_cpp(prob_ycluster);
  }
  
  // Store indices of non-empty y-clusters 
  // (some y-clusters could be empty and there could be a new cluster)
  uvec unique_Sy_mc = unique(Sy_mc);
  
  // Store number of y-clusters including empty ones and possibly a new one
  int Ky_mc = unique_Sy_mc.max();
  
  // Store number of x-clusters in each y-cluster in a vector
  vec max_Kx_Sy_mc = zeros(Ky_mc);
  
  // If a new y-cluster was opened, draw parameters for this cluster from priors
  if (Ky_mc > Ky) {
    max_Kx_Sy_mc(Ky) = 1;
    
    for (int q = 0; q < p_Z; q++) {
      Zpi_prop = R::rbeta(a_pi,b_pi);
      Zpi_mc(q, 0, Ky) = Zpi_prop;
    }
    
    for(int q = 0; q < p_C1; q++) {
      Xpi_prop = R::rbeta(a_pi,b_pi);
      Xpi_mc(q, 0, Ky) = Xpi_prop;
    }
    
    for (int q = 0; q < p_C2; q++) {
      Xtau2_prop = rscainvchisq_cpp(1, a_tau2, b_tau2)[0];
      Xtau2_mc(q, 0, Ky) = Xtau2_prop;
      
      Xmu_prop = sqrt(Xtau2_prop/b_mu) * randn(1)[0] + a_mu;
      Xmu_mc(q, 0, Ky) = Xmu_prop;
    }
    
    Msig2_prop = rscainvchisq_cpp(1, a_Msig2, b_Msig2)[0];
    Mbeta_prop = mvnrnd(a_Mbeta, B_Mbeta);
    Mbeta_mc.insert_cols(Ky, Mbeta_prop);
    Msig2_mc.insert_cols(Ky, 1);
    Msig2_mc(Ky) = Msig2_prop;
    
    Ysig2_prop = rscainvchisq_cpp(1, a_Ysig2, b_Ysig2)[0];
    Ybeta_prop = mvnrnd(a_Ybeta, B_Ybeta);
    Ybeta_mc.insert_cols(Ky, Ybeta_prop);
    Ysig2_mc.insert_cols(Ky, 1);
    Ysig2_mc(Ky) = Ysig2_prop;
  }
  
  // Draw x-clusters within each y-cluster
  vec n_k_mc(Ky_mc);
  for (int k = 0; k < Ky_mc; k++) {
    // Number of observations (out of D total) in each y-cluster
    uvec ind_Sy_mc = find(Sy_mc == (k+1));
    int num_Sy_mc = ind_Sy_mc.n_elem;
    n_k_mc(k) = num_Sy_mc;
    
    if (k < Ky && num_Sy_mc > 0) {
      int max_Kx_Sy_temp = max_Kx_Sy(k);
      vec prob_xcluster(max_Kx_Sy_temp + 1);
      for(int r = 0; r < max_Kx_Sy_temp; r++) {
        prob_xcluster(r) = n_rk(k,r);
      }
      prob_xcluster(max_Kx_Sy_temp) = alpha_omega;
      
      uvec Sx_mc_temp(num_Sy_mc);
      for(int d = 0; d < num_Sy_mc; d++) {
        Sx_mc_temp(d) = rmultinom_cpp(prob_xcluster);
      }
      Sx_mc.rows(ind_Sy_mc) = Sx_mc_temp;
      
      // If a new x-cluster was opened, draw parameter from prior
      int max_Sx_mc_temp = Sx_mc_temp.max();
      max_Kx_Sy_mc(k) = max_Sx_mc_temp;
      if (max_Sx_mc_temp > max_Kx_Sy_temp) {
        for (int q = 0; q < p_Z; q++) {
          Zpi_prop = R::rbeta(a_pi,b_pi);
          Zpi_mc(q, max_Kx_Sy_temp, k) = Zpi_prop;
        }
        
        for (int q = 0; q < p_C1; q++) {
          Xpi_prop = R::rbeta(a_pi,b_pi);
          Xpi_mc(q, max_Kx_Sy_temp, k) = Xpi_prop;
        }
        
        for (int q = 0; q < p_C2; q++) {
          Xtau2_prop = rscainvchisq_cpp(1, a_tau2, b_tau2)[0];
          Xtau2_mc(q, max_Kx_Sy_temp, k) = Xtau2_prop;
          
          Xmu_prop = sqrt(Xtau2_prop/b_mu) * randn(1)[0] + a_mu;
          Xmu_mc(q, max_Kx_Sy_temp, k) = Xmu_prop;
        }
      }
    }
  }
  
  // Draw confounders for each of the D observations-------------------
  int Kyx_mc = sum(max_Kx_Sy_mc);
  mat n_rk_mc(Ky_mc, Kyx_mc);
  for (int k = 0; k < Ky_mc; k++) {
    int max_Kx_Sy_temp = max_Kx_Sy_mc(k);
    
    for (int r = 0; r < max_Kx_Sy_temp; r++) {
      uvec ind_Syx_mc = find(Sy_mc == (k+1) && Sx_mc == (r+1));
      int num_Syx_mc = ind_Syx_mc.n_elem;
      n_rk_mc(k,r) = num_Syx_mc;
      
      for (int q = 0; q < p_C1; q++) {
        Xpi_prop = Xpi_mc(q, r, k);
        for (int d = 0; d < num_Syx_mc; d++) {
          C_mc(ind_Syx_mc(d), q) = R::rbinom(1, Xpi_prop);
        }
      }
      for (int q = 0; q < p_C2; q++) {
        Xmu_prop = Xmu_mc(q, r, k);
        Xtau_prop = sqrt(Xtau2_mc(q, r, k));
        for (int d = 0; d < num_Syx_mc; d++) {
          C_mc(ind_Syx_mc(d),(p_C1+q)) = Xtau_prop * randn(1)[0] + Xmu_prop;
        }
      }
    }
  }
  // End draw of confounders-----------------------------------
  
  // Calculate f0_x_z1/z0, f0_m_z1/z0 and f0_y_z1/z0 for use in calculating causal effect.
  // Average covariate distribution over prior for each x_i.
  double ab_tau = a_tau2*b_tau2;
  double a_tau_new = (a_tau2+1)/2;
  double a_tau_half = a_tau2/2;
  double b_mu_ratio = b_mu/(b_mu + 1);
  double margin_part1_mc = (tgamma(a_tau_new)/tgamma(a_tau_half)) * 
    sqrt(b_mu_ratio/M_PI) * pow(ab_tau,a_tau_half);
  
  mat f0_c_all_mc(D, (p_C1+p_C2));
  for (int d = 0; d < D; d++) {
    // Binary confounders (not including treatment)
    // Beta-Binomial
    for (int q = 0; q < p_C1; q++) {
      c_mc = C_mc(d,q);
      Xpi_prop = R::beta(a_pi + c_mc, b_pi + 1 - c_mc); // R::beta(a_pi, b_pi) = 1
      f0_c_all_mc(d, q) = Xpi_prop;
    }
    // Continuous confounders
    for (int q = 0; q < p_C2; q++) {
      c_mc = C_mc(d, (p_C1+q));
      double margin_part2_mc = pow((ab_tau + b_mu_ratio * pow((c_mc-a_mu),2)),-a_tau_new);
      double margin_mc = margin_part1_mc * margin_part2_mc;
      f0_c_all_mc(d, (p_C1+q)) = margin_mc;
    }
  }
  // loop
  
  // Take product (confounders are assumed to be locally independent)
  // f0_z0_mc and f0_z1_mc were calculated before Gibbs sampling
  // Result is vector of size D
  vec f0_c_mc = prod(f0_c_all_mc, 1);
  vec f0_x0_mc = f0_c_mc % f0_z0_mc;
  vec f0_x1_mc = f0_c_mc % f0_z1_mc;
  
  // Set treatment to 0 and 1, and add column of 1'Syx.
  mat X0_temp = join_rows(z0_D, C_mc); // set treatment to 0
  mat X1_temp = join_rows(z1_D, C_mc); // set treatment to 1
  
  mat matX0_mc = join_rows(ones_D, X0_temp);  
  mat matX1_mc = join_rows(ones_D, X1_temp);
  
  // Calculate P(m|z,u;eta_k^*,omega_r|k^*) for mediation
  // Draw mediation for each of the D observations--------------------
  // Initialize
  mat lambda_M_x0_mc = zeros(D, (Ky_mc+1));
  mat lambda_M_x1_mc = zeros(D, (Ky_mc+1));
  
  // (Use unique_Sy_mc rather than 1:Ky_mc in case some y-clusters are empty)
  for (int k = 0; k < Ky_mc; k++) {
    int n_k_temp = n_k_mc(k);
    int max_Kx_Sy_temp = max_Kx_Sy_mc(k);
    
    vec sum_f_c_mc = zeros_D;
    vec sum_f_x0_mc = zeros_D;
    vec sum_f_x1_mc = zeros_D;
    
    for (int r = 0; r < max_Kx_Sy_temp; r++) {
      int n_rk_temp = n_rk_mc(k,r);
      
      vec prob_c_mc = ones_D;
      vec prob_z0_mc = ones_D;
      vec prob_z1_mc = ones_D;
      // Calculate f(z,u;omega*_r|k) for binary confounders
      for (int q = 0; q < p_C1; q++) {
        for (int d = 0; d < D; d++) {
          double dbinom_temp = R::dbinom(C_mc(d,q), 1, Xpi_mc(q, r, k), false);
          prob_c_mc(d) = prob_c_mc(d) * dbinom_temp;
        }
      }
      // Calculate f(z,u;omega*_r|k) for continuous confounders
      for (int q = 0; q < p_C2; q++) {
        vec dnorm_temp = normpdf(C_mc.col(p_C1+q), Xmu_mc(q, r, k), sqrt(Xtau2_mc(q, r, k)));
        prob_c_mc = prob_c_mc % dnorm_temp;
      }
      // Calculate f(z,u;omega*_r|k) for binary treatment
      for (int q = 0; q < p_Z; q++) {
        double dbinom_temp = Zpi_mc(q, r, k);
        prob_z0_mc = prob_z0_mc * (1 - dbinom_temp);
        prob_z1_mc = prob_z1_mc * dbinom_temp;
      }
      
      // Calculate sum_{r=1}^{K_r} \frac{n_{rk}}{alpha_omega+n_k}*f(z,u;omega*_r|k)
      // (This is part of the denominator of E[ Y | M=m, Z=z, U=u, theta^star, omega^star, Syx])
      double prob_omega_rk_mc = n_rk_temp/(alpha_omega + n_k_temp);
      sum_f_c_mc = sum_f_c_mc + prob_omega_rk_mc * prob_c_mc;
      sum_f_x0_mc = sum_f_x0_mc + prob_omega_rk_mc * prob_z0_mc % prob_c_mc;
      sum_f_x1_mc = sum_f_x1_mc + prob_omega_rk_mc * prob_z1_mc % prob_c_mc;
    }
    // Calculate sum_(k=1)^{Ky} lambda_k^m (z,u) f( m | z, r, theta_j*)
    // (This is part of the numerator of P(M | Z=z, U=u, theta^star, omega^star, Syx)
    // double prob_theta_k_mc = n_k_temp/(alpha_theta + D);
    double prob_theta_k_mc = n_k_temp;
    double prob_omega_k_mc = alpha_omega/(alpha_omega + n_k_temp);
    vec lambda_M_x0_k_mc = prob_theta_k_mc * (prob_omega_k_mc * f0_x0_mc + sum_f_x0_mc);
    vec lambda_M_x1_k_mc = prob_theta_k_mc * (prob_omega_k_mc * f0_x1_mc + sum_f_x1_mc);
    
    lambda_M_x0_mc.col(k) = lambda_M_x0_k_mc;
    lambda_M_x1_mc.col(k) = lambda_M_x1_k_mc;
  }
  double prob_theta_K_mc = alpha_theta;
  lambda_M_x0_mc.col(Ky_mc) = prob_theta_K_mc * f0_x0_mc;
  lambda_M_x1_mc.col(Ky_mc) = prob_theta_K_mc * f0_x1_mc;
  
  // Check 0 denominator
  vec sum_lambda_M_x0_mc = sum(lambda_M_x0_mc,1);
  vec sum_lambda_M_x1_mc = sum(lambda_M_x1_mc,1);
  uvec ind_M_x0_mc = find(sum_lambda_M_x0_mc == 0);
  uvec ind_M_x1_mc = find(sum_lambda_M_x1_mc == 0);
  sum_lambda_M_x0_mc(ind_M_x0_mc) = ones(ind_M_x0_mc.n_elem);
  sum_lambda_M_x1_mc(ind_M_x1_mc) = ones(ind_M_x1_mc.n_elem);
  lambda_M_x0_mc.each_col() /= sum_lambda_M_x0_mc;
  lambda_M_x1_mc.each_col() /= sum_lambda_M_x1_mc;
  vec lambda_M_x0_K_mc = lambda_M_x0_mc.col(Ky_mc);
  vec lambda_M_x1_K_mc = lambda_M_x1_mc.col(Ky_mc);
  
  uvec Sm0_mc = uones_D;
  uvec Sm1_mc = uones_D;
  for(int d = 0; d < D; d++) {
    vec prob_m0_temp = lambda_M_x0_mc.row(d).t();
    Sm0_mc(d) = rmultinom_cpp(prob_m0_temp);
    
    vec prob_m1_temp = lambda_M_x1_mc.row(d).t();
    Sm1_mc(d) = rmultinom_cpp(prob_m1_temp);
  }
  uvec unique_Sm0_mc = unique(Sm0_mc);
  uvec unique_Sm1_mc = unique(Sm1_mc);
  
  // Store number of m-clusters including empty ones and possibly a new one
  int Km_mc = join_cols(unique_Sm0_mc,unique_Sm1_mc).max();
  
  // If a new m-cluster was opened, draw parameters for this cluster from priors
  if (Km_mc > Ky_mc) {
    Msig2_prop = rscainvchisq_cpp(1, a_Msig2, b_Msig2)[0];
    Mbeta_prop = mvnrnd(a_Mbeta, B_Mbeta);
    Mbeta_mc.insert_cols(Ky_mc, Mbeta_prop);
    Msig2_mc.insert_cols(Ky_mc, 1);
    Msig2_mc(Ky_mc) = Msig2_prop;
    
    Ysig2_prop = rscainvchisq_cpp(1, a_Ysig2, b_Ysig2)[0];
    Ybeta_prop = mvnrnd(a_Ybeta, B_Ybeta);
    Ybeta_mc.insert_cols(Ky_mc, Ybeta_prop);
    Ysig2_mc.insert_cols(Ky_mc, 1);
    Ysig2_mc(Ky_mc) = Ysig2_prop;
  }
  
  for (int k = 0; k < Km_mc; k++) {
    vec Mbeta_temp = Mbeta_mc.col(k);
    double Msig_temp = sqrt(Msig2_mc(k));
    
    uvec ind_Sm0_mc = find(Sm0_mc == (k+1));
    int num_Sm0_mc = ind_Sm0_mc.n_elem;
    M0_mc(ind_Sm0_mc) = randn(num_Sm0_mc) * Msig_temp + (matX0_mc.rows(ind_Sm0_mc) * Mbeta_temp);
    
    uvec ind_Sm1_mc = find(Sm1_mc == (k+1));
    int num_Sm1_mc = ind_Sm1_mc.n_elem;
    M1_mc(ind_Sm1_mc) = randn(num_Sm1_mc) * Msig_temp + (matX1_mc.rows(ind_Sm1_mc) * Mbeta_temp);
  }
  
  mat M0x0_temp = join_rows(z0_D, z0_D % M0_mc,  M0_mc, C_mc);
  mat M0x1_temp = join_rows(z1_D, z1_D % M0_mc, M0_mc, C_mc);
  mat M1x1_temp = join_rows(z1_D, z1_D % M1_mc, M1_mc, C_mc);
  mat matM_x0m0_mc = join_rows(ones_D, M0x0_temp);
  mat matM_x1m0_mc = join_rows(ones_D, M0x1_temp);
  mat matM_x1m1_mc = join_rows(ones_D, M1x1_temp);
  
  // End draw of mediation------------------------------------
  
  // Compute f(M | Z=z, U=u, theta, omega, Syx) using D new draws
  // Compute E[Y|M=m,Z=z,U=u,theta,omega,Syx] using D new draws
  // Find f0( M | Z, U, theta) integrating over G_0^theta.
  // Mbeta_prior and Msig2_prior was calculated before Gibbs sampling loop
  // Ybeta_prior and sig2_prior_mc was calculated before Gibbs sampling loop
  // Find E0[Y|M,Z,U,theta] integrating over G_0^theta.
  // Initialize
  vec f0_M0_x0_mc = zeros_D;
  vec f0_M0_x1_mc = zeros_D;
  vec f0_M1_x1_mc = zeros_D;
  vec E0_Y_x0m0_mc = zeros_D;
  vec E0_Y_x1m0_mc = zeros_D;
  vec E0_Y_x1m1_mc = zeros_D;
  
  for(int num = 0; num < num_MC; num++) {
    vec Mbeta_temp = Mbeta_prior.col(num);
    vec Msig_temp = sqrt(Msig2_prior(num)) * ones_D;
    vec meanM_x0_temp = matX0_mc * Mbeta_temp;
    vec meanM_x1_temp = matX1_mc * Mbeta_temp;

    f0_M0_x0_mc = f0_M0_x0_mc + normpdf(M0_mc, meanM_x0_temp, Msig_temp);
    f0_M0_x1_mc = f0_M0_x1_mc + normpdf(M0_mc, meanM_x1_temp, Msig_temp);
    f0_M1_x1_mc = f0_M1_x1_mc + normpdf(M1_mc, meanM_x1_temp, Msig_temp);

    vec Ybeta_temp = Ybeta_prior.col(num);
    E0_Y_x0m0_mc = E0_Y_x0m0_mc + (matM_x0m0_mc * Ybeta_temp);
    E0_Y_x1m0_mc = E0_Y_x1m0_mc + (matM_x1m0_mc * Ybeta_temp);
    E0_Y_x1m1_mc = E0_Y_x1m1_mc + (matM_x1m1_mc * Ybeta_temp);
  }
  f0_M0_x0_mc = f0_M0_x0_mc/num_MC;
  f0_M0_x1_mc = f0_M0_x1_mc/num_MC;
  f0_M1_x1_mc = f0_M1_x1_mc/num_MC;

  E0_Y_x0m0_mc = E0_Y_x0m0_mc/num_MC;
  E0_Y_x1m0_mc = E0_Y_x1m0_mc/num_MC;
  E0_Y_x1m1_mc = E0_Y_x1m1_mc/num_MC;
  
  // Calculate E[Y^{z,M_z'}] using MC integration to integrate
  // E[Y | M=m, Z=z, U=u, theta^star, omega^star, Syx] over covariate distribution
  // Initialize
  vec numer_Y_x0m0_mc = zeros_D;
  vec denom_Y_x0m0_mc = zeros_D;
  vec numer_Y_x1m0_mc = zeros_D;
  vec denom_Y_x1m0_mc = zeros_D;
  vec numer_Y_x1m1_mc = zeros_D;
  vec denom_Y_x1m1_mc = zeros_D;
  
  for (int k = 0; k < Km_mc; k++) {
    vec Mbeta_temp = Mbeta_mc.col(k);
    vec Msig_temp = sqrt(Msig2_mc(k)) * ones_D;
    vec meanM_x0_temp = matX0_mc * Mbeta_temp;
    vec meanM_x1_temp = matX1_mc * Mbeta_temp;
    
    vec Ybeta_temp = Ybeta_mc.col(k);
    vec meanY_x0m0_temp = matM_x0m0_mc * Ybeta_temp;
    vec meanY_x1m0_temp = matM_x1m0_mc * Ybeta_temp;
    vec meanY_x1m1_temp = matM_x1m1_mc * Ybeta_temp;
    
    // Calculate lambda_k^y (m,z,u)
    vec lambda_Y_x0m0_k_mc = lambda_M_x0_mc.col(k) % normpdf(M0_mc, meanM_x0_temp, Msig_temp);
    vec lambda_Y_x1m0_k_mc = lambda_M_x1_mc.col(k) % normpdf(M0_mc, meanM_x1_temp, Msig_temp);
    vec lambda_Y_x1m1_k_mc = lambda_M_x1_mc.col(k) % normpdf(M1_mc, meanM_x1_temp, Msig_temp);
    
    // Calculate sum_(k=1)^{Ky} lambda_k^y (m,z,u)
    // (This is part of the numerator of E[Y | M=m, Z=z, U=u, theta^star, omega^star, Syx])
    numer_Y_x0m0_mc += lambda_Y_x0m0_k_mc % meanY_x0m0_temp;
    numer_Y_x1m0_mc += lambda_Y_x1m0_k_mc % meanY_x1m0_temp;
    numer_Y_x1m1_mc += lambda_Y_x1m1_k_mc % meanY_x1m1_temp;
    
    // Calculate sum_(k=1)^{Ky} lambda_k^y (m,z,u)
    // (This is part of the denominator of E[Y | M=m, Z=z, U=u, theta^star, omega^star, Syx])
    denom_Y_x0m0_mc += lambda_Y_x0m0_k_mc;
    denom_Y_x1m0_mc += lambda_Y_x1m0_k_mc;
    denom_Y_x1m1_mc += lambda_Y_x1m1_k_mc;
  }
  // Calculate lambda_{Ky+1}^{y} (m,z,u)
  vec lambda_Y_x0m0_K_mc = lambda_M_x0_K_mc % f0_M0_x0_mc;
  vec lambda_Y_x1m0_K_mc = lambda_M_x1_K_mc % f0_M0_x1_mc;
  vec lambda_Y_x1m1_K_mc = lambda_M_x1_K_mc % f0_M1_x1_mc;
  
  // Numerator of E[Y | M=m, Z=z, U=u, theta^star, omega^star, Syx]
  numer_Y_x0m0_mc += lambda_Y_x0m0_K_mc % E0_Y_x0m0_mc;
  numer_Y_x1m0_mc += lambda_Y_x1m0_K_mc % E0_Y_x1m0_mc;
  numer_Y_x1m1_mc += lambda_Y_x1m1_K_mc % E0_Y_x1m1_mc;
  
  // Denominator of E[Y | M=m, Z=z, U=u, theta^star, omega^star, Syx]
  denom_Y_x0m0_mc += lambda_Y_x0m0_K_mc;
  denom_Y_x1m0_mc += lambda_Y_x1m0_K_mc;
  denom_Y_x1m1_mc += lambda_Y_x1m1_K_mc;
  
  // ---------------------------------------------------------------------------------------
  // ---------------------------------------------------------------------------------------
  // ---------------------------------------------------------------------------------------
  // Check 0 denominator
  uvec ind_Y_x0m0_mc = find(denom_Y_x0m0_mc == 0);
  uvec ind_Y_x1m0_mc = find(denom_Y_x1m0_mc == 0);
  uvec ind_Y_x1m1_mc = find(denom_Y_x1m1_mc == 0);
  denom_Y_x0m0_mc(ind_Y_x0m0_mc) = ones(ind_Y_x0m0_mc.n_elem);
  denom_Y_x1m0_mc(ind_Y_x1m0_mc) = ones(ind_Y_x1m0_mc.n_elem);
  denom_Y_x1m1_mc(ind_Y_x1m1_mc) = ones(ind_Y_x1m1_mc.n_elem);
  
  // Calculate E[Y_{(d)}^{0,M0}]
  vec E_Y_x0m0_mc = numer_Y_x0m0_mc/denom_Y_x0m0_mc;
  // Calculate E[Y_{(d)}^{1,M0}]
  vec E_Y_x1m0_mc = numer_Y_x1m0_mc/denom_Y_x1m0_mc;
  // Calculate E[Y_{(d)}^{1,M1}]
  vec E_Y_x1m1_mc = numer_Y_x1m1_mc/denom_Y_x1m1_mc;
  
  // Check NaN
  uvec ind_E_Y_x0m0_mc = find_nan(E_Y_x0m0_mc);
  uvec ind_E_Y_x1m0_mc = find_nan(E_Y_x1m0_mc);
  uvec ind_E_Y_x1m1_mc = find_nan(E_Y_x1m1_mc);
  E_Y_x0m0_mc.shed_rows(ind_E_Y_x0m0_mc);
  E_Y_x1m0_mc.shed_rows(ind_E_Y_x1m0_mc);
  E_Y_x1m1_mc.shed_rows(ind_E_Y_x1m1_mc);
  
  // Calculate E[Y_{0,M0}]
  double E_Y_z0m0 = mean(E_Y_x0m0_mc);
  // Calculate E[Y_{1,Mr}]
  double E_Y_z1m0 = mean(E_Y_x1m0_mc);
  // Calculate E[Y_{1,M1}]
  double E_Y_z1m1 = mean(E_Y_x1m1_mc);
  
  // Calculate E[Y_{1,M1}] - E[Y_{1,Mr}]
  double NIE = E_Y_z1m1 - E_Y_z1m0;
  // Calculate E[Y_{1,Mr}] - E[Y_{0,M0}]
  double NDE = E_Y_z1m0 - E_Y_z0m0;
  // Calculate E[Y_{1,M1}] - E[Y_{0,M0}]
  double ATE = E_Y_z1m1 - E_Y_z0m0;
  
  return List::create(_["NIE"] = NIE,
                      _["NDE"] = NDE,
                      _["ATE"] = ATE);
  
  // return List::create(_["NIE"] = E_Y_x0m0,
  //                     _["NDE"] = E_Y_z1mr,
  //                     _["ATE"] = E_Y_x1m1);
}