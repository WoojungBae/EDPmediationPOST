# Define functions --------------------------------------------------------------
invlogit = function(x_temp) { 
  return((1/(1+exp(-x_temp))))
}

# extract from package "lcmix"
rmvgamma = function(n, shape, rate, rho){
  p = 2
  corrmat = matrix(c(1, rho, rho, 1), ncol=p)
  
  # generate standard multivariate normal matrix, convert to CDF
  Z = t(rmvn_cpp(N, MU = rep(0,p), SIG = corrmat))
  cdf = pnorm(Z)
  
  # convert to gamma, return
  rgam = sapply(1:n, function(l) qgamma(cdf[l,], shape[l,], rate))
  
  return(rgam)
}

# extract from package "sn"
rsn = function (n, xi, omega, alpha) {
  delta = alpha/sqrt(1 + alpha^{2})
  tn = matrix(rnorm(2*n), nrow=2)
  chi = abs(tn[1,])
  nrv = tn[2,]
  z = delta * chi + sqrt(1 - delta^{2}) * nrv
  
  y = xi + omega * z
  return(y)
}

dsn = function(x, xi, omega, alpha, logt = FALSE){
  z = (x-xi)/omega
  dx = (2/omega)*dnorm(z)*pnorm(alpha*z)
  if (logt == TRUE) {dx = log(dx)}
  return(dx)
}

dscainvchisq = function(tau, a_tau2, b_tau2, logt = FALSE){
  ab_tau = a_tau2*b_tau2
  a_tau_half = a_tau2/2
  ab_tau_half = ab_tau/2
  
  dtau = a_tau_half * log(ab_tau_half) - ab_tau_half/tau -
    lgamma(a_tau_half) - (1+a_tau_half) * log(tau)
  if (logt == FALSE) {dtau = exp(dtau)}
  
  return(dtau)
}

# Define function to update mean parameter for bivariate confounders
update_conf_bin = function(X_temp, a_bin, b_bin) {
  n_temp = length(X_temp)
  sum_X_temp = sum(X_temp)
  
  a_bin_new = a_bin + sum_X_temp
  b_bin_new = b_bin + n_temp - sum_X_temp
  
  bin_prop = c(rbeta(1, a_bin_new, b_bin_new))
  
  return(bin_prop)
  # return(list(bin_prop = bin_prop, a_bin_new = a_bin_new, b_bin_new = b_bin_new))
}

# Define function to update mean parameter for categorical confounders
update_conf_cat = function(X_temp, a_cat) {
  
  n_cat = length(a_cat)
  
  a_cat_new = count_cpp(X_temp,n_cat) + a_cat
  cat_prop  = c(rdirichlet_cpp(1,a_cat_new))
  
  return(cat_prop)
}

# Define function to update mean parameter and variance 
# parameter for continuous confounders, N-Inv-chi^2 model
update_conf_con = function(X_temp, a_var, b_var, a_mean, b_mean) {
  n_temp = length(X_temp)
  sum_X_temp = sum(X_temp)
  mean_X_temp = sum_X_temp/n_temp
  
  if (n_temp == 1) { 
    sum2_temp = 0
    weighted_mean2_temp = (b_mean * n_temp/(b_mean + n_temp)) * (X_temp - a_mean)^{2}
  } else {
    sum2_temp = (n_temp - 1) * var(X_temp)
    weighted_mean2_temp = (b_mean * n_temp/(b_mean + n_temp)) * (mean_X_temp - a_mean)^{2}
  }
  
  # Update the continuous confounders variance
  a_var_new = a_var + n_temp
  b_var_new = (a_var * b_var + sum2_temp + weighted_mean2_temp)/a_var_new
  var_prop = c(rscainvchisq_cpp(1, a_var_new, b_var_new))
  
  # Update the continuous confounders mean
  b_mean_new = b_mean + n_temp
  a_mean_new = (a_mean * b_mean + sum_X_temp)/b_mean_new
  var_new = var_prop/b_mean_new
  mean_prop = c(rnorm(1, a_mean_new, sqrt(var_new)))
  
  return(list(mean_par = mean_prop, var_par = var_prop))
}

# Define function to update outcome regression coefficients
update_reg_con = function(Y_temp,matX_temp,beta_curr,
                           a_sig2,b_sig2,Binv_beta,aBinv_beta) {
  n_temp = length(Y_temp)
  p_temp = length(beta_curr)
  matX_temp = matrix(matX_temp, ncol = p_temp)
  matXtX_temp = t(matX_temp) %*% matX_temp
  matXtY_temp = t(matX_temp) %*% Y_temp

  # Update the outcome regression parameters (sig2, beta)
  SSE_sig2 = Y_temp - matX_temp %*% beta_curr
  SSE_sig2 = sum(SSE_sig2^{2})
  a_sig2_new = a_sig2 + n_temp
  b_sig2_new = (a_sig2 * b_sig2 + SSE_sig2)/a_sig2_new
  sig2_prop = c(rscainvchisq_cpp(1, a_sig2_new, b_sig2_new))

  Binv_beta_new = matXtX_temp/sig2_prop + Binv_beta
  B_beta_new = inv_cpp(Binv_beta_new)
  a_beta_new = B_beta_new %*% (matXtY_temp/sig2_prop + aBinv_beta)
  beta_prop = c(rmvn_cpp(1, a_beta_new, B_beta_new))

  return(list(beta_par = beta_prop, sig2_par = sig2_prop))
}

# Define function to generate datasets
generate_data = function(N,Scenario){
  
  # Generate treatment, Z
  PI_z=0.5
  Z=rbinom(N,1,PI_z)
  
  z0=0
  z1=1
  ind_z0=(Z==z0)
  ind_z1=(Z==z1)
  n0=sum(ind_z0)
  n1=sum(ind_z1)
  
  prob_Y1=0.6
  I_Y = numeric(N)
  I_Y[ind_z0] = rbinom(n0, 1, prob_Y1)
  I_Y[ind_z1] = rbinom(n1, 1, prob_Y1)
  
  cond123=(sum(Scenario==c(1,2,3))==1)
  cond456=(sum(Scenario==c(4,5,6))==1)
  cond789=(sum(Scenario==c(7,8,9))==1)
  cond1012=(sum(Scenario==c(10,11,12))==1)
  
  if (cond123||cond456){
    # Generate confounders, C
    mu_C1=0;sig_C1=(3)
    mu_C2=0;sig_C2=(4)
    C1=rnorm(N,mu_C1,sig_C1)
    C2=rnorm(N,mu_C2,sig_C2)
    C=cbind(C1,C2)
    
    # Generate a post-treatment confounder, V
    if (cond123){
      MU_V0=1.3+0.6*C1-0.7*C2
      MU_V1=1.4-0.5*C1+0.3*C2
      MU_V01=cbind(MU_V0,MU_V1)
      rho_V01=0.3
      SIG_V01=matrix(c(3,rho_V01*sqrt(3*10),rho_V01*sqrt(3*10),10),nrow=2)
      V01=t(sapply(1:N,function(l)rmvn_cpp(1,MU_V01[l,],SIG_V01)))
    } else if (cond456){
      MU_V0=1.3+0.6*C1-0.7*C2
      MU_V1=1.4-0.5*C1+0.3*C2
      MU_V01=cbind(MU_V0,MU_V1)
      rho_V01=0.3
      V01=t(rmvgamma(N,log(1+exp(MU_V01)),1,rho_V01))
    }
    V0=V01[,1]
    V1=V01[,2]
    V=ifelse(Z==z0,V0,V1)
    
    if (Scenario==1||Scenario==4){
      # ------------------------------- Scenario 1 & 4 -------------------------------
      # Generate mediation, M
      omega_M=3
      alpha_M=10
      xi_M0=1.0+1.7*z0+0.5*V0+0.4*C1+0.9*C2
      xi_M1=1.0+1.7*z1+0.5*V1+0.4*C1+0.9*C2
      M0=rsn(N,xi=xi_M0,omega=omega_M,alpha=alpha_M)
      M1=rsn(N,xi=xi_M1,omega=omega_M,alpha=alpha_M)
      
      # Generate outcome, Y
      sig_Y1=(1.5)
      sig_Y2=(0.5)
      mu1_Y_z0m0=5+2.5*z0+1.8*M0+0.3*V0-1.2*C1+0.3*C2
      mu2_Y_z0m0=-5-1.5*z0-1.0*M0-0.7*V0+0.4*C1+0.3*C2
      
      mu1_Y_z1m0=5+2.5*z1+1.8*M0+0.3*V1-1.2*C1+0.3*C2
      mu2_Y_z1m0=-5-1.5*z1-1.0*M0-0.7*V1+0.4*C1+0.3*C2
      
      mu1_Y_z1m1=5+2.5*z1+1.8*M1+0.3*V1-1.2*C1+0.3*C2
      mu2_Y_z1m1=-5-1.5*z1-1.0*M1-0.7*V1+0.4*C1+0.3*C2
      
      Y_z0m0=ifelse((I_Y==1),rnorm(N,mu1_Y_z0m0,sig_Y1),rnorm(N,mu2_Y_z0m0,sig_Y2))
      Y_z1m0=ifelse((I_Y==1),rnorm(N,mu1_Y_z1m0,sig_Y1),rnorm(N,mu2_Y_z1m0,sig_Y2))
      Y_z1m1=ifelse((I_Y==1),rnorm(N,mu1_Y_z1m1,sig_Y1),rnorm(N,mu2_Y_z1m1,sig_Y2))
      
    } else if (Scenario==2||Scenario==5){
      # ------------------------------- Scenario 2 & 5 -------------------------------
      # Generate mediation, M
      omega_M=3
      alpha_M=10
      xi_M0=1.0+1.7*z0+0.5*V0+0.4*C1+0.3*C2
      xi_M1=1.0+1.7*z1+0.5*V1+0.4*C1+0.3*C2
      M0=rsn(N,xi=xi_M0,omega=omega_M,alpha=alpha_M)
      M1=rsn(N,xi=xi_M1,omega=omega_M,alpha=alpha_M)
      
      # Generate outcome, Y
      sig_Y1=(1.5)
      sig_Y2=(0.5)
      
      mu1_Y_z0m0=5+2.5*z0+1.8*M0+1.0*z0*M0+0.3*V0-0.4*C1+0.3*C2
      mu2_Y_z0m0=-5-1.5*z0-1.0*M0-0.5*z0*M0-0.7*V0+0.4*C1+0.3*C2
      
      mu1_Y_z1m0=5+2.5*z1+1.8*M0+1.0*z1*M0+0.3*V1-0.4*C1+0.3*C2
      mu2_Y_z1m0=-5-1.5*z1-1.0*M0-0.5*z1*M0-0.7*V1+0.4*C1+0.3*C2
      
      mu1_Y_z1m1=5+2.5*z1+1.8*M1+1.0*z1*M1+0.3*V1-0.4*C1+0.3*C2
      mu2_Y_z1m1=-5-1.5*z1-1.0*M1-0.5*z1*M1-0.7*V1+0.4*C1+0.3*C2
      
      Y_z0m0=ifelse((I_Y==1),rnorm(N,mu1_Y_z0m0,sig_Y1),rnorm(N,mu2_Y_z0m0,sig_Y2))
      Y_z1m0=ifelse((I_Y==1),rnorm(N,mu1_Y_z1m0,sig_Y1),rnorm(N,mu2_Y_z1m0,sig_Y2))
      Y_z1m1=ifelse((I_Y==1),rnorm(N,mu1_Y_z1m1,sig_Y1),rnorm(N,mu2_Y_z1m1,sig_Y2))
      
    } else if (Scenario==3||Scenario==6){
      # ------------------------------- Scenario 3 & 6 -------------------------------
      # Generate mediation, M
      omega_M=1
      alpha_M=7
      xi_M0=-1.0+0.2*z0+0.2*V0+0.1*C1+0.3*C2
      xi_M1=-1.0+0.2*z1+0.2*V1+0.1*C1+0.3*C2
      M0=rsn(N,xi=xi_M0,omega=omega_M,alpha=alpha_M)
      M1=rsn(N,xi=xi_M1,omega=omega_M,alpha=alpha_M)
      
      # Generate outcome, Y
      step=function(x,knot){ifelse(x>knot,x-knot,0)}
      sig_Y=(0.2)
      
      mu_Y_z0m0=5+2.5*z0+0.2*step(M0,1.0)+1.0*M0^{2}+0.3*V0+0.4*C1+0.3*C2
      mu_Y_z1m0=5+2.5*z1+0.2*step(M0,1.0)+1.0*M0^{2}+0.3*V1+0.4*C1+0.3*C2
      mu_Y_z1m1=5+2.5*z1+0.2*step(M1,1.0)+1.0*M1^{2}+0.3*V1+0.4*C1+0.3*C2
      
      Y_z0m0=rnorm(N,mu_Y_z0m0,sig_Y)
      Y_z1m0=rnorm(N,mu_Y_z1m0,sig_Y)
      Y_z1m1=rnorm(N,mu_Y_z1m1,sig_Y)
    }
    
  } else if (cond789||cond1012){
    # ------------------------------------------------------------------------------
    # Generate confounders, C
    pC11=3
    pC12=3
    pC13=3
    p_C1=pC11+pC12+pC13
    p_C2=6
    
    one_p2=rep(1,p_C2)
    zero_p2=rep(0,p_C2)
    DIAG_p2=diag(p_C2)
    MU_p2=zero_p2
    SIG_p2=0.7*DIAG_p2+0.3*one_p2%*%t(one_p2)
    C2=t(rmvn_cpp(N,MU_p2,SIG_p2))
    
    c_C=0.7*C2[,3]-0.4*C2[,5]*C2[,6]
    b_C=0.6*C2[,3]*C2[,4]-0.2*(C2[,5])^{2}
    a_C=invlogit(2*(C2[,1]-2)^{2}-2*(C2[,2]+1)^{2})
    p_C=invlogit(a_C*invlogit(b_C)+(1-a_C)*invlogit(c_C))
    
    PI1=0.05
    PI2=0.5
    PI3=p_C
    C1=matrix(NA,nrow=N,ncol=p_C1)
    for(q in 1:pC11){
      C1[,q]=rbinom(N,1,PI1)
    }
    for(q in 1:pC12){
      C1[,(pC11+q)]=rbinom(N,1,PI2)
    }
    for(q in 1:pC13){
      C1[,(pC11+pC12+q)]=rbinom(N,1,PI3)
    }
    C=cbind(C1,C2)
    
    # ------------------------------------------------------------------------------
    # Generate a post-treatment confounder, V
    if (cond789){
      MU_V0=1.3+0.5*C[,10]-0.7*C[,11]+0.3*C[,12]
      MU_V1=1.4+0.3*C[,13]-0.2*C[,14]-0.4*C[,15]
      MU_V01=cbind(MU_V0,MU_V1)
      rho_V01=0.3
      SIG_V01=matrix(c(3,rho_V01*sqrt(3*10),rho_V01*sqrt(3*10),10),nrow=2)
      V01=t(sapply(1:N,function(l)rmvn_cpp(1,MU_V01[l,],SIG_V01)))
    } else if (cond1012){
      MU_V0=0.5+0.5*C[,10]-0.7*C[,11]+0.3*C[,12]
      MU_V1=1.4+0.3*C[,13]-0.2*C[,14]-0.4*C[,15]
      MU_V01=cbind(MU_V0,MU_V1)
      rho_V01=0.3
      V01=t(rmvgamma(N,log(1+exp(MU_V01)),1,rho_V01))
    }
    V0=V01[,1]
    V1=V01[,2]
    V=ifelse(Z==z0,V0,V1)
    
    if (Scenario==7||Scenario==10){
      # ------------------------------- Scenario 7 & 10 -------------------------------
      # Generate mediation, M
      omega_M=3
      alpha_M=10
      xi_M0=-1.0+1.7*z0+0.5*V0+0.3*C[,1]+0.1*C[,4]-0.2*C[,7]-0.4*C[,10]+0.6*C[,13]
      xi_M1=-1.0+1.7*z1+0.5*V1+0.3*C[,1]+0.1*C[,4]-0.2*C[,7]-0.4*C[,10]+0.6*C[,13]
      M0=rsn(N,xi=xi_M0,omega=omega_M,alpha=alpha_M)
      M1=rsn(N,xi=xi_M1,omega=omega_M,alpha=alpha_M)
      
      # Generate outcome, Y
      sig_Y1=(1.5)
      sig_Y2=(0.5)
      
      mu1_Y_z0m0=5+2.5*z0+1.8*M0+0.3*V0+0.1*C[,2]+0.3*C[,5]-0.4*C[,8]-0.2*C[,11]+0.6*C[,14]
      mu2_Y_z0m0=-5-1.5*z0-1.0*M0-0.7*V0+0.3*C[,3]+0.1*C[,6]-0.2*C[,9]+0.6*C[,12]-0.4*C[,15]
      
      mu1_Y_z1m0=5+2.5*z1+1.8*M0+0.3*V1+0.1*C[,2]+0.3*C[,5]-0.4*C[,8]-0.2*C[,11]+0.6*C[,14]
      mu2_Y_z1m0=-5-1.5*z1-1.0*M0-0.7*V1+0.3*C[,3]+0.1*C[,6]-0.2*C[,9]+0.6*C[,12]-0.4*C[,15]
      
      mu1_Y_z1m1=5+2.5*z1+1.8*M1+0.3*V1+0.1*C[,2]+0.3*C[,5]-0.4*C[,8]-0.2*C[,11]+0.6*C[,14]
      mu2_Y_z1m1=-5-1.5*z1-1.0*M1-0.7*V1+0.3*C[,3]+0.1*C[,6]-0.2*C[,9]+0.6*C[,12]-0.4*C[,15]
      
      Y_z0m0=ifelse((I_Y==1),rnorm(N,mu1_Y_z0m0,sig_Y1),rnorm(N,mu2_Y_z0m0,sig_Y2))
      Y_z1m0=ifelse((I_Y==1),rnorm(N,mu1_Y_z1m0,sig_Y1),rnorm(N,mu2_Y_z1m0,sig_Y2))
      Y_z1m1=ifelse((I_Y==1),rnorm(N,mu1_Y_z1m1,sig_Y1),rnorm(N,mu2_Y_z1m1,sig_Y2))
      
    } else if (Scenario==8||Scenario==11){
      # ------------------------------- Scenario 8 & 11 -------------------------------
      # Generate mediation, M
      omega_M=3
      alpha_M=10
      xi_M0=-1.0+1.7*z0+0.5*V0+0.3*C[,1]+0.1*C[,4]-0.2*C[,7]-0.4*C[,10]+0.6*C[,13]
      xi_M1=-1.0+1.7*z1+0.5*V1+0.3*C[,1]+0.1*C[,4]-0.2*C[,7]-0.4*C[,10]+0.6*C[,13]
      M0=rsn(N,xi=xi_M0,omega=omega_M,alpha=alpha_M)
      M1=rsn(N,xi=xi_M1,omega=omega_M,alpha=alpha_M)
      
      # Generate outcome, Y
      sig_Y1=(1.5)
      sig_Y2=(0.5)
      
      mu1_Y_z0m0=5+2.5*z0+1.8*M0+1.0*z0*M0+0.3*V0+0.1*C[,2]+0.3*C[,5]-0.4*C[,8]-0.2*C[,11]+0.6*C[,14]
      mu2_Y_z0m0=-5-1.5*z0-1.0*M0-1.5*z0*M0+0.7*V0+0.3*C[,3]+0.1*C[,6]-0.2*C[,9]+0.6*C[,12]-0.4*C[,15]
      
      mu1_Y_z1m0=5+2.5*z1+1.8*M0+1.0*z1*M0+0.3*V1+0.1*C[,2]+0.3*C[,5]-0.4*C[,8]-0.2*C[,11]+0.6*C[,14]
      mu2_Y_z1m0=-5-1.5*z1-1.0*M0-0.5*z1*M0+0.7*V1+0.3*C[,3]+0.1*C[,6]-0.2*C[,9]+0.6*C[,12]-0.4*C[,15]
      
      mu1_Y_z1m1=5+2.5*z1+1.8*M1+1.0*z1*M1+0.3*V1+0.1*C[,2]+0.3*C[,5]-0.4*C[,8]-0.2*C[,11]+0.6*C[,14]
      mu2_Y_z1m1=-5-1.5*z1-1.0*M1-0.5*z1*M1+0.7*V1+0.3*C[,3]+0.1*C[,6]-0.2*C[,9]+0.6*C[,12]-0.4*C[,15]
      
      Y_z0m0=ifelse((I_Y==1),rnorm(N,mu1_Y_z0m0,sig_Y1),rnorm(N,mu2_Y_z0m0,sig_Y2))
      Y_z1m0=ifelse((I_Y==1),rnorm(N,mu1_Y_z1m0,sig_Y1),rnorm(N,mu2_Y_z1m0,sig_Y2))
      Y_z1m1=ifelse((I_Y==1),rnorm(N,mu1_Y_z1m1,sig_Y1),rnorm(N,mu2_Y_z1m1,sig_Y2))
      
    } else if (Scenario==9||Scenario==12){
      # ------------------------------- Scenario 9 & 12 -------------------------------
      # Generate mediation, M
      omega_M=1
      alpha_M=7
      xi_M0=-1.0+0.2*z0+0.2*V0+0.5*C[,1]+0.6*C[,4]-0.2*C[,7]-0.4*C[,10]+0.6*C[,13]
      xi_M1=-1.0+0.2*z1+0.2*V1+0.5*C[,1]+0.6*C[,4]-0.2*C[,7]-0.4*C[,10]+0.6*C[,13]
      M0=rsn(N,xi=xi_M0,omega=omega_M,alpha=alpha_M)
      M1=rsn(N,xi=xi_M1,omega=omega_M,alpha=alpha_M)
      
      # Generate outcome, Y
      step=function(x,knot){ifelse(x>knot,x-knot,0)}
      sig_Y=(0.2)
      
      mu_Y_z0m0=5+2.5*z0+0.2*step(M0,1.0)+1.0*M0^{2}+0.3*V0+0.2*C[,2]+0.3*C[,5]-0.4*C[,8]-0.2*C[,11]+0.6*C[,14]
      mu_Y_z1m0=5+2.5*z1+0.2*step(M0,1.0)+1.0*M0^{2}+0.3*V1+0.2*C[,2]+0.3*C[,5]-0.4*C[,8]-0.2*C[,11]+0.6*C[,14]
      mu_Y_z1m1=5+2.5*z1+0.2*step(M1,1.0)+1.0*M1^{2}+0.3*V1+0.2*C[,2]+0.3*C[,5]-0.4*C[,8]-0.2*C[,11]+0.6*C[,14]
      
      Y_z0m0=rnorm(N,mu_Y_z0m0,sig_Y)
      Y_z1m0=rnorm(N,mu_Y_z1m0,sig_Y)
      Y_z1m1=rnorm(N,mu_Y_z1m1,sig_Y)
    }
  }
  
  M = ifelse(Z==z0,M0,M1)
  Y = ifelse(Z==z0,Y_z0m0,Y_z1m1)
  
  # NIE_true = mean(Y_z1m1) - mean(Y_z1m0)
  # NDE_true = mean(Y_z1m0) - mean(Y_z0m0)
  # ATE_true = mean(Y_z1m1) - mean(Y_z0m0)
  NIE_true = mean(Y_z1m1)
  NDE_true = mean(Y_z1m0)
  ATE_true = mean(Y_z0m0)
  CE_true = c(NIE_true,NDE_true,ATE_true)
  
  return(list(Y = Y, M = M, V = V, Z = Z, C = C, CE_true = CE_true))
}

# Define function to fit the EDP model
EDPmediationMCMC = function(Y, M, V, Z, C, 
                            gibbs_iter = 2e4, gibbs_burnin = 2e4, gibbs_thin = 100,
                            D = 1e4, num_MC = 1e4, num_MC_prior = 1e5){
  
  # Define of constants
  
  # Define the number of observations for each dataset
  N = length(Y)
  
  # Define interation check
  iter_check = floor(gibbs_iter/20)
  
  # Define the number of Markov Chain Monte Carlo (MCMC) draws in Gibbs Sampler
  gibbs_total = gibbs_iter + gibbs_burnin
  
  # Define design matrix
  X = cbind(Z, C)
  matX = cbind(1, Z, C)
  matV = cbind(1, Z, V, C)
  # matM = cbind(1, Z, M, V, C)
  matM = cbind(1, Z, M*Z, M, V, C)
  
  matXtX = t(matX) %*% matX
  INVmatXtX = inv_cpp(matXtX)
  matVtV = t(matV) %*% matV
  INVmatVtV = inv_cpp(matVtV)
  matMtM = t(matM) %*% matM
  INVmatMtM = inv_cpp(matMtM)
  
  # Calculate the number of parameters
  py = ncol(matM)
  pm = ncol(matV)
  pv = ncol(matX)
  px = ncol(X)
  
  # Hyperparameters
  # for modeling binary confounders
  a_pi = 1
  b_pi = 1
  
  # for modeling continuous confounders
  a_tau2 = 1
  b_tau2 = 1
  b_mu = 0.5
  a_mu = 0
  
  # Parameters in posttreatment confounder regression
  # for posttreatment confounder regression coefficients
  glm_Vbeta = glm(V ~ matX-1)
  dof_Vbeta = glm_Vbeta$df.residual
  dev_Vbeta = glm_Vbeta$deviance
  var_Vbeta = dev_Vbeta/dof_Vbeta
  a_Vbeta = ifelse(is.na(glm_Vbeta$coef),0,glm_Vbeta$coef)
  B_Vbeta = var_Vbeta * INVmatXtX
  c_Vbeta = N/5
  B_Vbeta = c_Vbeta * B_Vbeta
  Binv_Vbeta = matXtX/(var_Vbeta * c_Vbeta)
  aBinv_Vbeta = Binv_Vbeta %*% a_Vbeta
  
  # for posttreatment confounder regression coefficients
  a_Vsig2 = 1
  b_Vsig2 = var_Vbeta/a_Vsig2
  
  # Parameters in mediator regression
  # for mediator regression coefficients
  glm_Mbeta = glm(M ~ matV-1)
  dof_Mbeta = glm_Mbeta$df.residual
  dev_Mbeta = glm_Mbeta$deviance
  var_Mbeta = dev_Mbeta/dof_Mbeta
  a_Mbeta = ifelse(is.na(glm_Mbeta$coef),0,glm_Mbeta$coef)
  B_Mbeta = var_Mbeta * INVmatVtV
  c_Mbeta = N/5
  B_Mbeta = c_Mbeta * B_Mbeta
  Binv_Mbeta = matVtV/(var_Mbeta * c_Mbeta)
  aBinv_Mbeta = Binv_Mbeta %*% a_Mbeta
  
  # for mediator regression coefficients
  a_Msig2 = 1
  b_Msig2 = var_Mbeta/a_Msig2
  
  # Parameters in outcome regression
  # for outcome regression coefficients
  glm_Ybeta = glm(Y ~ matM-1)
  dof_Ybeta = glm_Ybeta$df.residual
  dev_Ybeta = glm_Ybeta$deviance
  var_Ybeta = dev_Ybeta/dof_Ybeta
  a_Ybeta = ifelse(is.na(glm_Ybeta$coef),0,glm_Ybeta$coef)
  B_Ybeta = var_Ybeta * INVmatMtM
  c_Ybeta = N/5
  B_Ybeta = c_Ybeta * B_Ybeta
  Binv_Ybeta = matMtM/(var_Ybeta * c_Ybeta)
  aBinv_Ybeta = Binv_Ybeta %*% a_Ybeta
  
  # for outcome regression coefficients
  a_Ysig2 = 1
  b_Ysig2 = var_Ybeta/a_Ysig2
  
  # Set initial values ------------------------------------------------------------
  YbetaPars = c(rmvn_cpp(1, a_Ybeta, B_Ybeta))
  MbetaPars = c(rmvn_cpp(1, a_Mbeta, B_Mbeta))
  VbetaPars = c(rmvn_cpp(1, a_Vbeta, B_Vbeta))
  
  # Update values of betas in outcome regressions
  # Update values of betas in mediator regressions
  # Update values of betas in posttreatment regressions
  Yregression_list = update_reg_con(Y,matM,YbetaPars,a_Ysig2,b_Ysig2,Binv_Ybeta,aBinv_Ybeta)
  Mregression_list = update_reg_con(M,matV,MbetaPars,a_Msig2,b_Msig2,Binv_Mbeta,aBinv_Mbeta)
  Vregression_list = update_reg_con(V,matX,VbetaPars,a_Vsig2,b_Vsig2,Binv_Vbeta,aBinv_Vbeta)
  Ysig2Pars = Yregression_list$sig2_par
  YbetaPars = Yregression_list$beta_par
  Msig2Pars = Mregression_list$sig2_par
  MbetaPars = Mregression_list$beta_par
  Vsig2Pars = Vregression_list$sig2_par
  VbetaPars = Vregression_list$beta_par
  
  # parameter of baseline confounders (CrhoPars) - BB model
  CrhoPars = c(rdirichlet_cpp(1, rep(1,N)))
  
  # parameter of binary treatment (ZpiPars)
  ZpiPars = update_conf_bin(Z, a_pi, b_pi)
  
  # Make vectors to store draws from Gibbs Sampler
  n_store = floor(gibbs_iter/gibbs_thin)
  
  YbetaLists = list(NA)
  Ysig2Lists = list(NA)
  MbetaLists = list(NA)
  Msig2Lists = list(NA)
  VbetaLists = list(NA)
  Vsig2Lists = list(NA)
  CrhoLists  = list(NA)
  ZpiLists   = list(NA)
  
  count_it = 1
  # End initial values ------------------------------------------------------------
  for (gibbs_reps in 1:gibbs_total) {
    # Update values of betas in outcome regressions
    # Update values of betas in mediator regressions
    # Update values of betas in posttreatment regressions
    Yregression_list = update_reg_con(Y,matM,YbetaPars,a_Ysig2,b_Ysig2,Binv_Ybeta,aBinv_Ybeta)
    Mregression_list = update_reg_con(M,matV,MbetaPars,a_Msig2,b_Msig2,Binv_Mbeta,aBinv_Mbeta)
    Vregression_list = update_reg_con(V,matX,VbetaPars,a_Vsig2,b_Vsig2,Binv_Vbeta,aBinv_Vbeta)
    Ysig2Pars = Yregression_list$sig2_par
    YbetaPars = Yregression_list$beta_par
    Msig2Pars = Mregression_list$sig2_par
    MbetaPars = Mregression_list$beta_par
    Vsig2Pars = Vregression_list$sig2_par
    VbetaPars = Vregression_list$beta_par
    
    # parameter of binary treatment (ZpiPars)
    ZpiPars = update_conf_bin(Z, a_pi, b_pi)
    
    # End update of all parameters in Parametric-BB model--------------------------------
    
    if (gibbs_reps < gibbs_burnin) {
    } else if (gibbs_reps == gibbs_burnin) {
      cat("Bur-In End",gibbs_reps,"Time:",date(),"\n")
    } else if (gibbs_reps > gibbs_burnin) {
      if (gibbs_reps %% gibbs_thin == 0) {
        # parameter of baseline confounders (CrhoPars)
        CrhoPars = c(rdirichlet_cpp(1, rep(1,N)))
        
        YbetaLists[[count_it]] = YbetaPars
        Ysig2Lists[[count_it]] = Ysig2Pars
        MbetaLists[[count_it]] = MbetaPars
        Msig2Lists[[count_it]] = Msig2Pars
        VbetaLists[[count_it]] = VbetaPars
        Vsig2Lists[[count_it]] = Vsig2Pars
        CrhoLists[[count_it]]  = CrhoPars
        ZpiLists[[count_it]]   = ZpiPars
        
        count_it = count_it + 1
      }
      
      if (gibbs_reps %% iter_check == 0) {
        cat("Gibbs Iteration",(gibbs_reps-gibbs_burnin),"(",(gibbs_reps-gibbs_burnin)/gibbs_iter*100,"%)","Time:",date(),"\n")
      }
    }
  }
  
  # constants
  constants = list(N = N, D = D, num_MC = num_MC, n_MCMC = n_store)
  
  # priors
  priors = list(a_Ysig2 = a_Ysig2, b_Ysig2 = b_Ysig2, 
                a_Ybeta = a_Ybeta, B_Ybeta = B_Ybeta,
                a_Msig2 = a_Msig2, b_Msig2 = b_Msig2, 
                a_Mbeta = a_Mbeta, B_Mbeta = B_Mbeta,
                a_Vsig2 = a_Vsig2, b_Vsig2 = b_Vsig2, 
                a_Vbeta = a_Vbeta, B_Vbeta = B_Vbeta)
  
  # MCMC Posteriors
  MCMCposteriors = list(YbetaLists = YbetaLists,
                        Ysig2Lists = Ysig2Lists,
                        MbetaLists = MbetaLists,
                        Msig2Lists = Msig2Lists,
                        VbetaLists = VbetaLists,
                        Vsig2Lists = Vsig2Lists,
                        CrhoLists  = CrhoLists,
                        ZpiLists   = ZpiLists)
  
  # MCMC results
  MCMCresult = list(priors = priors,
                    constants = constants,
                    MCMCposteriors = MCMCposteriors)
  
  return(MCMCresult)
}

EDPmediationPOST = function(MCMCresult, esttype, save_cluster){
  
  # 
  z0 = 0
  z1 = 1
  
  # Significance Level alpha
  level = 0.05
  quantile_alpha = c(level/2,1-level/2)
  
  # the number of Gibbs Sampler stored
  n_MCMC = MCMCresult$constants$n_MCMC
  
  # Define interation check
  iter_check = floor(n_MCMC/10)
  
  NIE = numeric(n_MCMC)
  NDE = numeric(n_MCMC)
  ATE = numeric(n_MCMC)
  
  N = MCMCresult$constants$N
  D = MCMCresult$constants$D
  
  a_Ysig2 = MCMCresult$priors$a_Ysig2
  b_Ysig2 = MCMCresult$priors$b_Ysig2
  a_Ybeta = MCMCresult$priors$a_Ybeta
  B_Ybeta = MCMCresult$priors$B_Ybeta
  a_Msig2 = MCMCresult$priors$a_Msig2
  b_Msig2 = MCMCresult$priors$b_Msig2
  a_Mbeta = MCMCresult$priors$a_Mbeta
  B_Mbeta = MCMCresult$priors$B_Mbeta
  a_Vsig2 = MCMCresult$priors$a_Vsig2
  b_Vsig2 = MCMCresult$priors$b_Vsig2
  a_Vbeta = MCMCresult$priors$a_Vbeta
  B_Vbeta = MCMCresult$priors$B_Vbeta
  
  rho_mc = runif(n_MCMC)
  
  # End initial values ------------------------------------------------------------
  
  # G-computation -----------------------------------------------------------------
  for (post_reps in 1:n_MCMC) {
    rho = rho_mc[post_reps]
    
    YbetaPars = MCMCresult$MCMCposteriors$YbetaLists[[post_reps]]
    Ysig2Pars = MCMCresult$MCMCposteriors$Ysig2Lists[[post_reps]]
    MbetaPars = MCMCresult$MCMCposteriors$MbetaLists[[post_reps]]
    Msig2Pars = MCMCresult$MCMCposteriors$Msig2Lists[[post_reps]]
    VbetaPars = MCMCresult$MCMCposteriors$VbetaLists[[post_reps]]
    Vsig2Pars = MCMCresult$MCMCposteriors$Vsig2Lists[[post_reps]]
    CrhoPars  = MCMCresult$MCMCposteriors$CrhoLists[[post_reps]]
    ZpiPars   = MCMCresult$MCMCposteriors$ZpiLists[[post_reps]]
    
    # -------------------------------------------------------------------------------
    C_mc = C[sample(N,D,replace=TRUE,prob=CrhoPars),]
    
    # -------------------------------------------------------------------------------
    matX0_mc = cbind(1, z0, C_mc)
    matX1_mc = cbind(1, z1, C_mc)
    V0_mc = rnorm(D, matX0_mc %*% VbetaPars, sqrt(Vsig2Pars))
    V1_mc = rnorm(D, matX1_mc %*% VbetaPars, sqrt(Vsig2Pars))
    
    Vr_condMU_mc = (matX0_mc %*% VbetaPars) + rho * (V1_mc - matX1_mc %*% VbetaPars)
    Vr_condSIG2_mc = (1 - rho^{2}) * Vsig2Pars
    
    Vr_mc = rnorm(N, Vr_condMU_mc, sqrt(Vr_condSIG2_mc))
    
    # -------------------------------------------------------------------------------
    matV0_x0_mc = cbind(1, z0, V0_mc, C_mc)
    matV1_x1_mc = cbind(1, z1, V1_mc, C_mc)
    matVr_x0_mc = cbind(1, z0, Vr_mc, C_mc)
    
    # -------------------------------------------------------------------------------
    M0_mc = rnorm(D, matV0_x0_mc %*% MbetaPars, sqrt(Msig2Pars))
    Mr_mc = rnorm(D, matVr_x0_mc %*% MbetaPars, sqrt(Msig2Pars))
    M1_mc = rnorm(D, matV1_x1_mc %*% MbetaPars, sqrt(Msig2Pars))
    
    matM0_x0v0_mc = cbind(1, z0, M0_mc*z0, M0_mc, V0_mc, C_mc)
    matMr_x1v1_mc = cbind(1, z1, Mr_mc*z1, Mr_mc, V1_mc, C_mc)
    matM1_x1v1_mc = cbind(1, z1, M1_mc*z1, M1_mc, V1_mc, C_mc)
    
    # -------------------------------------------------------------------------------
    E_Y_x0v0m0_mc = matM0_x0v0_mc %*% YbetaPars
    E_Y_x1v1mr_mc = matMr_x1v1_mc %*% YbetaPars
    E_Y_x1v1m1_mc = matM1_x1v1_mc %*% YbetaPars
    
    E_Y_z0m0 = sum(E_Y_x0v0m0_mc)/D
    E_Y_z1mr = sum(E_Y_x1v1mr_mc)/D
    E_Y_z1m1 = sum(E_Y_x1v1m1_mc)/D
    
    NIE[post_reps] = E_Y_z1m1 - E_Y_z1mr
    NDE[post_reps] = E_Y_z1mr - E_Y_z0m0
    ATE[post_reps] = E_Y_z1m1 - E_Y_z0m0
    
    if (post_reps %% iter_check == 0){
      cat("Post-Processing",post_reps,"(",(post_reps/n_MCMC)*100,"%)","Time:",date(),"\n")
    }
  }
  # for (post_reps in 1:n_MCMC) {
  #   rho = rho_mc[post_reps]
  #   
  #   YbetaPars = MCMCresult$MCMCposteriors$YbetaLists[[post_reps]]
  #   Ysig2Pars = MCMCresult$MCMCposteriors$Ysig2Lists[[post_reps]]
  #   MbetaPars = MCMCresult$MCMCposteriors$MbetaLists[[post_reps]]
  #   Msig2Pars = MCMCresult$MCMCposteriors$Msig2Lists[[post_reps]]
  #   VbetaPars = MCMCresult$MCMCposteriors$VbetaLists[[post_reps]]
  #   Vsig2Pars = MCMCresult$MCMCposteriors$Vsig2Lists[[post_reps]]
  #   CrhoPars  = MCMCresult$MCMCposteriors$CrhoLists[[post_reps]]
  #   ZpiPars   = MCMCresult$MCMCposteriors$ZpiLists[[post_reps]]
  #   
  #   # -------------------------------------------------------------------------------
  #   matX0_mc = cbind(1, z0, scaC)
  #   matX1_mc = cbind(1, z1, scaC)
  #   V0_mc = rnorm(N, matX0_mc %*% VbetaPars, sqrt(Vsig2Pars))
  #   V1_mc = rnorm(N, matX1_mc %*% VbetaPars, sqrt(Vsig2Pars))
  #   
  #   Vr_condMU_mc = (matX0_mc %*% VbetaPars) + rho * (V1_mc - matX1_mc %*% VbetaPars)
  #   Vr_condSIG2_mc = (1 - rho^{2}) * Vsig2Pars
  #   
  #   Vr_mc = rnorm(N, Vr_condMU_mc, sqrt(Vr_condSIG2_mc))
  #   
  #   # -------------------------------------------------------------------------------
  #   matV0_x0_mc = cbind(1, z0, V0_mc, scaC)
  #   matV1_x1_mc = cbind(1, z1, V1_mc, scaC)
  #   matVr_x0_mc = cbind(1, z0, Vr_mc, scaC)
  #   
  #   # -------------------------------------------------------------------------------
  #   M0_mc = rnorm(N, matV0_x0_mc %*% MbetaPars, sqrt(Msig2Pars))
  #   Mr_mc = rnorm(N, matVr_x0_mc %*% MbetaPars, sqrt(Msig2Pars))
  #   M1_mc = rnorm(N, matV1_x1_mc %*% MbetaPars, sqrt(Msig2Pars))
  #   
  #   matM0_x0v0_mc = cbind(1, z0, M0_mc*z0, M0_mc, V0_mc, scaC)
  #   matMr_x1v1_mc = cbind(1, z1, Mr_mc*z1, Mr_mc, V1_mc, scaC)
  #   matM1_x1v1_mc = cbind(1, z1, M1_mc*z1, M1_mc, V1_mc, scaC)
  #   
  #   # -------------------------------------------------------------------------------
  #   E_Y_x0v0m0_mc = matM0_x0v0_mc %*% YbetaPars
  #   E_Y_x1v1mr_mc = matMr_x1v1_mc %*% YbetaPars
  #   E_Y_x1v1m1_mc = matM1_x1v1_mc %*% YbetaPars
  #   
  #   E_Y_z0m0 = sum(E_Y_x0v0m0_mc * CrhoPars)
  #   E_Y_z1mr = sum(E_Y_x1v1mr_mc * CrhoPars)
  #   E_Y_z1m1 = sum(E_Y_x1v1m1_mc * CrhoPars)
  #   
  #   NIE[post_reps] = E_Y_z1m1 - E_Y_z1mr
  #   NDE[post_reps] = E_Y_z1mr - E_Y_z0m0
  #   ATE[post_reps] = E_Y_z1m1 - E_Y_z0m0
  #   
  #   if (post_reps %% iter_check == 0){
  #     cat("Post-Processing",post_reps,"(",(post_reps/n_MCMC)*100,"%)","Time:",date(),"\n")
  #   }
  # }
  # End G-computation -------------------------------------------------------------
  
  # Estimates: mean or median
  if (esttype == "median"){
    NIE_est_mc = median(NIE)
    NDE_est_mc = median(NDE)
    ATE_est_mc = median(ATE)
  } else {
    NIE_est_mc = mean(NIE)
    NDE_est_mc = mean(NDE)
    ATE_est_mc = mean(ATE)
  }
  
  # Calculate median of posterior for NIE
  NIE_sd_mc = sd(NIE)
  NIE_quantile_mc = quantile(NIE, quantile_alpha)
  NIE_CIlength_mc = abs(diff(NIE_quantile_mc))
  NIE_quantile025_mc = min(NIE_quantile_mc)
  NIE_quantile975_mc = max(NIE_quantile_mc)
  
  # Calculate median of posterior for NDE
  NDE_sd_mc = sd(NDE)
  NDE_quantile_mc = quantile(NDE, quantile_alpha)
  NDE_CIlength_mc = abs(diff(NDE_quantile_mc))
  NDE_quantile025_mc = min(NDE_quantile_mc)
  NDE_quantile975_mc = max(NDE_quantile_mc)
  
  # Calculate median of posterior for ATE
  ATE_sd_mc = sd(ATE)
  ATE_quantile_mc = quantile(ATE, quantile_alpha)
  ATE_CIlength_mc = abs(diff(ATE_quantile_mc))
  ATE_quantile025_mc = min(ATE_quantile_mc)
  ATE_quantile975_mc = max(ATE_quantile_mc)
  
  # Save the results
  NIE_result_mc = cbind(NIE_est_mc, NIE_sd_mc, NIE_quantile025_mc, NIE_quantile975_mc, NIE_CIlength_mc)
  colnames(NIE_result_mc) = c("estimates","sd", "quantile025", "quantile975", "CIlength95")
  rownames(NIE_result_mc) = c("")
  
  NDE_result_mc = cbind(NDE_est_mc, NDE_sd_mc, NDE_quantile025_mc, NDE_quantile975_mc, NDE_CIlength_mc)
  colnames(NDE_result_mc) = c("estimates","sd", "quantile025", "quantile975", "CIlength95")
  rownames(NDE_result_mc) = c("")
  
  ATE_result_mc = cbind(ATE_est_mc, ATE_sd_mc, ATE_quantile025_mc, ATE_quantile975_mc, ATE_CIlength_mc)
  colnames(ATE_result_mc) = c("estimates","sd", "quantile025", "quantile975", "CIlength95")
  rownames(ATE_result_mc) = c("")
  
  POSTresult = list(NIE_result_mc = NIE_result_mc,
                    NDE_result_mc = NDE_result_mc,
                    ATE_result_mc = ATE_result_mc,
                    NIE = NIE, NDE = NDE, ATE = ATE)
  
  return(POSTresult)
}

EDPmediation = function(Y, M, V, Z, C, 
                        gibbs_iter = 2e4, gibbs_burnin = 2e4, gibbs_thin = 100,
                        D = 2e4, num_MC = 1e4, num_MC_prior = 1e5, 
                        esttype = c("median","mean"), save_cluster = FALSE){
  
  MCMCresult = EDPmediationMCMC(Y, M, V, Z, C, 
                                gibbs_iter, gibbs_burnin, gibbs_thin,
                                D, num_MC, num_MC_prior)
  
  POSTresult = EDPmediationPOST(MCMCresult, esttype, save_cluster)
  
  if (save_cluster == TRUE){
    result = c(MCMCresult, POSTresult)
  } else {
    result = POSTresult
  }
  
  return(result)
}

# End function definitions