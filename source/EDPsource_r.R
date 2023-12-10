# Define functions --------------------------------------------------------------
VARprior = function(x){
  varprior = (diff(range(x)/4)^{2})
  return(varprior)
}

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
      mu1_Y_z0m0=5+2.5*z0+1.8*M0+1.3*V0-1.2*C1+0.3*C2
      mu2_Y_z0m0=-5-1.5*z0-1.0*M0-0.7*V0+0.4*C1+0.3*C2
      
      mu1_Y_z1m0=5+2.5*z1+1.8*M0+1.3*V1-1.2*C1+0.3*C2
      mu2_Y_z1m0=-5-1.5*z1-1.0*M0-0.7*V1+0.4*C1+0.3*C2
      
      mu1_Y_z1m1=5+2.5*z1+1.8*M1+1.3*V1-1.2*C1+0.3*C2
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
      
      mu1_Y_z0m0=5+2.5*z0+1.8*M0+1.0*z0*M0+1.3*V0-0.4*C1+0.3*C2
      mu2_Y_z0m0=-5-1.5*z0-1.0*M0-0.5*z0*M0-0.7*V0+0.4*C1+0.3*C2
      
      mu1_Y_z1m0=5+2.5*z1+1.8*M0+1.0*z1*M0+1.3*V1-0.4*C1+0.3*C2
      mu2_Y_z1m0=-5-1.5*z1-1.0*M0-0.5*z1*M0-0.7*V1+0.4*C1+0.3*C2
      
      mu1_Y_z1m1=5+2.5*z1+1.8*M1+1.0*z1*M1+1.3*V1-0.4*C1+0.3*C2
      mu2_Y_z1m1=-5-1.5*z1-1.0*M1-0.5*z1*M1-0.7*V1+0.4*C1+0.3*C2
      
      Y_z0m0=ifelse((I_Y==1),rnorm(N,mu1_Y_z0m0,sig_Y1),rnorm(N,mu2_Y_z0m0,sig_Y2))
      Y_z1m0=ifelse((I_Y==1),rnorm(N,mu1_Y_z1m0,sig_Y1),rnorm(N,mu2_Y_z1m0,sig_Y2))
      Y_z1m1=ifelse((I_Y==1),rnorm(N,mu1_Y_z1m1,sig_Y1),rnorm(N,mu2_Y_z1m1,sig_Y2))
      
    } else if (Scenario==3||Scenario==6){
      # ------------------------------- Scenario 3 & 6 -------------------------------
      # Generate mediation, M
      omega_M=1
      alpha_M=7
      xi_M0=-1.5+0.5*z0+0.1*V0+0.1*C1+0.3*C2
      xi_M1=-1.5+0.5*z1+0.1*V1+0.1*C1+0.3*C2
      M0=rsn(N,xi=xi_M0,omega=omega_M,alpha=alpha_M)
      M1=rsn(N,xi=xi_M1,omega=omega_M,alpha=alpha_M)
      
      # Generate outcome, Y
      step=function(x,knot){ifelse(x>knot,x-knot,0)}
      sig_Y=(0.2)
      
      mu_Y_z0m0=5+2.5*z0+0.2*step(M0,0.4)+0.6*M0^{2}+0.3*V0+0.4*C1+0.3*C2
      mu_Y_z1m0=5+2.5*z1+0.2*step(M0,0.4)+0.6*M0^{2}+0.3*V1+0.4*C1+0.3*C2
      mu_Y_z1m1=5+2.5*z1+0.2*step(M1,0.4)+0.6*M1^{2}+0.3*V1+0.4*C1+0.3*C2
      
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
      MU_V0=1.3+0.5*C[,10]-0.7*C[,11]+0.3*C[,12]
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
      
      mu1_Y_z0m0=5+2.5*z0+1.8*M0+1.3*V0+0.1*C[,2]+0.3*C[,5]-0.4*C[,8]-0.2*C[,11]+0.6*C[,14]
      mu2_Y_z0m0=-5-1.5*z0-1.0*M0-0.7*V0+0.3*C[,3]+0.1*C[,6]-0.2*C[,9]+0.6*C[,12]-0.4*C[,15]
      
      mu1_Y_z1m0=5+2.5*z1+1.8*M0+1.3*V1+0.1*C[,2]+0.3*C[,5]-0.4*C[,8]-0.2*C[,11]+0.6*C[,14]
      mu2_Y_z1m0=-5-1.5*z1-1.0*M0-0.7*V1+0.3*C[,3]+0.1*C[,6]-0.2*C[,9]+0.6*C[,12]-0.4*C[,15]
      
      mu1_Y_z1m1=5+2.5*z1+1.8*M1+1.3*V1+0.1*C[,2]+0.3*C[,5]-0.4*C[,8]-0.2*C[,11]+0.6*C[,14]
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
      
      mu1_Y_z0m0=5+2.5*z0+1.8*M0+1.0*z0*M0+1.3*V0+0.1*C[,2]+0.3*C[,5]-0.4*C[,8]-0.2*C[,11]+0.6*C[,14]
      mu2_Y_z0m0=-5-1.5*z0-1.0*M0-0.5*z0*M0-0.7*V0+0.3*C[,3]+0.1*C[,6]-0.2*C[,9]+0.6*C[,12]-0.4*C[,15]
      
      mu1_Y_z1m0=5+2.5*z1+1.8*M0+1.0*z1*M0+1.3*V1+0.1*C[,2]+0.3*C[,5]-0.4*C[,8]-0.2*C[,11]+0.6*C[,14]
      mu2_Y_z1m0=-5-1.5*z1-1.0*M0-0.5*z1*M0-0.7*V1+0.3*C[,3]+0.1*C[,6]-0.2*C[,9]+0.6*C[,12]-0.4*C[,15]
      
      mu1_Y_z1m1=5+2.5*z1+1.8*M1+1.0*z1*M1+1.3*V1+0.1*C[,2]+0.3*C[,5]-0.4*C[,8]-0.2*C[,11]+0.6*C[,14]
      mu2_Y_z1m1=-5-1.5*z1-1.0*M1-0.5*z1*M1-0.7*V1+0.3*C[,3]+0.1*C[,6]-0.2*C[,9]+0.6*C[,12]-0.4*C[,15]
      
      Y_z0m0=ifelse((I_Y==1),rnorm(N,mu1_Y_z0m0,sig_Y1),rnorm(N,mu2_Y_z0m0,sig_Y2))
      Y_z1m0=ifelse((I_Y==1),rnorm(N,mu1_Y_z1m0,sig_Y1),rnorm(N,mu2_Y_z1m0,sig_Y2))
      Y_z1m1=ifelse((I_Y==1),rnorm(N,mu1_Y_z1m1,sig_Y1),rnorm(N,mu2_Y_z1m1,sig_Y2))
      
    } else if (Scenario==9||Scenario==12){
      # ------------------------------- Scenario 9 & 12 -------------------------------
      # Generate mediation, M
      omega_M=1
      alpha_M=7
      xi_M0=-0.7+0.2*z0+0.1*V0+0.5*C[,1]+0.3*C[,4]+0.2*C[,7]-0.1*C[,10]+0.4*C[,13]
      xi_M1=-0.7+0.2*z1+0.1*V1+0.5*C[,1]+0.3*C[,4]+0.2*C[,7]-0.1*C[,10]+0.4*C[,13]
      M0=rsn(N,xi=xi_M0,omega=omega_M,alpha=alpha_M)
      M1=rsn(N,xi=xi_M1,omega=omega_M,alpha=alpha_M)
      
      # Generate outcome, Y
      step=function(x,knot){ifelse(x>knot,x-knot,0)}
      sig_Y=(0.2)
      
      mu_Y_z0m0=5+2.5*z0+0.2*step(M0,0.4)+0.6*M0^{2}+0.3*V0+0.2*C[,2]+0.1*C[,5]-0.4*C[,8]-0.2*C[,11]+0.3*C[,14]
      mu_Y_z1m0=5+2.5*z1+0.2*step(M0,0.4)+0.6*M0^{2}+0.3*V1+0.2*C[,2]+0.1*C[,5]-0.4*C[,8]-0.2*C[,11]+0.3*C[,14]
      mu_Y_z1m1=5+2.5*z1+0.2*step(M1,0.4)+0.6*M1^{2}+0.3*V1+0.2*C[,2]+0.1*C[,5]-0.4*C[,8]-0.2*C[,11]+0.3*C[,14]
      
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
  
  # Define the number of posttreatment confounders
  p_V = length(V)/N
  
  # Define the number of treatment variables
  p_Z = length(Z)/N
  
  # Define the number of confounders
  p_C = ncol(C)
  
  # Define the number of binary confounders
  p_C1 = sum(apply(C, 2, function(x) { all(x %in% 0:1) }))
  
  # Define the number of continuous confounders
  p_C2 = p_C - p_C1
  
  # Define interation check
  iter_check = floor(gibbs_iter/20)
  
  # Define the number of Markov Chain Monte Carlo (MCMC) draws in Gibbs Sampler
  gibbs_total = gibbs_iter + gibbs_burnin
  
  # Standardize continuous confounders (important when choosing priors)
  if (p_C1>0 && p_C2>0){
    C1 = C[,1:p_C1]
    C2 = C[,(p_C1+1):(p_C1+p_C2)]
    scaC2 = apply(C2, 2, scale)
    scaC = cbind(C1, scaC2)
  } else if (p_C2>0){
    scaC = apply(C, 2, scale)
  } else {
    scaC = C
  }
  
  # Define design matrix
  X = cbind(Z, scaC)
  matX = cbind(1, Z, scaC)
  matV = cbind(1, Z, V, scaC)
  # matM = cbind(1, Z, M, V, scaC)
  matM = cbind(1, Z, M*Z, M, V, scaC)
  
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
  # a_Vbeta = numeric(pv)
  # B_Vbeta = diag(pv)
  # c_Vbeta = N/5
  # B_Vbeta = c_Vbeta * B_Vbeta
  # Binv_Vbeta = diag(pv)/c_Vbeta
  # aBinv_Vbeta = Binv_Vbeta %*% a_Vbeta
  
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
  # a_Mbeta = numeric(pm)
  # B_Mbeta = diag(pm)
  # c_Mbeta = N/5
  # B_Mbeta = c_Mbeta * B_Mbeta
  # Binv_Mbeta = diag(pm)/c_Mbeta
  # aBinv_Mbeta = Binv_Mbeta %*% a_Mbeta
  
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
  # a_Ybeta = numeric(py)
  # B_Ybeta = diag(py)
  # c_Ybeta = N/5
  # B_Ybeta = c_Ybeta * B_Ybeta
  # Binv_Ybeta = diag(py)/c_Ybeta
  # aBinv_Ybeta = Binv_Ybeta %*% a_Ybeta
  
  # for outcome regression coefficients
  a_Ysig2 = 1
  b_Ysig2 = var_Ybeta/a_Ysig2
  
  # for concentration parameters
  # randg: shape and scale at default in Rcpp & RcppArmadillo
  # rgamma: shape and rate at default in R
  # a_theta = (p_Y + p_M)
  a_theta = 1
  b_theta = 1
  # a_omega = (p_Z + p_C)
  a_omega = 1
  b_omega = 1
  
  # Set initial values ------------------------------------------------------------
  
  # Set initial values for cluster membership
  Syx = matrix(nrow = N, ncol = 2, 1)
  
  # Use Ky-means to find Ky = 2 whole data clusters
  Ky = 2
  Sy = kmeans(cbind(Y, M, V, Z, scaC), Ky)$cluster
  Syx[,1] = Sy
  
  # Make Kx = 2 x-clusters within each y-cluster
  Kx = 2
  Sx = numeric(N)
  for (k in 1:Ky) {
    ind_Sy = (Sy == k)
    Sx[ind_Sy] = kmeans(cbind(V, Z, scaC)[ind_Sy,], Kx)$cluster
  }
  Syx[,2] = Sx
  
  # Determine unique clusters from the cluster membership variable, Syx
  unique_Syx = unique(Syx)
  
  # Sort clusters
  unique_Syx = unique_Syx[order(unique_Syx[,1], unique_Syx[,2]), , drop = FALSE]
  
  # Make vector of y-clusters
  unique_Sy = sort(unique(Sy))
  
  # Calculate the number of y-clusters
  Ky = length(unique_Sy)
  
  # Update values of etas in mediator regressions  
  # Update values of betas in outcome regressions
  Msig2Pars = matrix(nrow = 1, ncol = Ky)
  MbetaPars = rmvn_cpp(Ky, a_Mbeta, B_Mbeta)
  Ysig2Pars = matrix(nrow = 1, ncol = Ky)
  YbetaPars = rmvn_cpp(Ky, a_Ybeta, B_Ybeta)
  
  for (k in unique_Sy) {
    ind_Sy = (Sy == k)
    
    Mregression_list = update_reg_con(M[ind_Sy],matV[ind_Sy,],MbetaPars[,k],
                                       a_Msig2,b_Msig2,Binv_Mbeta,aBinv_Mbeta)
    Yregression_list = update_reg_con(Y[ind_Sy],matM[ind_Sy,],YbetaPars[,k],
                                       a_Ysig2,b_Ysig2,Binv_Ybeta,aBinv_Ybeta)
    Msig2Pars[,k] = Mregression_list$sig2_par
    MbetaPars[,k] = Mregression_list$beta_par
    Ysig2Pars[,k] = Yregression_list$sig2_par
    YbetaPars[,k] = Yregression_list$beta_par
  }
  
  # Update confounders parameters
  # Within each cluster parameters are stored in a long vector n_k x 1
  # Binary confounders have 1 parameter (XpiPars)
  # Continuous confounders have 2 parameters (XmuPars, Xtau2Pars)
  
  # Initialize vectors n_k is the total number of clusters (x and y clusters)
  Kyx = nrow(unique(Syx))
  
  # Update parameters for binary confounders prior is beta(a_pi,b_pi), posterior is beta
  Vsig2Pars = matrix(nrow = 1, ncol = Kyx)
  VbetaPars = rmvn_cpp(Kyx, a_Vbeta, B_Vbeta)
  ZpiPars = matrix(nrow = p_Z, ncol = Kyx)
  XpiPars = matrix(nrow = p_C1, ncol = Kyx)
  XmuPars = matrix(nrow = p_C2, ncol = Kyx)
  Xtau2Pars = matrix(nrow = p_C2, ncol = Kyx)
  
  # update parameters for continuous posttreatment confounders
  {
    count = 1
    for (k in unique_Sy) {
      ind_Sy = which(Sy == k)
      
      v_temp_k = V[ind_Sy]
      matX_temp_k = matrix(matX[ind_Sy,],ncol=pv)
      
      Sx_temp = Sx[ind_Sy]
      unique_Sx_temp = sort(unique(Sx_temp))
      for (r in unique_Sx_temp) {
        ind_Syx = (Sx_temp == r)
        v_temp_kr = v_temp_k[ind_Syx]
        matX_temp_kr = matrix(matX_temp_k[ind_Syx,],ncol=pv)
        
        Vregression_list = update_reg_con(v_temp_kr,matX_temp_kr,VbetaPars[,count],
                                           a_Vsig2,b_Vsig2,Binv_Vbeta,aBinv_Vbeta)
        Vsig2Pars[,count] = Vregression_list$sig2_par
        VbetaPars[,count] = Vregression_list$beta_par
        
        count = count + 1
      }
    }
  }
  
  # update parameters for binary treatment
  for (q in 1:p_Z) {
    count = 1
    z_temp_q = X[,q]
    for (k in unique_Sy) {
      ind_Sy = (Sy == k)
      
      Sx_temp = Sx[ind_Sy]
      z_temp_qk = z_temp_q[ind_Sy]
      unique_Sx_temp = sort(unique(Sx_temp))
      for (r in unique_Sx_temp) {
        z_temp_qkr = z_temp_qk[Sx_temp == r]
        
        # posterior is beta
        bivariate_pars = update_conf_bin(z_temp_qkr, a_pi, b_pi)
        ZpiPars[q, count] = bivariate_pars
        
        count = count + 1
      }
    }
  }
  
  # update parameters for binary confounders
  if (p_C1 > 0) {
    for (q in 1:p_C1) {
      # beta(1,1) prior
      count = 1
      x_temp_q = X[,(p_Z + q)]
      for (k in unique_Sy) {
        ind_Sy = (Sy == k)
        
        Sx_temp = Sx[ind_Sy]
        x_temp_qk = x_temp_q[ind_Sy]
        unique_Sx_temp = sort(unique(Sx_temp))
        for (r in unique_Sx_temp) {
          x_temp_qkr = x_temp_qk[Sx_temp == r]
          
          # posterior is beta
          bivariate_pars = update_conf_bin(x_temp_qkr, a_pi, b_pi)
          XpiPars[q, count] = bivariate_pars
          
          count = count + 1
        }
      }
    }
  }
  
  # update parameters for continuous confounders
  if (p_C2 > 0) {
    for (q in 1:p_C2) {
      # beta(1,1) prior
      count = 1
      x_temp_q = X[,(p_Z + p_C1 + q)]
      for (k in unique_Sy) {
        ind_Sy = (Sy == k)
        
        Sx_temp = Sx[ind_Sy]
        x_temp_qk = x_temp_q[ind_Sy]
        unique_Sx_temp = sort(unique(Sx_temp))
        for (r in unique_Sx_temp) {
          x_temp_qkr = x_temp_qk[Sx_temp == r]
          
          # posterior for mu. prior for mu|sigma^2 mean 0 prior sample size2
          continuous_pars = update_conf_con(x_temp_qkr, a_tau2, b_tau2, a_mu, b_mu)
          Xtau2Pars[q, count] = continuous_pars$var_par
          XmuPars[q, count] = continuous_pars$mean_par
          
          count = count + 1
        }
      }
    }
  }
  
  # Calculate f0_x and f0_y for use in cluster function in Gibbs Sampler
  # Average confounders distribution over prior for each x_i
  ab_tau = a_tau2 * b_tau2
  a_tau_new = (a_tau2 + 1)/2
  a_tau_half = a_tau2/2
  b_mu_ratio = b_mu/(b_mu+1)
  margin_part1 = (gamma(a_tau_new)/gamma(a_tau_half)) * sqrt(b_mu_ratio/pi) * (ab_tau)^{a_tau_half}
  
  f0_x_all = matrix(0, nrow = N, ncol = px)
  # Binary confounders (including treatment)
  # Beta-Binomial
  for (q in 1:p_Z) {
    f0_x_all[, q] = beta(a_pi + X[, q], b_pi - X[, q] + 1) # beta(a_pi,b_pi) = 1 in this case
  }
  
  if (p_C1>0){
    for (q in (p_Z + 1):(p_Z + p_C1)) {
      f0_x_all[, q] = beta(a_pi + X[, q], b_pi - X[, q] + 1) # beta(a_pi,b_pi) = 1 in this case
    }
  }
  # Continuous confounders
  if (p_C2>0){
    for (q in 1:p_C2) {
      x_temp = X[, (p_Z + p_C1 + q)]
      margin_part2 = (ab_tau + b_mu_ratio * (x_temp-a_mu)^{2})^(-a_tau_new)
      f0_x_all[, (p_Z + p_C1 + q)] = margin_part1 * margin_part2
    }
  }
  
  # Take product (confounders are assumed to be locally independent).
  # Result is vector of size N
  f0_x = apply(f0_x_all, 1, prod)
  
  # Find f0(V|Z,C,theta) integrating over prior for theta coefficient using
  # Monte Carlo integration
  f0_v = numeric(N)
  for (num in 1:num_MC_prior) {
    Vsig2_temp = rscainvchisq_cpp(1, a_Vsig2, b_Vsig2)
    Vbeta_temp = rmvn_cpp(1, a_Vbeta, B_Vbeta)
    matXVbeta_temp = matX %*% Vbeta_temp
    f0_v = f0_v + dnorm(V, matXVbeta_temp, sqrt(Vsig2_temp))
  }
  f0_v = f0_v/num_MC_prior
  
  # Find f0(M|V,Z,C,theta) integrating over prior for theta coefficient using
  # Monte Carlo integration
  f0_m = numeric(N)
  for (num in 1:num_MC_prior) {
    Msig2_temp = rscainvchisq_cpp(1, a_Msig2, b_Msig2)
    Mbeta_temp = rmvn_cpp(1, a_Mbeta, B_Mbeta)
    matVMbeta_temp = matV %*% Mbeta_temp
    f0_m = f0_m + dnorm(M, matVMbeta_temp, sqrt(Msig2_temp))
  }
  f0_m = f0_m/num_MC_prior
  # f0_m = f0_m * 2
  
  # Find f0(Y|M,V,Z,C,theta] integrating over prior for theta coefficient using
  # Monte Carlo integration
  f0_y = numeric(N)
  for (num in 1:num_MC_prior) {
    Ysig2_temp = rscainvchisq_cpp(1, a_Ysig2, b_Ysig2)
    Ybeta_temp = rmvn_cpp(1, a_Ybeta, B_Ybeta)
    matMYbeta_temp = matM %*% Ybeta_temp
    f0_y = f0_y + dnorm(Y, matMYbeta_temp, sqrt(Ysig2_temp))
  }
  f0_y = f0_y/num_MC_prior
  
  # Will need the following for calculating causal effects using MC integration.
  # Any calculations involving X are done in the Gibbs Sampler
  pi_prior = a_pi/(a_pi+b_pi)
  f0_z0_mc = matrix(1-pi_prior, nrow = D, ncol = p_Z)
  f0_z1_mc = matrix(pi_prior, nrow = D, ncol = p_Z)
  f0_z0_mc = apply(f0_z0_mc, 1, prod)
  f0_z1_mc = apply(f0_z1_mc, 1, prod)
  
  Vsig2Prior_mc = rscainvchisq_cpp(num_MC, a_Vsig2, b_Vsig2)
  VbetaPrior_mc = rmvn_cpp(num_MC, a_Vbeta, B_Vbeta)
  Msig2Prior_mc = rscainvchisq_cpp(num_MC, a_Msig2, b_Msig2)
  MbetaPrior_mc = rmvn_cpp(num_MC, a_Mbeta, B_Mbeta)
  Ysig2Prior_mc = rscainvchisq_cpp(num_MC, a_Ysig2, b_Ysig2)
  YbetaPrior_mc = rmvn_cpp(num_MC, a_Ybeta, B_Ybeta)
  
  # Initialize alpha parameters
  alpha_theta = 2
  alpha_omega = 2
  
  # Make vectors to store draws from Gibbs Sampler
  n_store = floor(gibbs_iter/gibbs_thin)
  
  # Make lists to store draws from Gibbs Sampler
  alpha_theta_draws = numeric(n_store)
  alpha_omega_draws = numeric(n_store)
  
  YbetaLists = list(NA)
  Ysig2Lists = list(NA)
  MbetaLists = list(NA)
  Msig2Lists = list(NA)
  VbetaLists = list(NA)
  Vsig2Lists = list(NA)
  ZpiLists   = list(NA)
  XpiLists   = list(NA)
  XmuLists   = list(NA)
  Xtau2Lists = list(NA)
  
  n_kLists       = list(NA)
  n_rkLists      = list(NA)
  KyLists        = list(NA)
  KyxLists       = list(NA)
  max_Kx_SyLists = list(NA)
  
  count_it = 1
  # End initial values ------------------------------------------------------------
  for (gibbs_reps in 1:gibbs_total) {
    # gibbs_burnin
    # (gibbs_burnin+1)
    # gibbs_reps = 1
    # First draw each parameter for BNP model. Then calculate causal effect.
    
    # Update cluster membership --------------------------------------------
    cluster_res = EDPcluster_concon(N, p_Z, p_C1, p_C2, 
                                    Y, M, V, X, matM, matV, matX,
                                    Sy, Sx, unique_Syx,
                                    YbetaPars, Ysig2Pars, 
                                    MbetaPars, Msig2Pars,
                                    VbetaPars, Vsig2Pars,
                                    ZpiPars, XpiPars, XmuPars, Xtau2Pars,
                                    alpha_theta, alpha_omega,
                                    f0_y, f0_m, f0_v, f0_x,
                                    a_Ysig2, b_Ysig2, a_Ybeta, 
                                    B_Ybeta, Binv_Ybeta, aBinv_Ybeta,
                                    a_Msig2, b_Msig2, a_Mbeta, 
                                    B_Mbeta, Binv_Mbeta, aBinv_Mbeta,
                                    a_Vsig2, b_Vsig2, a_Vbeta, 
                                    B_Vbeta, Binv_Vbeta, aBinv_Vbeta, 
                                    a_pi, b_pi, a_mu, b_mu, a_tau2, b_tau2)
    
    # Store cluster membership output from cluster function
    Syx = cluster_res$Syx
    Sy  = Syx[,1]
    Sx  = Syx[,2]
    
    YbetaPars = cluster_res$YbetaPars
    Ysig2Pars = cluster_res$Ysig2Pars
    
    MbetaPars = cluster_res$MbetaPars
    Msig2Pars = cluster_res$Msig2Pars
    
    VbetaPars = cluster_res$VbetaPars
    Vsig2Pars = cluster_res$Vsig2Pars
    
    ZpiPars   = cluster_res$ZpiPars
    XpiPars   = cluster_res$XpiPars
    XmuPars   = cluster_res$XmuPars
    Xtau2Pars = cluster_res$Xtau2Pars
    
    # Make matrix of clusters
    unique_Syx = cluster_res$unique_Syx
    
    # Make vector of y-clusters
    unique_Sy = cluster_res$unique_Sy
    
    # Calculate the number of clusters
    Kyx = cluster_res$Kyx
    
    # Calculate the number of y-clusters
    Ky = cluster_res$Ky
    
    # Find the largest number of x clusters
    max_Sx = cluster_res$max_Sx
    
    # Store the number of x-clusters in each y-cluster in a vector
    max_Kx_Sy = cluster_res$max_Kx_Sy
    
    # Calculate the number of subjects in each y-cluster and store in vector n_k.
    # Calculate the number of subjects in each x-cluster and store in matrix n_rk.
    # Use k to index y-clusters and r to index x-clusters.
    n_k = cluster_res$n_k
    n_rk = cluster_res$n_rk
    
    # End update of cluster membership --------------------------------------------
    
    # Update values of theta coefficients in outcome regressions--------------------
    for (k in unique_Sy) {
      ind_Sy = (Sy == k)
      
      Mregression_list = update_reg_con(M[ind_Sy],matV[ind_Sy,],MbetaPars[,k],
                                         a_Msig2,b_Msig2,Binv_Mbeta,aBinv_Mbeta)
      Yregression_list = update_reg_con(Y[ind_Sy],matM[ind_Sy,],YbetaPars[,k],
                                         a_Ysig2,b_Ysig2,Binv_Ybeta,aBinv_Ybeta)
      Msig2Pars[,k] = Mregression_list$sig2_par
      MbetaPars[,k] = Mregression_list$beta_par
      Ysig2Pars[,k] = Yregression_list$sig2_par
      YbetaPars[,k] = Yregression_list$beta_par
    }
    # End update of theta coefficients in outcome regressions-------------------
    
    # Update cluster specific parameters for confounders distributions.----------
    # Binary confounders have 1 parameter (XpiPars).
    # Continuous confounders have2 parameters (XmuPars,Xtau2Pars).
    # Also store parameters in a matrix for use in MC integration later.
    # Vectors of parameters (ZpiPars, XpiPars, Xtau2Pars, XmuPars) is for cluster function
    VbetaPars_array = array(0, c(pv, (max_Sx + 2), (Ky + 2)))
    Vsig2Pars_array = array(0, c(p_V, (max_Sx + 2), (Ky + 2)))
    
    ZpiPars_array = array(0, c(p_Z, (max_Sx + 1), (Ky + 1)))
    XpiPars_array = array(0, c(p_C1, (max_Sx + 1), (Ky + 1)))
    XmuPars_array = array(0, c(p_C2, (max_Sx + 1), (Ky + 1)))
    Xtau2Pars_array = array(0, c(p_C2, (max_Sx + 1), (Ky + 1)))
    
    # update parameters for continuous posttreatment confounders
    {
      count = 1
      for (k in unique_Sy) {
        ind_Sy = which(Sy == k)
        
        v_temp_k = V[ind_Sy]
        matX_temp_k = matrix(matX[ind_Sy,],ncol=pv)
        
        Sx_temp = Sx[ind_Sy]
        unique_Sx_temp = sort(unique(Sx_temp))
        for (r in unique_Sx_temp) {
          ind_Syx = (Sx_temp == r)
          v_temp_kr = v_temp_k[ind_Syx]
          matX_temp_kr = matrix(matX_temp_k[ind_Syx,],ncol=pv)
          
          Vregression_list = update_reg_con(v_temp_kr,matX_temp_kr,VbetaPars[,count],
                                             a_Vsig2,b_Vsig2,Binv_Vbeta,aBinv_Vbeta)
          Vsig2Pars[,count] = Vregression_list$sig2_par
          VbetaPars[,count] = Vregression_list$beta_par
          Vsig2Pars_array[,r,k] = Vregression_list$sig2_par
          VbetaPars_array[,r,k] = Vregression_list$beta_par
          
          count = count + 1
        }
      }
    }
    
    # update parameters for binary treatment
    for (q in 1:p_Z) {
      count = 1
      x_temp_q = X[,q]
      for (k in unique_Sy) {
        ind_Sy = (Sy == k)
        
        Sx_temp = Sx[ind_Sy]
        unique_Sx_temp = sort(unique(Sx_temp))
        
        x_temp_qk = x_temp_q[ind_Sy]
        for (r in unique_Sx_temp) {
          x_temp_qkr = x_temp_qk[Sx_temp == r]
          
          # posterior is beta
          bivariate_pars = update_conf_bin(x_temp_qkr, a_pi, b_pi)
          ZpiPars_array[q, r, k] = bivariate_pars
          ZpiPars[q, count] = bivariate_pars
          
          count = count + 1
        }
      }
    }
    
    # update parameters for binary confounders
    if (p_C1 > 0) {
      for (q in 1:p_C1) {
        count = 1
        x_temp_q = X[,(p_Z + q)]
        for (k in unique_Sy) {
          ind_Sy = (Sy == k)
          
          Sx_temp = Sx[ind_Sy]
          unique_Sx_temp = sort(unique(Sx_temp))
          
          x_temp_qk = x_temp_q[ind_Sy]
          for (r in unique_Sx_temp) {
            x_temp_qkr = x_temp_qk[Sx_temp == r]
            
            # posterior is beta
            bivariate_pars = update_conf_bin(x_temp_qkr, a_pi, b_pi)
            XpiPars_array[q, r, k] = bivariate_pars
            XpiPars[q, count] = bivariate_pars
            
            count = count + 1
          }
        }
      }
    }
    
    # update parameters for continuous confounders
    if (p_C2 > 0) {
      for (q in 1:p_C2) {
        count = 1
        x_temp_q = X[,(p_Z + p_C1 + q)]
        for (k in unique_Sy) {
          ind_Sy = (Sy == k)
          
          Sx_temp = Sx[ind_Sy]
          unique_Sx_temp = sort(unique(Sx_temp))
          
          x_temp_qk = x_temp_q[ind_Sy]
          for (r in unique_Sx_temp) {
            x_temp_qkr = x_temp_qk[Sx_temp == r]
            
            # posterior for mu. prior for mu|sigma^2 mean 0 prior sample size2
            continuous_pars = update_conf_con(x_temp_qkr, a_tau2, b_tau2, a_mu, b_mu)
            Xtau2Pars_array[q, r, k] = continuous_pars$var_par
            XmuPars_array[q, r, k] = continuous_pars$mean_par
            
            Xtau2Pars[q, count] = continuous_pars$var_par
            XmuPars[q, count] = continuous_pars$mean_par
            
            count = count + 1
          }
        }
      }
    }
    # End update of cluster-specific parameters for confounders---------------------
    
    # Update concentration parameters--------------------------------------------
    alpha_theta = update_alpha_theta_cpp(N, Ky, alpha_theta, a_theta, b_theta)
    alpha_omega = update_alpha_omega_cpp(Kyx, Ky, Sy, unique_Sy, alpha_omega, a_omega, b_omega)
    
    # End update of concentration parameters-----------------------------------
    
    # End update of all parameters in BNP model--------------------------------
    
    if (gibbs_reps < gibbs_burnin) {
    } else if (gibbs_reps == gibbs_burnin) {
      cat("Bur-In End",gibbs_reps,"Time:",date(),"\n")
    } else if (gibbs_reps > gibbs_burnin) {
      if (gibbs_reps %% gibbs_thin == 0) {
        n_kLists[[count_it]]       = n_k
        n_rkLists[[count_it]]      = n_rk
        KyLists[[count_it]]        = Ky
        KyxLists[[count_it]]       = Kyx
        max_Kx_SyLists[[count_it]] = max_Kx_Sy
        
        YbetaLists[[count_it]] = YbetaPars
        Ysig2Lists[[count_it]] = Ysig2Pars
        MbetaLists[[count_it]] = MbetaPars
        Msig2Lists[[count_it]] = Msig2Pars
        VbetaLists[[count_it]] = VbetaPars_array
        Vsig2Lists[[count_it]] = Vsig2Pars_array
        ZpiLists[[count_it]]   = ZpiPars_array
        XpiLists[[count_it]]   = XpiPars_array
        XmuLists[[count_it]]   = XmuPars_array
        Xtau2Lists[[count_it]] = Xtau2Pars_array
        
        alpha_theta_draws[count_it] = alpha_theta
        alpha_omega_draws[count_it] = alpha_omega
        
        count_it = count_it + 1
      }
      
      if (gibbs_reps %% iter_check == 0) {
        cat("Gibbs Iteration",(gibbs_reps-gibbs_burnin),"(",(gibbs_reps-gibbs_burnin)/gibbs_iter*100,"%)","Time:",date(),"\n")
      }
    }
  }
  
  # constants
  constants = list(N = N, D = D, num_MC = num_MC, n_MCMC = n_store,
                   p_V = p_V, p_Z = p_Z, p_C1 = p_C1, p_C2 = p_C2)
  
  # priors
  priors = list(a_Ysig2 = a_Ysig2, b_Ysig2 = b_Ysig2, 
                a_Ybeta = a_Ybeta, B_Ybeta = B_Ybeta,
                a_Msig2 = a_Msig2, b_Msig2 = b_Msig2, 
                a_Mbeta = a_Mbeta, B_Mbeta = B_Mbeta,
                a_Vsig2 = a_Vsig2, b_Vsig2 = b_Vsig2, 
                a_Vbeta = a_Vbeta, B_Vbeta = B_Vbeta,
                a_pi = a_pi, b_pi = b_pi, 
                a_mu = a_mu, b_mu = b_mu, 
                a_tau2 = a_tau2, b_tau2 = b_tau2)
  
  # Precalculated values for Post-Processing steps
  MCpriors = list(f0_z0_mc = f0_z0_mc, 
                  f0_z1_mc = f0_z1_mc,
                  YbetaPrior_mc = YbetaPrior_mc, 
                  MbetaPrior_mc = MbetaPrior_mc, 
                  Msig2Prior_mc = Msig2Prior_mc,
                  VbetaPrior_mc = VbetaPrior_mc, 
                  Vsig2Prior_mc = Vsig2Prior_mc)
  
  # MCMC Posteriors
  MCMCposteriors = list(YbetaLists = YbetaLists,
                        Ysig2Lists = Ysig2Lists,
                        MbetaLists = MbetaLists,
                        Msig2Lists = Msig2Lists,
                        VbetaLists = VbetaLists,
                        Vsig2Lists = Vsig2Lists,
                        XpiLists   = XpiLists,
                        XmuLists   = XmuLists,
                        Xtau2Lists = Xtau2Lists,
                        ZpiLists   = ZpiLists,
                        n_kLists   = n_kLists,
                        n_rkLists  = n_rkLists,
                        KyLists    = KyLists,
                        KyxLists   = KyxLists,
                        max_Kx_SyLists    = max_Kx_SyLists,
                        alpha_theta_draws = alpha_theta_draws,
                        alpha_omega_draws = alpha_omega_draws)
  
  # MCMC results
  MCMCresult = list(priors = priors,
                    MCpriors = MCpriors,
                    constants = constants,
                    MCMCposteriors = MCMCposteriors)
  
  return(MCMCresult)
}

EDPmediationPOST = function(MCMCresult, esttype, save_cluster){
  
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
  num_MC = MCMCresult$constants$num_MC
  p_V = MCMCresult$constants$p_V
  p_Z = MCMCresult$constants$p_Z
  p_C1 = MCMCresult$constants$p_C1
  p_C2 = MCMCresult$constants$p_C2
  
  f0_z0_mc = MCMCresult$MCpriors$f0_z0_mc
  f0_z1_mc = MCMCresult$MCpriors$f0_z1_mc
  YbetaPrior_mc = MCMCresult$MCpriors$YbetaPrior_mc
  MbetaPrior_mc = MCMCresult$MCpriors$MbetaPrior_mc
  Msig2Prior_mc = MCMCresult$MCpriors$Msig2Prior_mc
  VbetaPrior_mc = MCMCresult$MCpriors$VbetaPrior_mc
  Vsig2Prior_mc = MCMCresult$MCpriors$Vsig2Prior_mc
  
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
  a_pi   = MCMCresult$priors$a_pi
  b_pi   = MCMCresult$priors$b_pi
  a_mu   = MCMCresult$priors$a_mu
  b_mu   = MCMCresult$priors$b_mu
  a_tau2 = MCMCresult$priors$a_tau2
  b_tau2 = MCMCresult$priors$b_tau2
  
  alpha_theta_draws = MCMCresult$MCMCposteriors$alpha_theta_draws
  alpha_omega_draws = MCMCresult$MCMCposteriors$alpha_omega_draws
  
  rho = runif(n_MCMC)
  
  # End initial values ------------------------------------------------------------
  for (post_reps in 1:n_MCMC) {
    Ky = MCMCresult$MCMCposteriors$KyLists[[post_reps]]
    Kyx = MCMCresult$MCMCposteriors$KyxLists[[post_reps]]
    n_k = MCMCresult$MCMCposteriors$n_kLists[[post_reps]]
    n_rk = MCMCresult$MCMCposteriors$n_rkLists[[post_reps]]
    max_Kx_Sy = MCMCresult$MCMCposteriors$max_Kx_SyLists[[post_reps]]
    
    rho_mc = rho[post_reps] * rep(1,D)
    YbetaPars = MCMCresult$MCMCposteriors$YbetaLists[[post_reps]]
    Ysig2Pars = MCMCresult$MCMCposteriors$Ysig2Lists[[post_reps]]
    MbetaPars = MCMCresult$MCMCposteriors$MbetaLists[[post_reps]]
    Msig2Pars = MCMCresult$MCMCposteriors$Msig2Lists[[post_reps]]
    VbetaPars_array = MCMCresult$MCMCposteriors$VbetaLists[[post_reps]]
    Vsig2Pars_array = MCMCresult$MCMCposteriors$Vsig2Lists[[post_reps]]
    ZpiPars_array   = MCMCresult$MCMCposteriors$ZpiLists[[post_reps]]
    XpiPars_array   = MCMCresult$MCMCposteriors$XpiLists[[post_reps]]
    XmuPars_array   = MCMCresult$MCMCposteriors$XmuLists[[post_reps]]
    Xtau2Pars_array = MCMCresult$MCMCposteriors$Xtau2Lists[[post_reps]]
    alpha_theta = alpha_theta_draws[post_reps]
    alpha_omega = alpha_omega_draws[post_reps]
    
    MCresult = EDPpostmc_concon(p_C1, p_C2,
                                n_k, n_rk, Ky, Kyx, max_Kx_Sy, 
                                D, num_MC, rho_mc,
                                YbetaPars, Ysig2Pars,
                                MbetaPars, Msig2Pars,
                                VbetaPars_array, Vsig2Pars_array,
                                ZpiPars_array, XpiPars_array, 
                                XmuPars_array, Xtau2Pars_array,
                                alpha_theta, alpha_omega, 
                                f0_z0_mc, f0_z1_mc,
                                YbetaPrior_mc, 
                                MbetaPrior_mc, Msig2Prior_mc,
                                VbetaPrior_mc, Vsig2Prior_mc,
                                a_Ysig2, b_Ysig2, a_Ybeta, B_Ybeta,
                                a_Msig2, b_Msig2, a_Mbeta, B_Mbeta,
                                a_Vsig2, b_Vsig2, a_Vbeta, B_Vbeta,
                                a_pi, b_pi, a_mu, b_mu, a_tau2, b_tau2)
    
    NIE[post_reps] = MCresult$NIE
    NDE[post_reps] = MCresult$NDE
    ATE[post_reps] = MCresult$ATE
    
    if (post_reps %% iter_check == 0){
      cat("Post-Processing",post_reps,"(",(post_reps/n_MCMC)*100,"%)","Time:",date(),"\n")
    }
  }
  
  # Calculate median of posterior for concentration parameters
  median_alpha_theta = median(alpha_theta_draws)
  median_alpha_omega = median(alpha_omega_draws)
  
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
  alpha_result = cbind(median_alpha_theta, median_alpha_omega)
  colnames(alpha_result) = c("alpha_theta","alpha_omega")
  rownames(alpha_result) = c("")
  
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

# EDPmediationPOST = function(MCMCresult, esttype, rho, save_cluster){
#   
#   # Significance Level alpha
#   level = 0.05
#   quantile_alpha = c(level/2,1-level/2)
#   
#   # the number of Gibbs Sampler stored
#   n_MCMC = MCMCresult$constants$n_MCMC
#   
#   # Define interation check
#   iter_check = floor(n_MCMC/10)
#   
#   NIE = numeric(n_MCMC)
#   NDE = numeric(n_MCMC)
#   ATE = numeric(n_MCMC)
#   
#   N = MCMCresult$constants$N
#   D = MCMCresult$constants$D
#   num_MC = MCMCresult$constants$num_MC
#   p_V = MCMCresult$constants$p_V
#   p_Z = MCMCresult$constants$p_Z
#   p_C1 = MCMCresult$constants$p_C1
#   p_C2 = MCMCresult$constants$p_C2
#   
#   f0_z0_mc = MCMCresult$MCpriors$f0_z0_mc
#   f0_z1_mc = MCMCresult$MCpriors$f0_z1_mc
#   YbetaPrior_mc = MCMCresult$MCpriors$YbetaPrior_mc
#   MbetaPrior_mc = MCMCresult$MCpriors$MbetaPrior_mc
#   Msig2Prior_mc = MCMCresult$MCpriors$Msig2Prior_mc
#   VbetaPrior_mc = MCMCresult$MCpriors$VbetaPrior_mc
#   Vsig2Prior_mc = MCMCresult$MCpriors$Vsig2Prior_mc
#   
#   a_Ysig2 = MCMCresult$priors$a_Ysig2
#   b_Ysig2 = MCMCresult$priors$b_Ysig2
#   a_Ybeta = MCMCresult$priors$a_Ybeta
#   B_Ybeta = MCMCresult$priors$B_Ybeta
#   a_Msig2 = MCMCresult$priors$a_Msig2
#   b_Msig2 = MCMCresult$priors$b_Msig2
#   a_Mbeta = MCMCresult$priors$a_Mbeta
#   B_Mbeta = MCMCresult$priors$B_Mbeta
#   a_Vsig2 = MCMCresult$priors$a_Vsig2
#   b_Vsig2 = MCMCresult$priors$b_Vsig2
#   a_Vbeta = MCMCresult$priors$a_Vbeta
#   B_Vbeta = MCMCresult$priors$B_Vbeta
#   a_pi   = MCMCresult$priors$a_pi
#   b_pi   = MCMCresult$priors$b_pi
#   a_mu   = MCMCresult$priors$a_mu
#   b_mu   = MCMCresult$priors$b_mu
#   a_tau2 = MCMCresult$priors$a_tau2
#   b_tau2 = MCMCresult$priors$b_tau2
#   
#   alpha_theta_draws = MCMCresult$MCMCposteriors$alpha_theta_draws
#   alpha_omega_draws = MCMCresult$MCMCposteriors$alpha_omega_draws
#   
#   if (rho == "trim10") {
#     rho = rtriangle(n_MCMC,-1,0,-1)
#   } else if (rho == "unifm10") {
#     rho = runif(n_MCMC,-1,0)
#   } else if (rho == "unifm11") {
#     rho = runif(n_MCMC,-1,1)
#   } else if (rho == "unif01") {
#     rho = runif(n_MCMC,0,1)
#   } else if (rho == "tri01") {
#     rho = rtriangle(n_MCMC,0,1,1)
#   } else {
#     rho = rep(rho, n_MCMC)
#   }
#   
#   # End initial values ------------------------------------------------------------
#   for (post_reps in 1:n_MCMC) {
#     Ky = MCMCresult$MCMCposteriors$KyLists[[post_reps]]
#     Kyx = MCMCresult$MCMCposteriors$KyxLists[[post_reps]]
#     n_k = MCMCresult$MCMCposteriors$n_kLists[[post_reps]]
#     n_rk = MCMCresult$MCMCposteriors$n_rkLists[[post_reps]]
#     max_Kx_Sy = MCMCresult$MCMCposteriors$max_Kx_SyLists[[post_reps]]
#     
#     rho_mc = rho[post_reps] * rep(1,D)
#     YbetaPars = MCMCresult$MCMCposteriors$YbetaLists[[post_reps]]
#     Ysig2Pars = MCMCresult$MCMCposteriors$Ysig2Lists[[post_reps]]
#     MbetaPars = MCMCresult$MCMCposteriors$MbetaLists[[post_reps]]
#     Msig2Pars = MCMCresult$MCMCposteriors$Msig2Lists[[post_reps]]
#     VbetaPars_array = MCMCresult$MCMCposteriors$VbetaLists[[post_reps]]
#     Vsig2Pars_array = MCMCresult$MCMCposteriors$Vsig2Lists[[post_reps]]
#     ZpiPars_array   = MCMCresult$MCMCposteriors$ZpiLists[[post_reps]]
#     XpiPars_array   = MCMCresult$MCMCposteriors$XpiLists[[post_reps]]
#     XmuPars_array   = MCMCresult$MCMCposteriors$XmuLists[[post_reps]]
#     Xtau2Pars_array = MCMCresult$MCMCposteriors$Xtau2Lists[[post_reps]]
#     alpha_theta = alpha_theta_draws[post_reps]
#     alpha_omega = alpha_omega_draws[post_reps]
#     
#     MCresult = EDPpostmc_concon(p_C1, p_C2,
#                                 n_k, n_rk, Ky, Kyx, max_Kx_Sy, 
#                                 D, num_MC, rho_mc,
#                                 YbetaPars, Ysig2Pars,
#                                 MbetaPars, Msig2Pars,
#                                 VbetaPars_array, Vsig2Pars_array,
#                                 ZpiPars_array, XpiPars_array, 
#                                 XmuPars_array, Xtau2Pars_array,
#                                 alpha_theta, alpha_omega, 
#                                 f0_z0_mc, f0_z1_mc,
#                                 YbetaPrior_mc, 
#                                 MbetaPrior_mc, Msig2Prior_mc,
#                                 VbetaPrior_mc, Vsig2Prior_mc,
#                                 a_Ysig2, b_Ysig2, a_Ybeta, B_Ybeta,
#                                 a_Msig2, b_Msig2, a_Mbeta, B_Mbeta,
#                                 a_Vsig2, b_Vsig2, a_Vbeta, B_Vbeta,
#                                 a_pi, b_pi, a_mu, b_mu, a_tau2, b_tau2)
#     
#     NIE[post_reps] = MCresult$NIE
#     NDE[post_reps] = MCresult$NDE
#     ATE[post_reps] = MCresult$ATE
#     
#     if (post_reps %% iter_check == 0){
#       cat("Post-Processing",post_reps,"(",(post_reps/n_MCMC)*100,"%)","Time:",date(),"\n")
#     }
#   }
#   
#   # Calculate median of posterior for concentration parameters
#   median_alpha_theta = median(alpha_theta_draws)
#   median_alpha_omega = median(alpha_omega_draws)
#   
#   # Estimates: mean or median
#   if (esttype == "median"){
#     NIE_est_mc = median(NIE)
#     NDE_est_mc = median(NDE)
#     ATE_est_mc = median(ATE)
#   } else {
#     NIE_est_mc = mean(NIE)
#     NDE_est_mc = mean(NDE)
#     ATE_est_mc = mean(ATE)
#   }
#   
#   # Calculate median of posterior for NIE
#   NIE_sd_mc = sd(NIE)
#   NIE_quantile_mc = quantile(NIE, quantile_alpha)
#   NIE_CIlength_mc = abs(diff(NIE_quantile_mc))
#   NIE_quantile025_mc = min(NIE_quantile_mc)
#   NIE_quantile975_mc = max(NIE_quantile_mc)
#   
#   # Calculate median of posterior for NDE
#   NDE_sd_mc = sd(NDE)
#   NDE_quantile_mc = quantile(NDE, quantile_alpha)
#   NDE_CIlength_mc = abs(diff(NDE_quantile_mc))
#   NDE_quantile025_mc = min(NDE_quantile_mc)
#   NDE_quantile975_mc = max(NDE_quantile_mc)
#   
#   # Calculate median of posterior for ATE
#   ATE_sd_mc = sd(ATE)
#   ATE_quantile_mc = quantile(ATE, quantile_alpha)
#   ATE_CIlength_mc = abs(diff(ATE_quantile_mc))
#   ATE_quantile025_mc = min(ATE_quantile_mc)
#   ATE_quantile975_mc = max(ATE_quantile_mc)
#   
#   # Save the results
#   alpha_result = cbind(median_alpha_theta, median_alpha_omega)
#   colnames(alpha_result) = c("alpha_theta","alpha_omega")
#   rownames(alpha_result) = c("")
#   
#   NIE_result_mc = cbind(NIE_est_mc, NIE_sd_mc, NIE_quantile025_mc, NIE_quantile975_mc, NIE_CIlength_mc)
#   colnames(NIE_result_mc) = c("estimates","sd", "quantile025", "quantile975", "CIlength95")
#   rownames(NIE_result_mc) = c("")
#   
#   NDE_result_mc = cbind(NDE_est_mc, NDE_sd_mc, NDE_quantile025_mc, NDE_quantile975_mc, NDE_CIlength_mc)
#   colnames(NDE_result_mc) = c("estimates","sd", "quantile025", "quantile975", "CIlength95")
#   rownames(NDE_result_mc) = c("")
#   
#   ATE_result_mc = cbind(ATE_est_mc, ATE_sd_mc, ATE_quantile025_mc, ATE_quantile975_mc, ATE_CIlength_mc)
#   colnames(ATE_result_mc) = c("estimates","sd", "quantile025", "quantile975", "CIlength95")
#   rownames(ATE_result_mc) = c("")
#   
#   POSTresult = list(NIE_result_mc = NIE_result_mc,
#                     NDE_result_mc = NDE_result_mc,
#                     ATE_result_mc = ATE_result_mc,
#                     NIE = NIE, NDE = NDE, ATE = ATE)
#   return(POSTresult)
# }

EDPmediationPOST_condC = function(MCMCresult, condC, q_condC, esttype, rho, save_cluster){
  
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
  num_MC = MCMCresult$constants$num_MC
  p_V = MCMCresult$constants$p_V
  p_Z = MCMCresult$constants$p_Z
  p_C1 = MCMCresult$constants$p_C1
  p_C2 = MCMCresult$constants$p_C2
  
  f0_z0_mc = MCMCresult$MCpriors$f0_z0_mc
  f0_z1_mc = MCMCresult$MCpriors$f0_z1_mc
  YbetaPrior_mc = MCMCresult$MCpriors$YbetaPrior_mc
  MbetaPrior_mc = MCMCresult$MCpriors$MbetaPrior_mc
  Msig2Prior_mc = MCMCresult$MCpriors$Msig2Prior_mc
  VbetaPrior_mc = MCMCresult$MCpriors$VbetaPrior_mc
  Vsig2Prior_mc = MCMCresult$MCpriors$Vsig2Prior_mc
  
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
  a_pi   = MCMCresult$priors$a_pi
  b_pi   = MCMCresult$priors$b_pi
  a_mu   = MCMCresult$priors$a_mu
  b_mu   = MCMCresult$priors$b_mu
  a_tau2 = MCMCresult$priors$a_tau2
  b_tau2 = MCMCresult$priors$b_tau2
  
  alpha_theta_draws = MCMCresult$MCMCposteriors$alpha_theta_draws
  alpha_omega_draws = MCMCresult$MCMCposteriors$alpha_omega_draws
  
  if (rho == "trim10") {
    rho = rtriangle(n_MCMC,-1,0,-1)
  } else if (rho == "unifm10") {
    rho = runif(n_MCMC,-1,0)
  } else if (rho == "unifm11") {
    rho = runif(n_MCMC,-1,1)
  } else if (rho == "unif01") {
    rho = runif(n_MCMC,0,1)
  } else if (rho == "tri01") {
    rho = rtriangle(n_MCMC,0,1,1)
  } else {
    rho = rep(rho, n_MCMC)
  }
  
  # conditional on confounder: C = condC
  condC = rep(condC, D)
  
  # End initial values ------------------------------------------------------------
  for (post_reps in 1:n_MCMC) {
    Ky = MCMCresult$MCMCposteriors$KyLists[[post_reps]]
    Kyx = MCMCresult$MCMCposteriors$KyxLists[[post_reps]]
    n_k = MCMCresult$MCMCposteriors$n_kLists[[post_reps]]
    n_rk = MCMCresult$MCMCposteriors$n_rkLists[[post_reps]]
    max_Kx_Sy = MCMCresult$MCMCposteriors$max_Kx_SyLists[[post_reps]]
    
    rho_mc = rho[post_reps] * rep(1,D)
    YbetaPars = MCMCresult$MCMCposteriors$YbetaLists[[post_reps]]
    Ysig2Pars = MCMCresult$MCMCposteriors$Ysig2Lists[[post_reps]]
    MbetaPars = MCMCresult$MCMCposteriors$MbetaLists[[post_reps]]
    Msig2Pars = MCMCresult$MCMCposteriors$Msig2Lists[[post_reps]]
    VbetaPars_array = MCMCresult$MCMCposteriors$VbetaLists[[post_reps]]
    Vsig2Pars_array = MCMCresult$MCMCposteriors$Vsig2Lists[[post_reps]]
    ZpiPars_array   = MCMCresult$MCMCposteriors$ZpiLists[[post_reps]]
    XpiPars_array   = MCMCresult$MCMCposteriors$XpiLists[[post_reps]]
    XmuPars_array   = MCMCresult$MCMCposteriors$XmuLists[[post_reps]]
    Xtau2Pars_array = MCMCresult$MCMCposteriors$Xtau2Lists[[post_reps]]
    alpha_theta = alpha_theta_draws[post_reps]
    alpha_omega = alpha_omega_draws[post_reps]
    
    MCresult = EDPpostmc_concon_condC(condC, q_condC, p_C1, p_C2,
                                      n_k, n_rk, Ky, Kyx, max_Kx_Sy, 
                                      D, num_MC, rho_mc,
                                      YbetaPars, Ysig2Pars,
                                      MbetaPars, Msig2Pars,
                                      VbetaPars_array, Vsig2Pars_array,
                                      ZpiPars_array, XpiPars_array, 
                                      XmuPars_array, Xtau2Pars_array,
                                      alpha_theta, alpha_omega, 
                                      f0_z0_mc, f0_z1_mc,
                                      YbetaPrior_mc, 
                                      MbetaPrior_mc, Msig2Prior_mc,
                                      VbetaPrior_mc, Vsig2Prior_mc,
                                      a_Ysig2, b_Ysig2, a_Ybeta, B_Ybeta,
                                      a_Msig2, b_Msig2, a_Mbeta, B_Mbeta,
                                      a_Vsig2, b_Vsig2, a_Vbeta, B_Vbeta,
                                      a_pi, b_pi, a_mu, b_mu, a_tau2, b_tau2)
    
    NIE[post_reps] = MCresult$NIE
    NDE[post_reps] = MCresult$NDE
    ATE[post_reps] = MCresult$ATE
    
    if (post_reps %% iter_check == 0){
      cat("Post-Processing",post_reps,"(",(post_reps/n_MCMC)*100,"%)","Time:",date(),"\n")
    }
  }
  
  # Calculate median of posterior for concentration parameters
  median_alpha_theta = median(alpha_theta_draws)
  median_alpha_omega = median(alpha_omega_draws)
  
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
  alpha_result = cbind(median_alpha_theta, median_alpha_omega)
  colnames(alpha_result) = c("alpha_theta","alpha_omega")
  rownames(alpha_result) = c("")
  
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