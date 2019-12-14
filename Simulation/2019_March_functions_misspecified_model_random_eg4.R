library(RobustGaSP)
library(rospca)
library(Rcpp)
library(RcppEigen)
library(nloptr)
library(FastGP)
library(FastGaSP)
library(rstiefel)
library(dlm)
sourceCpp(file='../src/functions.cpp') 


##several functions
neg_log_lik_shared_cov_FFBS<-function(param,kernel_type){

  G_log_det_cov=Get_G_log_det_cov(param, output, delta_x,d=1,kernel_type=kernel_type)
  G=output_2-G_log_det_cov[[1]]
  
  eigen_G=eigen(G)
  
  -(-(sum(G_log_det_cov[[2]]))*d/2-(n*k)/2*log(trace_output_2-sum(eigen_G$values[1:d]) ))
  

}

neg_log_lik_shared_cov_FFBS_fixed_A<-function(param,fixed_A,kernel_type){

  G_log_det_cov=Get_G_log_det_cov(param, output, delta_x,d=1,kernel_type=kernel_type)
  G=output_2-G_log_det_cov[[1]]

  
  -(-(sum(G_log_det_cov[[2]]))*d/2-(n*k)/2*log(trace_output_2-sum( (fixed_A%*%t(fixed_A))*G)))
  
}


pred_FFBS_FastGaSP<-function(param,A_hat,input,testing_input, output_here,d,var_data=F,kernel_type){
  beta=param[1]
  sigma_2=param[2]
  sigma_2_0=param[3]
  
  
  
  #beta=param[1:d]
  #sigma_2=param[(d+1):(2*d)]
  #sigma_2_0=param[2*d+1]
  num_testing=length(testing_input)
  num_obs=length(input)
  
  predict_all_dlm=matrix(0, num_obs,d)
  var_all_dlm=matrix(0, num_obs,d)
  output_t_A=t(output_here)%*%A_hat
  var_est=0
  
  #if(needcov==T){
  #  var_est=array(0,c(k,k,num_testing))
  #}
  
  for(i_d in 1:d){
    
    m.model=fgasp(input,output_t_A[,i_d],have_noise=TRUE,kernel_type=kernel_type)
    m.pred=predict(param=c(log(beta),log(sigma_2_0/sigma_2)),object=m.model,
                   testing_input=testing_input,var_data=FALSE,sigma_2=sigma_2)
    predict_all_dlm[,i_d]=m.pred@mean
    
    var_all_dlm[,i_d]= m.pred@var  
    
    var_est= var_est+(A_hat[,i_d]^2)%*%t(as.matrix(var_all_dlm[,i_d],num_testing,1))
  }
  
  return.list=as.list(1:2)
  return.list[[1]]=A_hat%*%t(predict_all_dlm)
  ###here I include the noise in the data
  if(var_data==F){
    return.list[[2]]=var_est
    #t(var_all_dlm)
  }else{
    return.list[[2]]=var_est+rep(sigma_2_0,k)
  }
  
  return.list
  
}


get_chol<-function(x,beta,kernel_type_real){
  
  ##kernel_type_real is from Gaussian or exponential
  
  R0_00=abs(outer(x,x,'-'))
  if(kernel_type_real=='Gauss'){
    R=pow_exp_funct(R0_00,beta,alpha_i=2)
  }else{
    R=pow_exp_funct(R0_00,beta,alpha_i=1)
  }
  
  rcppeigen_get_chol(R)
}








#############################start of simimulation####################
##beta is fixed
##try 8 cases
##sigma_0_real=sqrt(4) for Gaussian kernel and sigma_0_real=sqrt(0.25) for exponential
##d=4, k=20,  n=200, n=400, n=600

N=100
n=100
k=20
d=4

record_angle_1=matrix(0,N,4)
record_param_1=matrix(0,N,3)
record_param_pc_1=matrix(0,N,3)
record_param_QW_1_1=matrix(0,N,3)
record_param_QW_5_1=matrix(0,N,3)

record_MSE_1=matrix(0,N,4)
record_prop_diff_1=matrix(0,N,2)

beta_real=0.01


sigma_real=1
input=seq(1,n,1)


k_0=5


## fix the type of kernel
kernel_type='matern_5_2'

## have Gauss with sigma_0_real=2 or have exponential with sigma_0_real=0.25
#kernel_type_real='Gauss'
#sigma_0_real=sqrt(4)
kernel_type_real='exponential'
sigma_0_real=sqrt(0.25)

L_sample=get_chol(input,beta_real,kernel_type_real=kernel_type_real)

##Qiwei Bka set k_0=5 and their annals paper set k_0=1
for(i_angle in 1:N){
  set.seed(i_angle)
  print(i_angle)
  input=sort(input)
  delta_x=input[2:length(input)]-input[1:(length(input)-1)]
  
  ##if A is not from a Stiefel manifold then the estimation of variance has other meaning
  #A=rustiefel(k, d)  ##sample from Stiefel manifold
  ##here we use the uniform
  A=matrix(runif(k*d),k,d)
  
  Factor=matrix(0,d,n)
  
  for(i in 1: d){
    Factor[i,]=sigma_real^2*L_sample%*%rnorm(n)
  }
  
  output=A%*%Factor+matrix(rnorm(n*k,mean=0,sd=sigma_0_real),k,n)

  output_2=output%*%t(output)
  
  trace_output_2=sum(output^2)
  
  svd_output=svd(output)
  
  A_est_pc=svd_output$u[,1:d]
  #%*%diag(svd_output$d[1:d])/sqrt(n)
  
  
  ##qiwei
  output_lag_1=matrix(0,k,k)
  for(k_i in 1:k){
    for(k_j in 1:k){
      output_lag_1[k_i,k_j]=output_lag_1[k_i,k_j]+cov(output[k_i,1:(n-1)],output[k_j,(1+1):(n)])
    }
  }
  eigen_QW_1=eigen((output_lag_1)%*%t(output_lag_1))
  A_est_QW_1=eigen_QW_1$vectors[,1:d]
  
  ##qiwei
  output_lag_1_cov_5=matrix(0,k,k)
  for(k_t in 1:k_0){
    for(k_i in 1:k){
      for(k_j in 1:k){
        output_lag_1_cov_5[k_i,k_j]=output_lag_1_cov_5[k_i,k_j]+cov(output[k_i,1:(n-k_t)],output[k_j,(1+k_t):(n)])
      }
    }
  }
  eigen_QW_5=eigen((output_lag_1_cov_5)%*%t(output_lag_1_cov_5))
  A_est_QW_5=eigen_QW_5$vectors[,1:d]

  ###
  A_ini=svd_output$u[,1:d] 
  
  param_ini=c(log(.1),log(10))

  ##
  m=try(optim(param_ini,neg_log_lik_shared_cov_FFBS, kernel_type=kernel_type, method="L-BFGS-B"),silent=T)
  
  m_pc=try(optim(param_ini,neg_log_lik_shared_cov_FFBS_fixed_A, fixed_A=A_est_pc, kernel_type=kernel_type, method="L-BFGS-B"),silent=T)
  
  m_QW_1=try(optim(param_ini,neg_log_lik_shared_cov_FFBS_fixed_A,fixed_A=A_est_QW_1,kernel_type=kernel_type, method="L-BFGS-B"),silent=T)
  m_QW_5=try(optim(param_ini,neg_log_lik_shared_cov_FFBS_fixed_A, fixed_A=A_est_QW_5,kernel_type=kernel_type, method="L-BFGS-B"),silent=T)
  
  ##if they fail, find another start

  while(is.character(m)){
    param_ini=param_ini+runif(2)
    m=try(optim(param_ini,neg_log_lik_shared_cov_FFBS,kernel_type=kernel_type, method="L-BFGS-B"),silent=T)
  } 
  param_ini=c(log(.1),log(10))
  
  while(is.character(m_pc)){
    param_ini=param_ini+runif(2)
    m_pc=try(optim(param_ini,neg_log_lik_shared_cov_FFBS_fixed_A,fixed_A=A_est_pc,kernel_type=kernel_type,  method="L-BFGS-B"),silent=T)
  } 
  param_ini=c(log(.1),log(10))
  
  while(is.character(m_QW_1)){
    param_ini=param_ini+runif(2)
    m_QW_1=try(optim(param_ini,neg_log_lik_shared_cov_FFBS_fixed_A,fixed_A=A_est_QW_1, kernel_type=kernel_type, method="L-BFGS-B"),silent=T)
  } 
  param_ini=c(log(.1),log(10))
  
  while(is.character(m_QW_5)){
    param_ini=param_ini+runif(2)
    m_QW_5=try(optim(param_ini,neg_log_lik_shared_cov_FFBS_fixed_A,fixed_A=A_est_QW_5, kernel_type=kernel_type, method="L-BFGS-B"),silent=T)
  } 
  
  
  beta_hat=exp(m$par[1])
  tau_hat=exp(m$par[2])
  
  G_log_det_cov_hat=Get_G_log_det_cov( m$par, output, delta_x,d=1,kernel_type=kernel_type)

  G_hat=output_2- G_log_det_cov_hat[[1]]
  
  eigen_G_hat=eigen(G_hat)
  
  sigma_2_0_hat=(trace_output_2-sum(eigen_G_hat$values[1:d]) )/(n*k)
  
  ###pc 
  
  G_log_det_cov_pc_hat=Get_G_log_det_cov( m_pc$par, output, delta_x,d=1,kernel_type=kernel_type)
 
  G_hat_pc=output_2- G_log_det_cov_pc_hat[[1]]
  
  sigma_2_0_pc_hat=(trace_output_2-sum( (A_est_pc%*%t(A_est_pc))*G_hat_pc) )/(n*k)
  
  ##QW_1
  G_log_det_cov_QW_1_hat=Get_G_log_det_cov( m_QW_1$par, output, delta_x,d=1,kernel_type=kernel_type)
 
  G_hat_QW_1=output_2- G_log_det_cov_QW_1_hat[[1]]
  
  sigma_2_0_QW_1_hat=(trace_output_2-sum( (A_est_QW_1%*%t(A_est_QW_1))*G_hat_QW_1) )/(n*k)
  
  ##QW_5
  G_log_det_cov_QW_5_hat=Get_G_log_det_cov( m_QW_5$par, output, delta_x,d=1,kernel_type=kernel_type)

  G_hat_QW_5=output_2- G_log_det_cov_QW_5_hat[[1]]
  
  sigma_2_0_QW_5_hat=(trace_output_2-sum( (A_est_QW_5%*%t(A_est_QW_5))*G_hat_QW_5) )/(n*k)
  
  ###
  record_angle_1[i_angle,1]=rospca::angle(A,A_est_pc)
  
  #record_angle[i_angle,2]=rospca::angle(A,A_est_single)
  record_angle_1[i_angle,2]=rospca::angle(A,eigen_G_hat$vectors[,1:d])
  
  record_angle_1[i_angle,3]=rospca::angle(A,A_est_QW_1)
  
  record_angle_1[i_angle,4]=rospca::angle(A,A_est_QW_5)
  
  ###################################################################################
  ###about the estimation of the last principal angle. It is not related to the order
  #A_hat=eigen_G_hat$vectors[,1:d]
  #rospca::angle(A,cbind(A_hat[,2:d],A_hat[,1]))
  
  ##the principal angle can be computed through the SVD
  #svd_A_hat_A=svd(t(A_hat)%*%A)
  
  ##the following two are allow the same
  #min(svd_A_hat_A$d)
  #cos( rospca::angle(A,cbind(A_hat[,2:d],A_hat[,1]))*pi/2 )
  ###################################################################################
  
  record_param_1[i_angle,1]=beta_hat
  record_param_1[i_angle,2]=tau_hat*(sigma_2_0_hat)
  record_param_1[i_angle,3]=sigma_2_0_hat
  
  ##pc
  record_param_pc_1[i_angle,1]=exp(m_pc$par[1])
  record_param_pc_1[i_angle,2]=exp(m_pc$par[2])*(sigma_2_0_pc_hat)
  record_param_pc_1[i_angle,3]=sigma_2_0_pc_hat
  
  ##QW_1
  record_param_QW_1_1[i_angle,1]=exp(m_QW_1$par[1])
  record_param_QW_1_1[i_angle,2]=exp(m_QW_1$par[2])*(sigma_2_0_QW_1_hat)
  record_param_QW_1_1[i_angle,3]=sigma_2_0_QW_1_hat
  
  ##pc
  record_param_QW_5_1[i_angle,1]=exp(m_QW_5$par[1])
  record_param_QW_5_1[i_angle,2]=exp(m_QW_5$par[2])*(sigma_2_0_QW_5_hat)
  record_param_QW_5_1[i_angle,3]=sigma_2_0_QW_5_hat
  
  ###here let's compute the estimation of the mean (AZ)
  
  AZ_posterior=pred_FFBS_FastGaSP(param=record_param_1[i_angle,],A_hat=eigen_G_hat$vectors[,1:d],input=input,testing_input=input, output_here=output,d=d,kernel_type=kernel_type)
    
  Y_hat=AZ_posterior[[1]]
  Y_var=AZ_posterior[[2]]
  Y_95_lower=Y_hat+sqrt(Y_var)*qnorm(0.025)
  Y_95_upper=Y_hat+sqrt(Y_var)*qnorm(0.975)
  
  Y_hat_pc=A_est_pc%*%t(A_est_pc)%*%output
  Y_hat_QW_1=A_est_QW_1%*%t(A_est_QW_1)%*%output
  Y_hat_QW_5=A_est_QW_5%*%t(A_est_QW_5)%*%output
  
  Y_mean=A%*%Factor
  
  record_MSE_1[i_angle,1]=mean( (Y_hat_pc-Y_mean)^2)
  record_MSE_1[i_angle,2]=mean((Y_hat-Y_mean)^2)
  record_MSE_1[i_angle,3]=mean((Y_hat_QW_1-Y_mean)^2)
  record_MSE_1[i_angle,4]=mean((Y_hat_QW_5-Y_mean)^2)
  
  
  record_prop_diff_1[i_angle,1]=sum(Y_95_lower<=Y_mean & Y_95_upper>=Y_mean)/(n*k)
  record_prop_diff_1[i_angle,2]=mean(abs(Y_95_upper-Y_95_lower))
  
  #print(record_angle_1[i_angle,])
  #print(record_MSE_1[i_angle,])
  
  #print(record_prop_diff_1[i_angle,])
  
  
}


n=200

record_angle_2=matrix(0,N,4)
record_param_2=matrix(0,N,3)
record_param_pc_2=matrix(0,N,3)
record_param_QW_1_2=matrix(0,N,3)
record_param_QW_5_2=matrix(0,N,3)


record_MSE_2=matrix(0,N,4)
record_prop_diff_2=matrix(0,N,2)


##when beta get bigger this becomes another noise, then sigma_real is confounded with sigma_0_real
input=seq(1,n,1)

L_sample=get_chol(input,beta_real,kernel_type_real=kernel_type_real)

##Qiwei Bka set k_0=5 and their annals paper set k_0=1
for(i_angle in 1:N){
  set.seed(i_angle)
  print(i_angle)
  input=sort(input)
  delta_x=input[2:length(input)]-input[1:(length(input)-1)]
  
  ##if A is not from a Stiefel manifold then the estimation of variance has other meaning
  #A=rustiefel(k, d)  ##sample from Stiefel manifold
  ##here we use the uniform
  A=matrix(runif(k*d),k,d)
  
  Factor=matrix(0,d,n)
  
  for(i in 1: d){
    #set.seed(i+k)
    Factor[i,]=sigma_real^2*L_sample%*%rnorm(n)
    #  test_funct(input,beta_real[i],sigma_real[i])
    #output_2[i,]=rnorm(n)
  }
  
  output=A%*%Factor+matrix(rnorm(n*k,mean=0,sd=sigma_0_real),k,n)

  
  output_2=output%*%t(output)
  
  trace_output_2=sum(output^2)
  
  svd_output=svd(output)
  
  A_est_pc=svd_output$u[,1:d]

  
  ##qiwei
  output_lag_1=matrix(0,k,k)
  for(k_i in 1:k){
    for(k_j in 1:k){
      output_lag_1[k_i,k_j]=output_lag_1[k_i,k_j]+cov(output[k_i,1:(n-1)],output[k_j,(1+1):(n)])
    }
  }
  eigen_QW_1=eigen((output_lag_1)%*%t(output_lag_1))
  A_est_QW_1=eigen_QW_1$vectors[,1:d]
  
  ##qiwei
  output_lag_1_cov_5=matrix(0,k,k)
  for(k_t in 1:k_0){
    for(k_i in 1:k){
      for(k_j in 1:k){
        output_lag_1_cov_5[k_i,k_j]=output_lag_1_cov_5[k_i,k_j]+cov(output[k_i,1:(n-k_t)],output[k_j,(1+k_t):(n)])
      }
    }
  }
  eigen_QW_5=eigen((output_lag_1_cov_5)%*%t(output_lag_1_cov_5))
  A_est_QW_5=eigen_QW_5$vectors[,1:d]
  
  ###
  A_ini=svd_output$u[,1:d] 
  
  param_ini=c(log(.1),log(10))
  
  
  
  ##
  m=try(optim(param_ini,neg_log_lik_shared_cov_FFBS, kernel_type=kernel_type, method="L-BFGS-B"),silent=T)
  
  m_pc=try(optim(param_ini,neg_log_lik_shared_cov_FFBS_fixed_A, fixed_A=A_est_pc, kernel_type=kernel_type, method="L-BFGS-B"),silent=T)
  
  m_QW_1=try(optim(param_ini,neg_log_lik_shared_cov_FFBS_fixed_A,fixed_A=A_est_QW_1,kernel_type=kernel_type, method="L-BFGS-B"),silent=T)
  m_QW_5=try(optim(param_ini,neg_log_lik_shared_cov_FFBS_fixed_A, fixed_A=A_est_QW_5,kernel_type=kernel_type, method="L-BFGS-B"),silent=T)
  
  ##if they fail, find another start
  
  while(is.character(m)){
    param_ini=param_ini+runif(2)
    m=try(optim(param_ini,neg_log_lik_shared_cov_FFBS,kernel_type=kernel_type, method="L-BFGS-B"),silent=T)
  } 
  param_ini=c(log(.1),log(10))
  
  while(is.character(m_pc)){
    param_ini=param_ini+runif(2)
    m_pc=try(optim(param_ini,neg_log_lik_shared_cov_FFBS_fixed_A,fixed_A=A_est_pc,kernel_type=kernel_type,  method="L-BFGS-B"),silent=T)
  } 
  param_ini=c(log(.1),log(10))
  
  while(is.character(m_QW_1)){
    param_ini=param_ini+runif(2)
    m_QW_1=try(optim(param_ini,neg_log_lik_shared_cov_FFBS_fixed_A,fixed_A=A_est_QW_1, kernel_type=kernel_type, method="L-BFGS-B"),silent=T)
  } 
  param_ini=c(log(.1),log(10))
  
  while(is.character(m_QW_5)){
    param_ini=param_ini+runif(2)
    m_QW_5=try(optim(param_ini,neg_log_lik_shared_cov_FFBS_fixed_A,fixed_A=A_est_QW_5, kernel_type=kernel_type, method="L-BFGS-B"),silent=T)
  } 
  
  
  beta_hat=exp(m$par[1])
  tau_hat=exp(m$par[2])
  
  G_log_det_cov_hat=Get_G_log_det_cov( m$par, output, delta_x,d=1,kernel_type=kernel_type)
  
  G_hat=output_2- G_log_det_cov_hat[[1]]
  
  eigen_G_hat=eigen(G_hat)
  
  sigma_2_0_hat=(trace_output_2-sum(eigen_G_hat$values[1:d]) )/(n*k)
  
  ###pc 
  
  G_log_det_cov_pc_hat=Get_G_log_det_cov( m_pc$par, output, delta_x,d=1,kernel_type=kernel_type)
  
  G_hat_pc=output_2- G_log_det_cov_pc_hat[[1]]
  
  sigma_2_0_pc_hat=(trace_output_2-sum( (A_est_pc%*%t(A_est_pc))*G_hat_pc) )/(n*k)
  
  ##QW_1
  G_log_det_cov_QW_1_hat=Get_G_log_det_cov( m_QW_1$par, output, delta_x,d=1,kernel_type=kernel_type)
  
  G_hat_QW_1=output_2- G_log_det_cov_QW_1_hat[[1]]
  
  sigma_2_0_QW_1_hat=(trace_output_2-sum( (A_est_QW_1%*%t(A_est_QW_1))*G_hat_QW_1) )/(n*k)
  
  ##QW_5
  G_log_det_cov_QW_5_hat=Get_G_log_det_cov( m_QW_5$par, output, delta_x,d=1,kernel_type=kernel_type)
  
  G_hat_QW_5=output_2- G_log_det_cov_QW_5_hat[[1]]
  
  sigma_2_0_QW_5_hat=(trace_output_2-sum( (A_est_QW_5%*%t(A_est_QW_5))*G_hat_QW_5) )/(n*k)
  
  ###
  record_angle_2[i_angle,1]=rospca::angle(A,A_est_pc)
  
  #record_angle[i_angle,2]=rospca::angle(A,A_est_single)
  record_angle_2[i_angle,2]=rospca::angle(A,eigen_G_hat$vectors[,1:d])
  
  record_angle_2[i_angle,3]=rospca::angle(A,A_est_QW_1)
  
  record_angle_2[i_angle,4]=rospca::angle(A,A_est_QW_5)
  

  
  record_param_2[i_angle,1]=beta_hat
  record_param_2[i_angle,2]=tau_hat*(sigma_2_0_hat)
  record_param_2[i_angle,3]=sigma_2_0_hat
  
  ##pc
  record_param_pc_2[i_angle,1]=exp(m_pc$par[1])
  record_param_pc_2[i_angle,2]=exp(m_pc$par[2])*(sigma_2_0_pc_hat)
  record_param_pc_2[i_angle,3]=sigma_2_0_pc_hat
  
  ##QW_1
  record_param_QW_1_2[i_angle,1]=exp(m_QW_1$par[1])
  record_param_QW_1_2[i_angle,2]=exp(m_QW_1$par[2])*(sigma_2_0_QW_1_hat)
  record_param_QW_1_2[i_angle,3]=sigma_2_0_QW_1_hat
  
  ##pc
  record_param_QW_5_2[i_angle,1]=exp(m_QW_5$par[1])
  record_param_QW_5_2[i_angle,2]=exp(m_QW_5$par[2])*(sigma_2_0_QW_5_hat)
  record_param_QW_5_2[i_angle,3]=sigma_2_0_QW_5_hat
  
  
  ###here let's compute the estimation of the mean (AZ)
  
  AZ_posterior=pred_FFBS_FastGaSP(param=record_param_1[i_angle,],A_hat=eigen_G_hat$vectors[,1:d],input=input,testing_input=input, output_here=output,d=d,kernel_type=kernel_type)
  
  Y_hat=AZ_posterior[[1]]
  Y_var=AZ_posterior[[2]]
  Y_95_lower=Y_hat+sqrt(Y_var)*qnorm(0.025)
  Y_95_upper=Y_hat+sqrt(Y_var)*qnorm(0.975)
  
  Y_hat_pc=A_est_pc%*%t(A_est_pc)%*%output
  Y_hat_QW_1=A_est_QW_1%*%t(A_est_QW_1)%*%output
  Y_hat_QW_5=A_est_QW_5%*%t(A_est_QW_5)%*%output
  
  Y_mean=A%*%Factor
  
  record_MSE_2[i_angle,1]=mean( (Y_hat_pc-Y_mean)^2)
  record_MSE_2[i_angle,2]=mean((Y_hat-Y_mean)^2)
  record_MSE_2[i_angle,3]=mean((Y_hat_QW_1-Y_mean)^2)
  record_MSE_2[i_angle,4]=mean((Y_hat_QW_5-Y_mean)^2)
  
  record_prop_diff_2[i_angle,1]=sum(Y_95_lower<=Y_mean & Y_95_upper>=Y_mean)/(n*k)
  record_prop_diff_2[i_angle,2]=mean(abs(Y_95_upper-Y_95_lower))
  

  #print(record_angle_2[i_angle,])
  
  #print(record_MSE_2[i_angle,])
  
  #print(record_prop_diff_2[i_angle,])
  
}




n=400

record_angle_3=matrix(0,N,4)
record_param_3=matrix(0,N,3)
record_param_pc_3=matrix(0,N,3)
record_param_QW_1_3=matrix(0,N,3)
record_param_QW_5_3=matrix(0,N,3)


record_MSE_3=matrix(0,N,4)
record_prop_diff_3=matrix(0,N,2)


##when beta get bigger this becomes another noise, then sigma_real is confounded with sigma_0_real
input=seq(1,n,1)

L_sample=get_chol(input,beta_real,kernel_type_real=kernel_type_real)

##Qiwei Bka set k_0=5 and their annals paper set k_0=1
for(i_angle in 1:N){
  set.seed(i_angle)
  print(i_angle)
  input=sort(input)
  delta_x=input[2:length(input)]-input[1:(length(input)-1)]
  
  ##if A is not from a Stiefel manifold then the estimation of variance has other meaning
  #A=rustiefel(k, d)  ##sample from Stiefel manifold
  ##here we use the uniform
  A=matrix(runif(k*d),k,d)
  
  Factor=matrix(0,d,n)
  
  for(i in 1: d){
    #set.seed(i+k)
    Factor[i,]=sigma_real^2*L_sample%*%rnorm(n)
    #  test_funct(input,beta_real[i],sigma_real[i])
    #output_2[i,]=rnorm(n)
  }
  
  output=A%*%Factor+matrix(rnorm(n*k,mean=0,sd=sigma_0_real),k,n)
  
  
  output_2=output%*%t(output)
  
  trace_output_2=sum(output^2)
  
  svd_output=svd(output)
  
  A_est_pc=svd_output$u[,1:d]
  
  
  ##qiwei
  output_lag_1=matrix(0,k,k)
  for(k_i in 1:k){
    for(k_j in 1:k){
      output_lag_1[k_i,k_j]=output_lag_1[k_i,k_j]+cov(output[k_i,1:(n-1)],output[k_j,(1+1):(n)])
    }
  }
  eigen_QW_1=eigen((output_lag_1)%*%t(output_lag_1))
  A_est_QW_1=eigen_QW_1$vectors[,1:d]
  
  ##qiwei
  output_lag_1_cov_5=matrix(0,k,k)
  for(k_t in 1:k_0){
    for(k_i in 1:k){
      for(k_j in 1:k){
        output_lag_1_cov_5[k_i,k_j]=output_lag_1_cov_5[k_i,k_j]+cov(output[k_i,1:(n-k_t)],output[k_j,(1+k_t):(n)])
      }
    }
  }
  eigen_QW_5=eigen((output_lag_1_cov_5)%*%t(output_lag_1_cov_5))
  A_est_QW_5=eigen_QW_5$vectors[,1:d]
  
  ###
  A_ini=svd_output$u[,1:d] 
  
  param_ini=c(log(.1),log(10))
  
  
  
  ##
  m=try(optim(param_ini,neg_log_lik_shared_cov_FFBS, kernel_type=kernel_type, method="L-BFGS-B"),silent=T)
  
  m_pc=try(optim(param_ini,neg_log_lik_shared_cov_FFBS_fixed_A, fixed_A=A_est_pc, kernel_type=kernel_type, method="L-BFGS-B"),silent=T)
  
  m_QW_1=try(optim(param_ini,neg_log_lik_shared_cov_FFBS_fixed_A,fixed_A=A_est_QW_1,kernel_type=kernel_type, method="L-BFGS-B"),silent=T)
  m_QW_5=try(optim(param_ini,neg_log_lik_shared_cov_FFBS_fixed_A, fixed_A=A_est_QW_5,kernel_type=kernel_type, method="L-BFGS-B"),silent=T)
  
  ##if they fail, find another start
  
  while(is.character(m)){
    param_ini=param_ini+runif(2)
    m=try(optim(param_ini,neg_log_lik_shared_cov_FFBS,kernel_type=kernel_type, method="L-BFGS-B"),silent=T)
  } 
  param_ini=c(log(.1),log(10))
  
  while(is.character(m_pc)){
    param_ini=param_ini+runif(2)
    m_pc=try(optim(param_ini,neg_log_lik_shared_cov_FFBS_fixed_A,fixed_A=A_est_pc,kernel_type=kernel_type,  method="L-BFGS-B"),silent=T)
  } 
  param_ini=c(log(.1),log(10))
  
  while(is.character(m_QW_1)){
    param_ini=param_ini+runif(2)
    m_QW_1=try(optim(param_ini,neg_log_lik_shared_cov_FFBS_fixed_A,fixed_A=A_est_QW_1, kernel_type=kernel_type, method="L-BFGS-B"),silent=T)
  } 
  param_ini=c(log(.1),log(10))
  
  while(is.character(m_QW_5)){
    param_ini=param_ini+runif(2)
    m_QW_5=try(optim(param_ini,neg_log_lik_shared_cov_FFBS_fixed_A,fixed_A=A_est_QW_5, kernel_type=kernel_type, method="L-BFGS-B"),silent=T)
  } 
  
  
  beta_hat=exp(m$par[1])
  tau_hat=exp(m$par[2])
  
  G_log_det_cov_hat=Get_G_log_det_cov( m$par, output, delta_x,d=1,kernel_type=kernel_type)
  
  G_hat=output_2- G_log_det_cov_hat[[1]]
  
  eigen_G_hat=eigen(G_hat)
  
  sigma_2_0_hat=(trace_output_2-sum(eigen_G_hat$values[1:d]) )/(n*k)
  
  ###pc 
  
  G_log_det_cov_pc_hat=Get_G_log_det_cov( m_pc$par, output, delta_x,d=1,kernel_type=kernel_type)
  
  G_hat_pc=output_2- G_log_det_cov_pc_hat[[1]]
  
  sigma_2_0_pc_hat=(trace_output_2-sum( (A_est_pc%*%t(A_est_pc))*G_hat_pc) )/(n*k)
  
  ##QW_1
  G_log_det_cov_QW_1_hat=Get_G_log_det_cov( m_QW_1$par, output, delta_x,d=1,kernel_type=kernel_type)
  
  G_hat_QW_1=output_2- G_log_det_cov_QW_1_hat[[1]]
  
  sigma_2_0_QW_1_hat=(trace_output_2-sum( (A_est_QW_1%*%t(A_est_QW_1))*G_hat_QW_1) )/(n*k)
  
  ##QW_5
  G_log_det_cov_QW_5_hat=Get_G_log_det_cov( m_QW_5$par, output, delta_x,d=1,kernel_type=kernel_type)
  
  G_hat_QW_5=output_2- G_log_det_cov_QW_5_hat[[1]]
  
  sigma_2_0_QW_5_hat=(trace_output_2-sum( (A_est_QW_5%*%t(A_est_QW_5))*G_hat_QW_5) )/(n*k)
  
  ###
  record_angle_3[i_angle,1]=rospca::angle(A,A_est_pc)
  
  #record_angle[i_angle,2]=rospca::angle(A,A_est_single)
  record_angle_3[i_angle,2]=rospca::angle(A,eigen_G_hat$vectors[,1:d])
  
  record_angle_3[i_angle,3]=rospca::angle(A,A_est_QW_1)
  
  record_angle_3[i_angle,4]=rospca::angle(A,A_est_QW_5)
  
  
  
  record_param_3[i_angle,1]=beta_hat
  record_param_3[i_angle,2]=tau_hat*(sigma_2_0_hat)
  record_param_3[i_angle,3]=sigma_2_0_hat
  
  ##pc
  record_param_pc_3[i_angle,1]=exp(m_pc$par[1])
  record_param_pc_3[i_angle,2]=exp(m_pc$par[2])*(sigma_2_0_pc_hat)
  record_param_pc_3[i_angle,3]=sigma_2_0_pc_hat
  
  ##QW_1
  record_param_QW_1_3[i_angle,1]=exp(m_QW_1$par[1])
  record_param_QW_1_3[i_angle,2]=exp(m_QW_1$par[2])*(sigma_2_0_QW_1_hat)
  record_param_QW_1_3[i_angle,3]=sigma_2_0_QW_1_hat
  
  ##pc
  record_param_QW_5_3[i_angle,1]=exp(m_QW_5$par[1])
  record_param_QW_5_3[i_angle,2]=exp(m_QW_5$par[2])*(sigma_2_0_QW_5_hat)
  record_param_QW_5_3[i_angle,3]=sigma_2_0_QW_5_hat
  
  
  ###here let's compute the estimation of the mean (AZ)
  
  AZ_posterior=pred_FFBS_FastGaSP(param=record_param_1[i_angle,],A_hat=eigen_G_hat$vectors[,1:d],input=input,testing_input=input, output_here=output,d=d,kernel_type=kernel_type)
  
  Y_hat=AZ_posterior[[1]]
  Y_var=AZ_posterior[[2]]
  Y_95_lower=Y_hat+sqrt(Y_var)*qnorm(0.025)
  Y_95_upper=Y_hat+sqrt(Y_var)*qnorm(0.975)
  
  Y_hat_pc=A_est_pc%*%t(A_est_pc)%*%output
  Y_hat_QW_1=A_est_QW_1%*%t(A_est_QW_1)%*%output
  Y_hat_QW_5=A_est_QW_5%*%t(A_est_QW_5)%*%output
  
  Y_mean=A%*%Factor
  
  record_MSE_3[i_angle,1]=mean( (Y_hat_pc-Y_mean)^2)
  record_MSE_3[i_angle,2]=mean((Y_hat-Y_mean)^2)
  record_MSE_3[i_angle,3]=mean((Y_hat_QW_1-Y_mean)^2)
  record_MSE_3[i_angle,4]=mean((Y_hat_QW_5-Y_mean)^2)
  
  record_prop_diff_3[i_angle,1]=sum(Y_95_lower<=Y_mean & Y_95_upper>=Y_mean)/(n*k)
  record_prop_diff_3[i_angle,2]=mean(abs(Y_95_upper-Y_95_lower))
  
  
  #print(record_angle_3[i_angle,])
  
  #print(record_MSE_3[i_angle,])
  
  #print(record_prop_diff_3[i_angle,])
  
}




n
k
d

colMeans(record_MSE_1)
colMeans(record_prop_diff_1)

colMeans(record_MSE_2)
colMeans(record_prop_diff_2)

colMeans(record_MSE_3)
colMeans(record_prop_diff_3)

##this is from pi/2
record_angle=cbind(record_angle_1,record_angle_2,record_angle_3)*pi/2

pdf("simulation_shared_cov_exponential_kernel_misspecified.pdf",height=5,width=6)
boxplot(record_angle,las=2,names=c('PCA','GPPCA','LY1','LY5','PCA','GPPCA','LY1','LY5','PCA','GPPCA','LY1','LY5'),
        col=c('red','blue','yellow','green','red','blue','yellow','green','red','blue','yellow','green'),
        ylab='Largest Principal Angle',ylim=c(0,pi/2),at=c(1,2,3,4,6,7,8,9,11,12,13,14))
Lines <- list(bquote("Exponential Kernel"),bquote(list(k==20,~d==4~"and"~tau==0.25)))
mtext(do.call(expression, Lines),side=3,line=1:0)
dev.off()




# ####record of the results
# ##exponential kernel
# > colMeans(record_MSE_1)
# [1] 0.07407344 0.03117674 0.14596418 0.13234743
# > colMeans(record_prop_diff_1)
# [1] 0.6794850 0.3575094
# > 
#   > colMeans(record_MSE_2)
# [1] 0.06094725 0.02631910 0.08230548 0.07336631
# > colMeans(record_prop_diff_2)
# [1] 0.7217575 0.3551905
# > 
#   > colMeans(record_MSE_3)
# [1] 0.05404183 0.02376615 0.05704994 0.05562626
# > colMeans(record_prop_diff_3)
# [1] 0.7462237 0.3536855


# # ##Gaussian kernel
# > colMeans(record_MSE_1)
# [1] 1.060543 0.718913 1.309999 1.302126
# > colMeans(record_prop_diff_1)
# [1] 0.729145 1.882897
# > 
#   > colMeans(record_MSE_2)
# [1] 0.8949066 0.6568222 0.9998848 0.9923675
# > colMeans(record_prop_diff_2)
# [1] 0.751760 1.876184
# > 
#   > colMeans(record_MSE_3)
# [1] 0.8406704 0.6194822 0.8633699 0.8625833
# > colMeans(record_prop_diff_3)
# [1] 0.76579 1.87015
 
