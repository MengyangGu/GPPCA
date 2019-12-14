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

##for Ind GP
neg_log_lik_shared_cov_FFBS_single<-function(param,output_here,kernel_type){

  G_log_det_cov=Get_G_log_det_cov(param, output_here, delta_x,d=1,kernel_type=kernel_type)

  -(-(sum(G_log_det_cov[[2]]))*1/2-n/2*log(diag(G_log_det_cov[[1]])))
  
}

##for PP GP
neg_log_lik_shared_cov_FFBS_all<-function(param,output_here,kernel_type){
  G_log_det_cov=Get_G_log_det_cov(param, output_here, delta_x,d=1,kernel_type=kernel_type)

  -(-(sum(G_log_det_cov[[2]]))*k/2-(n*k)/2*log(sum(diag(G_log_det_cov[[1]]))))
  

}


##I can re-implement it using FastGaSP
pred_mean_independent_GP_FFBS_funct<-function(param_record,delta_x, output_here){
  beta=exp(param_record[,1])
  sigma_2=exp(param_record[,2])
  #sigma_2_0=param[3]
  
  #tau=sigma_2/sigma_2_0
  
  lambda_beta=sqrt(2*5/2)*beta # this is from beta
  VV=1
  #sigma_2_0
  
  predict_all_dlm=matrix(0, n,k)
  #var_all_dlm=matrix(0, n,d)
  
  #output_t_A=t(output_here)%*%A_hat
  FF=matrix(0,1,3) ###3 \times 1 matrix
  FF[1,1]=1   ####fixed
  
  
  for(i_k in 1:k){
    
    Q0=Construct_Q0(sigma_2[i_k],lambda_beta[i_k])
    Qi=Construct_Qi(sigma_2[i_k],delta_x,lambda_beta[i_k],Q0)
    GG=Construct_GG(delta_x,lambda_beta[i_k])
    
    time_varying_comp=cbind(GG[[1]],Qi[[1]] ) 
    
    dlm_matern_5_2=dlm(FF= FF, V=VV[1], GG=diag(3),W=diag(3), m0=c(0,0,0),C0=Q0[[1]], JGG=matrix(1:9, nr = 3), JW=matrix(10:18,nr=3),X=time_varying_comp )
    
    dlm_filter_matern_5_2 <- dlmFilter(output_here[i_k,], dlm_matern_5_2)###filter
    
    dlm_smooth_matern_5_2 <- dlmSmooth(dlm_filter_matern_5_2)
    
    predict_all_dlm[,i_k]=dropFirst(dlm_smooth_matern_5_2$s[,1])
  } 
  ##only
  predict_all_dlm
  
}





#############################start of simimulation####################
##beta is fixed
##try 4 cases
##sigma_0_real=sqrt(0.25) 
##d=4, k=20, n=100, n=200, n=400, n=800

N=100
n=100
k=20
d=4

record_angle_1=matrix(0,N,4)
record_param_1=matrix(0,N,3)
record_param_pc_1=matrix(0,N,3)
record_param_QW_1_1=matrix(0,N,3)
record_param_QW_5_1=matrix(0,N,3)

record_MSE_1=matrix(0,N,6)
record_prop_diff_1=matrix(0,N,2)
param_record_est=matrix(0,k,2)

sigma_0_real=sqrt(0.25)
sigma_real=1

input=seq(1,n,1)


k_0=5


## fix the type of kernel
kernel_type='matern_5_2'

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
    theta=runif(1)
    Factor[i,]=sigma_real^2*cos(.05*pi*theta*input)
    
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
  
  
  
  
  ###Ind GP and PP GP
  
  for(index_k in 1:k){
    #print(index_k)
    output_matrix=t(as.matrix(output[index_k,]))
    model=try(optim(param_ini,neg_log_lik_shared_cov_FFBS_single,output_here=output_matrix, kernel_type=kernel_type, method="L-BFGS-B"),silent=T)
    param_record_est[index_k,]=model$par
  }
  
  
  
  ## Ind GP and PP GP
  pred_ind_GP=pred_mean_independent_GP_FFBS_funct(param_record_est,delta_x,output)
  

  record_MSE_1[i_angle,5]=mean((t(pred_ind_GP)-Y_mean)^2)
  
  
  model_all=try(optim(param_ini,neg_log_lik_shared_cov_FFBS_all,output_here=output, kernel_type=kernel_type,method="L-BFGS-B"),silent=T)
  
  pred_PP_GP=pred_mean_independent_GP_FFBS_funct( t(matrix(model_all$par,2,k)),delta_x,output)
  
  record_MSE_1[i_angle,6]=mean((t(pred_PP_GP)-Y_mean)^2)
  
  

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


record_MSE_2=matrix(0,N,6)
record_prop_diff_2=matrix(0,N,2)


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
    theta=runif(1)
    Factor[i,]=sigma_real^2*cos(.05*pi*theta*input)
    
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
  

  
  
  ###Ind GP and PP GP
  
  for(index_k in 1:k){
    #print(index_k)
    output_matrix=t(as.matrix(output[index_k,]))
    model=try(optim(param_ini,neg_log_lik_shared_cov_FFBS_single,output_here=output_matrix, kernel_type=kernel_type, method="L-BFGS-B"),silent=T)
    param_record_est[index_k,]=model$par
  }
  
  
  
  ## Ind GP and PP GP
  pred_ind_GP=pred_mean_independent_GP_FFBS_funct(param_record_est,delta_x,output)
  
  
  record_MSE_2[i_angle,5]=mean((t(pred_ind_GP)-Y_mean)^2)
  
  
  model_all=try(optim(param_ini,neg_log_lik_shared_cov_FFBS_all,output_here=output, kernel_type=kernel_type,method="L-BFGS-B"),silent=T)
  
  pred_PP_GP=pred_mean_independent_GP_FFBS_funct( t(matrix(model_all$par,2,k)),delta_x,output)
  
  record_MSE_2[i_angle,6]=mean((t(pred_PP_GP)-Y_mean)^2)
  
  
  
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


record_MSE_3=matrix(0,N,6)
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
    theta=runif(1)
    Factor[i,]=sigma_real^2*cos(.05*pi*theta*input)
    
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
  
  
  ###Ind GP and PP GP
  
  for(index_k in 1:k){
    #print(index_k)
    output_matrix=t(as.matrix(output[index_k,]))
    model=try(optim(param_ini,neg_log_lik_shared_cov_FFBS_single,output_here=output_matrix, kernel_type=kernel_type, method="L-BFGS-B"),silent=T)
    param_record_est[index_k,]=model$par
  }
  
  
  
  ## Ind GP and PP GP
  pred_ind_GP=pred_mean_independent_GP_FFBS_funct(param_record_est,delta_x,output)
  
  
  record_MSE_3[i_angle,5]=mean((t(pred_ind_GP)-Y_mean)^2)
  
  
  model_all=try(optim(param_ini,neg_log_lik_shared_cov_FFBS_all,output_here=output, kernel_type=kernel_type,method="L-BFGS-B"),silent=T)
  
  pred_PP_GP=pred_mean_independent_GP_FFBS_funct( t(matrix(model_all$par,2,k)),delta_x,output)
  
  record_MSE_3[i_angle,6]=mean((t(pred_PP_GP)-Y_mean)^2)
  
  
  #print(record_angle_3[i_angle,])
  
  #print(record_MSE_3[i_angle,])
  
  #print(record_prop_diff_3[i_angle,])
  
}



n=800

record_angle_4=matrix(0,N,4)
record_param_4=matrix(0,N,3)
record_param_pc_4=matrix(0,N,3)
record_param_QW_1_4=matrix(0,N,3)
record_param_QW_5_4=matrix(0,N,3)


record_MSE_4=matrix(0,N,6)
record_prop_diff_4=matrix(0,N,2)


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
    theta=runif(1)
    Factor[i,]=sigma_real^2*cos(.05*pi*theta*input)
    
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
  record_angle_4[i_angle,1]=rospca::angle(A,A_est_pc)
  
  #record_angle[i_angle,2]=rospca::angle(A,A_est_single)
  record_angle_4[i_angle,2]=rospca::angle(A,eigen_G_hat$vectors[,1:d])
  
  record_angle_4[i_angle,3]=rospca::angle(A,A_est_QW_1)
  
  record_angle_4[i_angle,4]=rospca::angle(A,A_est_QW_5)
  
  
  
  record_param_4[i_angle,1]=beta_hat
  record_param_4[i_angle,2]=tau_hat*(sigma_2_0_hat)
  record_param_4[i_angle,3]=sigma_2_0_hat
  
  ##pc
  record_param_pc_4[i_angle,1]=exp(m_pc$par[1])
  record_param_pc_4[i_angle,2]=exp(m_pc$par[2])*(sigma_2_0_pc_hat)
  record_param_pc_4[i_angle,3]=sigma_2_0_pc_hat
  
  ##QW_1
  record_param_QW_1_4[i_angle,1]=exp(m_QW_1$par[1])
  record_param_QW_1_4[i_angle,2]=exp(m_QW_1$par[2])*(sigma_2_0_QW_1_hat)
  record_param_QW_1_4[i_angle,3]=sigma_2_0_QW_1_hat
  
  ##pc
  record_param_QW_5_4[i_angle,1]=exp(m_QW_5$par[1])
  record_param_QW_5_4[i_angle,2]=exp(m_QW_5$par[2])*(sigma_2_0_QW_5_hat)
  record_param_QW_5_4[i_angle,3]=sigma_2_0_QW_5_hat
  
  
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
  
  record_MSE_4[i_angle,1]=mean( (Y_hat_pc-Y_mean)^2)
  record_MSE_4[i_angle,2]=mean((Y_hat-Y_mean)^2)
  record_MSE_4[i_angle,3]=mean((Y_hat_QW_1-Y_mean)^2)
  record_MSE_4[i_angle,4]=mean((Y_hat_QW_5-Y_mean)^2)
  
  record_prop_diff_4[i_angle,1]=sum(Y_95_lower<=Y_mean & Y_95_upper>=Y_mean)/(n*k)
  record_prop_diff_4[i_angle,2]=mean(abs(Y_95_upper-Y_95_lower))
  
  
  ###Ind GP and PP GP
  
  for(index_k in 1:k){
    #print(index_k)
    output_matrix=t(as.matrix(output[index_k,]))
    model=try(optim(param_ini,neg_log_lik_shared_cov_FFBS_single,output_here=output_matrix, kernel_type=kernel_type, method="L-BFGS-B"),silent=T)
    param_record_est[index_k,]=model$par
  }
  
  
  
  ## Ind GP and PP GP
  pred_ind_GP=pred_mean_independent_GP_FFBS_funct(param_record_est,delta_x,output)
  
  
  record_MSE_4[i_angle,5]=mean((t(pred_ind_GP)-Y_mean)^2)
  
  
  model_all=try(optim(param_ini,neg_log_lik_shared_cov_FFBS_all,output_here=output, kernel_type=kernel_type,method="L-BFGS-B"),silent=T)
  
  pred_PP_GP=pred_mean_independent_GP_FFBS_funct( t(matrix(model_all$par,2,k)),delta_x,output)
  
  record_MSE_4[i_angle,6]=mean((t(pred_PP_GP)-Y_mean)^2)
  
  
  #print(record_angle_4[i_angle,])
  
  #print(record_MSE_4[i_angle,])
  
  #print(record_prop_diff_4[i_angle,])
  
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

colMeans(record_MSE_4)
colMeans(record_prop_diff_4)


##this is from pi/2
record_angle=cbind(record_angle_1,record_angle_2,record_angle_3,record_angle_4)*pi/2

pdf("simulation_shared_cov_misspecified_deterministic.pdf",height=5,width=10)
boxplot(record_angle,las=2,names=c('PCA','GPPCA','LY1','LY5','PCA','GPPCA','LY1','LY5','PCA','GPPCA','LY1','LY5',
                                   'PCA','GPPCA','LY1','LY5'),
        col=c('red','blue','yellow','green','red','blue','yellow','green','red','blue','yellow','green',
              'red','blue','yellow','green'),
        ylab='Largest Principal Angle',ylim=c(0,pi/2),at=c(1,2,3,4,6,7,8,9,11,12,13,14,16,17,18,19))
Lines <- list(bquote("Deterministic Factors"),bquote(list(k==20,~d==4~"and"~tau==4)))
mtext(do.call(expression, Lines),side=3,line=1:0)
dev.off()





# ####record of the results
# > colMeans(record_MSE_1)
# [1] 0.07033202 0.01400260 0.09777058 0.09323400 0.02057671 0.01978409
# > colMeans(record_prop_diff_1)
# [1] 0.7559750 0.2807816
# > 
#   > colMeans(record_MSE_2)
# [1] 0.059652454 0.009196822 0.075762631 0.072799074 0.018531694 0.018516994
# > colMeans(record_prop_diff_2)
# [1] 0.8369775 0.2757958
# > 
#   > colMeans(record_MSE_3)
# [1] 0.054334915 0.006732468 0.062745591 0.061839852 0.017403122 0.017812079
# > colMeans(record_prop_diff_3)
# [1] 0.8946013 0.2735651
# > 
#   > colMeans(record_MSE_4)
# [1] 0.052196165 0.005516745 0.056758797 0.055973770 0.016958573 0.017527411
# > colMeans(record_prop_diff_4)
# [1] 0.9292662 0.2723262

# ##Gaussian kernel
