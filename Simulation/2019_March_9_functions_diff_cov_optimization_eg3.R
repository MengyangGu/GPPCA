library(RobustGaSP)
library(FastGaSP)
library(FastGP)
library(Rcpp)
library(RcppEigen)
library(rospca)
library(rstiefel)

sourceCpp(file='../src/functions.cpp') 

neg_log_lik_diff_cov_FFBS<-function(param,A_ini,kernel_type='matern_5_2'){
  
  
  G_log_det_cov=Get_G_log_det_cov(param,output, delta_x,d,kernel_type)
  # G_log_det_cov=Get_G_log_det_cov(c(rep(param[1],d),param[-1]), output_sub_mean, delta_x,d)
  
  G=as.list(1:d)
  
  for(i in 1:d){
    G[[i]]=output_2-G_log_det_cov[[i]]
  }
  
  ##my method
  A_hat_here=Optimization_Stiefel_Manifold(A_ini, G=G,max_iter=200)
  
  S_2=trace_output_2+F_Funct(A_hat_here,G)
  neg_log_lik=-(-(sum(G_log_det_cov[[d+1]]))/2-(n*k)/2*log(S_2))
  # print(param)
  # print(S_2/(n*k))
  # print(neg_log_lik)
  return(neg_log_lik)
  
  #-(-(sum(G_log_det_cov[[d+1]]))/2-(n*k)/2*log(trace_output_sub_mean_2+F_Funct(A_hat_here,G)))
  
}

neg_log_lik_diff_cov_FFBS_fixed_A<-function(param,fixed_A,kernel_type='matern_5_2'){
  #print(param)
  
  
  G_log_det_cov=Get_G_log_det_cov(param, output, delta_x,d,kernel_type)
  G=as.list(1:d)
  
  for(i in 1:d){
    G[[i]]=output_2-G_log_det_cov[[i]]
  }
  
  -(-(sum(G_log_det_cov[[d+1]]))/2-(n*k)/2*log(trace_output_2+F_Funct(fixed_A,G)))
  
  
}
pred_FFBS_FastGaSP<-function(param,A_hat,input,testing_input, output_here,d,var_data=F,kernel_type='matern_5_2'){
  
  
  beta=param[1:d]
  sigma_2=param[(d+1):(2*d)]
  sigma_2_0=param[2*d+1]
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
    m.pred=predict(param=c(log(beta[i_d]),log(sigma_2_0/sigma_2[i_d])),object=m.model,
                   testing_input=testing_input,var_data=FALSE,sigma_2=sigma_2[i_d])
    predict_all_dlm[,i_d]=m.pred@mean
    
    var_all_dlm[,i_d]= m.pred@var  
    
    var_est= var_est+(A_hat[,i_d]^2)%*%t(as.matrix(var_all_dlm[,i_d],num_testing,1))
    
    ###this line seem to calculate only the var, but now I need the whole cov...
    #if(needcov==F){
    #   var_est= var_est+(A_hat[,i_d]^2)%*%t(as.matrix(var_all_dlm[,i_d],num_testing,1))
    #}else{
    #  var_est
    #}
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
  
  ##rowSums(A_hat^2)%*%t(as.matrix(var_all_dlm,n,1))
  
  return.list
  
}


N=100
n=200
k=8
d=4

record_angle_1=matrix(0,N,4)
record_param_1=matrix(0,N,2*d+1)
record_param_pc_1=matrix(0,N,2*d+1)
record_param_QW_1_1=matrix(0,N,2*d+1)
record_param_QW_5_1=matrix(0,N,2*d+1)


#beta_real=rep(0.01,d)
sigma_real=rep(1,d)
#rep(1,d)
#sigma_0_real=.2
#sigma_0_real=sqrt(.01)
sigma_0_real=sqrt(.25)


G=list()
input=seq(1,n,1)
R0_00=abs(outer(input,input,'-'))


k_0=5
G_hat=as.list(1:d)
G_hat_pc=as.list(1:d)
G_hat_QW_1=as.list(1:d)
G_hat_QW_5=as.list(1:d)


record_MSE_1=matrix(0,N,4)
record_prop_diff_1=matrix(0,N,2)
record_prop_diff_data_1=matrix(0,N,2)


for(i_angle in 1:N){
  set.seed(i_angle)
  print(i_angle)
  #set.seed(2)
  #input=seq(0,1,1/(n-1))
  ##
  #min_beta=0.02
  #max_beta=0.005
  beta_real=runif(d,min=0.001,max=.1)
  
  L_sample_list=as.list(1:d)
  
  for(i in 1:d){
    # L_sample_list[[i]]=get_chol(input,beta_real[i])
    # 
    # L_sample_list2=Chol_rcppeigen(R_real)
    # L_sample_list2- t(chol(R_real))
    R_real=matern_5_2_funct(R0_00,beta_real[i])
    L_sample_list[[i]]=Chol_rcppeigen(R_real)
    
  }
  
  
  
  input=sort(input)
  delta_x=input[2:length(input)]-input[1:(length(input)-1)]
  
  #A=matrix(runif(k*d),k,d)
  
  ##if A is not from a Stiefel manifold then the estimation of variance could be wrong
  A=rustiefel(k, d)  ##sample from Stiefel manifold
  
  Factor=matrix(0,d,n)
  
  for(i_fac in 1: d){
    Factor[i_fac,]=sigma_real[i_fac]^2*L_sample_list[[i_fac]]%*%rnorm(n)
  }
  
  output_ori=A%*%Factor+matrix(rnorm(n*k,mean=0,sd=sigma_0_real),k,n)
  
  output=output_ori

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
  
  
  A_ini=svd_output$u[,1:d] 
  
  param_ini=c(rep(log(.1),d),rep(log(10),d))
  
  
  #sink(tempfile())
  #options(warn = -1) 
  m=try(optim(param_ini,neg_log_lik_diff_cov_FFBS,A_ini=A_ini,lower=c(rep(log(10^{-10}),d),rep(log(10^{-10}),d)),
              upper=c(rep(log(10^10),d),rep(log(10^10),d)),control = list(maxit = 100), method="L-BFGS-B"),silent=T)
  
  
  
  m_pc=try(optim(param_ini,neg_log_lik_diff_cov_FFBS_fixed_A, fixed_A=A_est_pc,lower=c(rep(log(10^{-10}),d),rep(log(10^{-10}),d)),
                 upper=c(rep(log(10^10),d),rep(log(10^10),d)),control = list(maxit = 100), method="L-BFGS-B"),silent=T)
  
  m_QW_1=try(optim(param_ini,neg_log_lik_diff_cov_FFBS_fixed_A,fixed_A=A_est_QW_1,lower=c(rep(log(10^{-10}),d),rep(log(10^{-10}),d)),
                   upper=c(rep(log(10^10),d),rep(log(10^10),d)),control = list(maxit = 100), method="L-BFGS-B"),silent=T)
  m_QW_5=try(optim(param_ini,neg_log_lik_diff_cov_FFBS_fixed_A, fixed_A=A_est_QW_5,lower=c(rep(log(10^{-10}),d),rep(log(10^{-10}),d)),
                   upper=c(rep(log(10^10),d),rep(log(10^10),d)),control = list(maxit = 100), method="L-BFGS-B"),silent=T)
  
  
  while(is.character(m)){
    param_ini=param_ini+runif(2*d)
    m=try(optim(param_ini,neg_log_lik_diff_cov_FFBS,A_ini=A_ini,lower=c(rep(log(10^{-10}),d),rep(log(10^{-10}),d)),
                upper=c(rep(log(10^10),d),rep(log(10^10),d)),control = list(maxit = 100), method="L-BFGS-B"),silent=T)
  } 
  
  while(is.character(m_pc)){
    param_ini=param_ini+runif(2*d)
    m_pc=try(optim(param_ini,neg_log_lik_diff_cov_FFBS_fixed_A, fixed_A=A_est_pc,lower=c(rep(log(10^{-10}),d),rep(log(10^{-10}),d)),
                   upper=c(rep(log(10^10),d),rep(log(10^10),d)),control = list(maxit = 100), method="L-BFGS-B"),silent=T)
  } 
  
  while(is.character(m_QW_1)){
    param_ini=param_ini+runif(2*d)
    m_QW_1=try(optim(param_ini,neg_log_lik_diff_cov_FFBS_fixed_A,fixed_A=A_est_QW_1,lower=c(rep(log(10^{-10}),d),rep(log(10^{-10}),d)),
                     upper=c(rep(log(10^10),d),rep(log(10^10),d)),control = list(maxit = 100), method="L-BFGS-B"),silent=T)
  } 
  
  while(is.character(m_QW_5)){
    param_ini=param_ini+runif(2*d)
    m_QW_5=try(optim(param_ini,neg_log_lik_diff_cov_FFBS_fixed_A, fixed_A=A_est_QW_5,
                     upper=c(rep(log(10^10),d),rep(log(10^10),d)),lower=c(rep(log(10^{-10}),d),rep(log(10^{-10}),d)),control = list(maxit = 100), method="L-BFGS-B"),silent=T)
  } 
  
  #m=try(optim(param_ini,neg_log_lik_diff_cov_FFBS,A_ini=A_ini,lower=c(rep(log(10^{-10}),d),rep(log(10^{-10}),d))),silent=T)
  
  #sink()
  
  
  ##my method
  beta_hat=exp(m$par[1:d])
  tau_hat=exp(m$par[(d+1):(2*d)])
  
  
  
  G_log_det_cov_hat=Get_G_log_det_cov(m$par, output, delta_x,d,kernel_type = 'matern_5_2')
  
  for(i in 1:d){
    G_hat[[i]]=output_2-G_log_det_cov_hat[[i]]
  }
  
  A_hat=Optimization_Stiefel_Manifold(A_ini, G=G_hat,max_iter=100)
  
  
  sigma_2_0_hat=(trace_output_2+F_Funct(A_hat,G_hat))/(n*k)
  
  ###pc 
  
  G_log_det_cov_pc_hat=Get_G_log_det_cov( m_pc$par, output, delta_x,d,kernel_type = 'matern_5_2')
  
  for(i in 1:d){
    G_hat_pc[[i]]=output_2-G_log_det_cov_pc_hat[[i]]
  }
  
  sigma_2_0_pc_hat=(trace_output_2+F_Funct(A_est_pc,G_hat_pc) )/(n*k)
  
  ##QW_1
  G_log_det_cov_QW_1_hat=Get_G_log_det_cov( m_QW_1$par, output, delta_x,d,kernel_type = 'matern_5_2')
  
  for(i in 1:d){
    G_hat_QW_1[[i]]=output_2-G_log_det_cov_QW_1_hat[[i]]
  }
  
  sigma_2_0_QW_1_hat=(trace_output_2+F_Funct(A_est_QW_1,G_hat_QW_1) )/(n*k)
  
  ##QW_5
  G_log_det_cov_QW_5_hat=Get_G_log_det_cov( m_QW_5$par, output, delta_x,d,kernel_type = 'matern_5_2')
  
  for(i in 1:d){
    G_hat_QW_5[[i]]=output_2-G_log_det_cov_QW_5_hat[[i]]
  }
  
  sigma_2_0_QW_5_hat=(trace_output_2+F_Funct(A_est_QW_5,G_hat_QW_5) )/(n*k)
  
  
  record_angle_1[i_angle,1]=rospca::angle(A,A_est_pc)
  
  record_angle_1[i_angle,2]=rospca::angle(A,A_hat)
  
  
  record_angle_1[i_angle,3]=rospca::angle(A,A_est_QW_1)
  
  record_angle_1[i_angle,4]=rospca::angle(A,A_est_QW_5)
  
  
  
  record_param_1[i_angle,1:d]=beta_hat
  record_param_1[i_angle,(d+1):(2*d)]=tau_hat*(sigma_2_0_hat)
  record_param_1[i_angle,2*d+1]=sigma_2_0_hat
  
  ##pc
  record_param_pc_1[i_angle,1:d]=exp(m_pc$par[1:d])
  record_param_pc_1[i_angle,(d+1):(2*d)]=exp(m_pc$par[(d+1):(2*d) ])*(sigma_2_0_pc_hat)
  record_param_pc_1[i_angle,2*d+1]=sigma_2_0_pc_hat
  
  ##QW_1
  record_param_QW_1_1[i_angle,1:d]=exp(m_QW_1$par[1:d])
  record_param_QW_1_1[i_angle,(d+1):(2*d)]=exp(m_QW_1$par[(d+1):(2*d)])*as.numeric(sigma_2_0_QW_1_hat)
  record_param_QW_1_1[i_angle,2*d+1]=sigma_2_0_QW_1_hat
  
  ##QW_5
  record_param_QW_5_1[i_angle,1:d]=exp(m_QW_5$par[1:d])
  record_param_QW_5_1[i_angle,(d+1):(2*d)]=exp(m_QW_5$par[(d+1):(2*d)])*as.numeric(sigma_2_0_QW_5_hat)
  record_param_QW_5_1[i_angle,2*d+1]=sigma_2_0_QW_5_hat
  
  
  
  print(record_angle_1[i_angle,]) 
  
  ##prediction
  Y_hat_pc=A_est_pc%*%t(A_est_pc)%*%output
  Y_hat_QW_1=A_est_QW_1%*%t(A_est_QW_1)%*%output
  Y_hat_QW_5=A_est_QW_5%*%t(A_est_QW_5)%*%output
  
  pred_GPPCA=pred_FFBS_FastGaSP(param=record_param_1[i_angle,],A_hat=A_hat,input=input,testing_input=input, output_here=output,d=d,var_data=F)
    
  Y_hat=pred_GPPCA[[1]]
  Y_var=pred_GPPCA[[2]]
  Y_95_lower=Y_hat+sqrt(Y_var)*qnorm(0.025)
  Y_95_upper=Y_hat+sqrt(Y_var)*qnorm(0.975)
  
  
  Y_mean=A%*%Factor
  
  record_MSE_1[i_angle,1]=mean( (Y_hat_pc-Y_mean)^2)
  record_MSE_1[i_angle,2]=mean((Y_hat-Y_mean)^2)
  record_MSE_1[i_angle,3]=mean((Y_hat_QW_1-Y_mean)^2)
  record_MSE_1[i_angle,4]=mean((Y_hat_QW_5-Y_mean)^2)
  
  
  record_prop_diff_1[i_angle,1]=sum(Y_95_lower<=Y_mean & Y_95_upper>=Y_mean)/(n*k)
  record_prop_diff_1[i_angle,2]=mean(abs(Y_95_upper-Y_95_lower))
  
  Y_var_data=Y_var+rep(sigma_2_0_hat,k)
  Y_95_data_lower=Y_hat+sqrt(Y_var_data)*qnorm(0.025)
  Y_95_data_upper=Y_hat+sqrt(Y_var_data)*qnorm(0.975)
  
  record_prop_diff_data_1[i_angle,1]=sum(Y_95_data_lower<=output & Y_95_data_upper>=output)/(n*k)
  record_prop_diff_data_1[i_angle,2]=mean(abs(Y_95_data_upper-Y_95_data_lower))
  

  print(record_MSE_1[i_angle,])
  
}





par(mfrow=c(1,2))
index_d=1
plot(input, (A%*%Factor)[index_d,],type='l',ylim=c(min(output_ori[index_d,]),max(output_ori[index_d,])),cex=1.2,xlab='x',ylab=expression(hat(Y)[1]), mgp=c(2.5,1,0))
polygon(c(input,rev(input)),c(Y_95_lower[index_d,],rev(Y_95_upper[index_d,])),col = "grey80", border = F)
lines(input, (A%*%Factor)[index_d,],type='l',ylim=c(min(output_ori[index_d,]),max(output_ori[index_d,])),cex=1.2)
lines(input,Y_hat[index_d,],col='blue',type='l',lty=2,cex=1.5)
lines(input,Y_hat_pc[index_d,],type='l',col='red',lty=3,cex=1.5)
lines(input,Y_hat_QW_1[index_d,],type='l',col='yellow',lty=4,cex=1.5)
lines(input,Y_hat_QW_5[index_d,],type='l',col='green',lty=5,cex=1.5)
legend('topright',legend=c('PCA', 'GPPCA','LY1','LY5'),lty=c(3,2,1),col=c('red','blue','yellow','green'),cex=.6)
lines(input, output[index_d,], type='p',pch=20,cex=.4)

index_d=2
plot(input, (A%*%Factor)[index_d,],type='l',ylim=c(min(output_ori[index_d,]),max(output_ori[index_d,])),cex=1.2,xlab='x',ylab=expression(hat(Y)[2]), mgp=c(2.5,1,0))
polygon(c(input,rev(input)),c(Y_95_lower[index_d,],rev(Y_95_upper[index_d,])),col = "grey80", border = F)
lines(input, (A%*%Factor)[index_d,],type='l',ylim=c(min(output_ori[index_d,]),max(output_ori[index_d,])),cex=1.2)
lines(input,Y_hat[index_d,],col='blue',type='l',lty=2,cex=1.5)
lines(input,(A_est_pc%*%t(A_est_pc)%*%output_ori)[index_d,],type='l',col='red',lty=3,cex=1.5)
legend('topright',legend=c('PCA', 'GPPCA','Truth'),lty=c(3,2,1),col=c('red','blue','black'),cex=.6)
lines(input, output[index_d,], type='p',pch=20,cex=.4)
#dev.off()


par(mfrow=c(1,2))
index_d=3
plot(input, (A%*%Factor)[index_d,],type='l',ylim=c(min(output_ori[index_d,]),max(output_ori[index_d,]+.5)),cex=1.2,xlab='x',ylab=expression(hat(Y)[3]), mgp=c(2.5,1,0))
polygon(c(input,rev(input)),c(Y_95_data_lower[index_d,],rev(Y_95_data_upper[index_d,])),col = "grey90", border = F)
polygon(c(input,rev(input)),c(Y_95_lower[index_d,],rev(Y_95_upper[index_d,])),col = "grey70", border = F)
lines(input, (A%*%Factor)[index_d,],type='l',ylim=c(min(output_ori[index_d,]),max(output_ori[index_d,])),cex=1.2)
lines(input,Y_hat[index_d,],col='blue',type='l',lty=2,cex=1.5)
lines(input,(A_est_pc%*%t(A_est_pc)%*%output_ori)[index_d,],type='l',col='red',lty=3,cex=1.5)
legend('topright',legend=c('PCA', 'GPPCA','Truth'),lty=c(3,2,1),col=c('red','blue','black'),cex=.6)
lines(input, output[index_d,], type='p',pch=20,cex=.4)

index_d=4
plot(input, (A%*%Factor)[index_d,],type='l',ylim=c(min(output_ori[index_d,]),max(output_ori[index_d,])),cex=1.2,xlab='x',ylab=expression(hat(Y)[4]), mgp=c(2.5,1,0))
polygon(c(input,rev(input)),c(Y_95_data_lower[index_d,],rev(Y_95_data_upper[index_d,])),col = "grey90", border = F)

polygon(c(input,rev(input)),c(Y_95_lower[index_d,],rev(Y_95_upper[index_d,])),col = "grey70", border = F)
lines(input, (A%*%Factor)[index_d,],type='l',ylim=c(min(output_ori[index_d,]),max(output_ori[index_d,])),cex=1.2)
lines(input,Y_hat[index_d,],col='blue',type='l',lty=2,cex=1.5)
lines(input,(A_est_pc%*%t(A_est_pc)%*%output_ori)[index_d,],type='l',col='red',lty=3,cex=1.5)
legend('topright',legend=c('PCA', 'GPPCA','Truth'),lty=c(3,2,1),col=c('red','blue','black'),cex=.6)
lines(input, output[index_d,], type='p',pch=20,cex=.4)
#dev.off()










n=400

record_angle_2=matrix(0,N,4)
record_param_2=matrix(0,N,2*d+1)
record_param_pc_2=matrix(0,N,2*d+1)
record_param_QW_1_2=matrix(0,N,2*d+1)
record_param_QW_5_2=matrix(0,N,2*d+1)


#beta_real=rep(0.01,d)
sigma_real=rep(1,d)
#rep(1,d)
#sigma_0_real=.2
#sigma_0_real=sqrt(.01)
sigma_0_real=sqrt(.25)


G=list()
input=seq(1,n,1)
R0_00=abs(outer(input,input,'-'))


k_0=5
G_hat=as.list(1:d)
G_hat_pc=as.list(1:d)
G_hat_QW_1=as.list(1:d)
G_hat_QW_5=as.list(1:d)


record_MSE_2=matrix(0,N,4)
record_prop_diff_2=matrix(0,N,2)
record_prop_diff_data_2=matrix(0,N,2)


for(i_angle in 1:N){
  set.seed(i_angle)
  print(i_angle)
  #set.seed(2)
  #input=seq(0,1,1/(n-1))
  ##
  #min_beta=0.02
  #max_beta=0.005
  beta_real=runif(d,min=0.001,max=.1)
  
  L_sample_list=as.list(1:d)
  
  for(i in 1:d){
    # L_sample_list[[i]]=get_chol(input,beta_real[i])
    # 
    # L_sample_list2=Chol_rcppeigen(R_real)
    # L_sample_list2- t(chol(R_real))
    R_real=matern_5_2_funct(R0_00,beta_real[i])
    L_sample_list[[i]]=Chol_rcppeigen(R_real)
    
  }
  
  
  
  input=sort(input)
  delta_x=input[2:length(input)]-input[1:(length(input)-1)]
  
  #A=matrix(runif(k*d),k,d)
  
  ##if A is not from a Stiefel manifold then the estimation of variance could be wrong
  A=rustiefel(k, d)  ##sample from Stiefel manifold
  
  Factor=matrix(0,d,n)
  
  for(i_fac in 1: d){
    Factor[i_fac,]=sigma_real[i_fac]^2*L_sample_list[[i_fac]]%*%rnorm(n)
  }
  
  output_ori=A%*%Factor+matrix(rnorm(n*k,mean=0,sd=sigma_0_real),k,n)
  
  output=output_ori
  
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
  
  
  A_ini=svd_output$u[,1:d] 
  
  param_ini=c(rep(log(.1),d),rep(log(10),d))
  
  
  #sink(tempfile())
  #options(warn = -1) 
  m=try(optim(param_ini,neg_log_lik_diff_cov_FFBS,A_ini=A_ini,lower=c(rep(log(10^{-10}),d),rep(log(10^{-10}),d)),
              upper=c(rep(log(10^10),d),rep(log(10^10),d)),control = list(maxit = 100), method="L-BFGS-B"),silent=T)
  
  
  
  m_pc=try(optim(param_ini,neg_log_lik_diff_cov_FFBS_fixed_A, fixed_A=A_est_pc,lower=c(rep(log(10^{-10}),d),rep(log(10^{-10}),d)),
                 upper=c(rep(log(10^10),d),rep(log(10^10),d)),control = list(maxit = 100), method="L-BFGS-B"),silent=T)
  
  m_QW_1=try(optim(param_ini,neg_log_lik_diff_cov_FFBS_fixed_A,fixed_A=A_est_QW_1,lower=c(rep(log(10^{-10}),d),rep(log(10^{-10}),d)),
                   upper=c(rep(log(10^10),d),rep(log(10^10),d)),control = list(maxit = 100), method="L-BFGS-B"),silent=T)
  m_QW_5=try(optim(param_ini,neg_log_lik_diff_cov_FFBS_fixed_A, fixed_A=A_est_QW_5,lower=c(rep(log(10^{-10}),d),rep(log(10^{-10}),d)),
                   upper=c(rep(log(10^10),d),rep(log(10^10),d)),control = list(maxit = 100), method="L-BFGS-B"),silent=T)
  
  
  while(is.character(m)){
    param_ini=param_ini+runif(2*d)
    m=try(optim(param_ini,neg_log_lik_diff_cov_FFBS,A_ini=A_ini,lower=c(rep(log(10^{-10}),d),rep(log(10^{-10}),d)),
                upper=c(rep(log(10^10),d),rep(log(10^10),d)),control = list(maxit = 100), method="L-BFGS-B"),silent=T)
  } 
  
  while(is.character(m_pc)){
    param_ini=param_ini+runif(2*d)
    m_pc=try(optim(param_ini,neg_log_lik_diff_cov_FFBS_fixed_A, fixed_A=A_est_pc,lower=c(rep(log(10^{-10}),d),rep(log(10^{-10}),d)),
                   upper=c(rep(log(10^10),d),rep(log(10^10),d)),control = list(maxit = 100), method="L-BFGS-B"),silent=T)
  } 
  
  while(is.character(m_QW_1)){
    param_ini=param_ini+runif(2*d)
    m_QW_1=try(optim(param_ini,neg_log_lik_diff_cov_FFBS_fixed_A,fixed_A=A_est_QW_1,lower=c(rep(log(10^{-10}),d),rep(log(10^{-10}),d)),
                     upper=c(rep(log(10^10),d),rep(log(10^10),d)),control = list(maxit = 100), method="L-BFGS-B"),silent=T)
  } 
  
  while(is.character(m_QW_5)){
    param_ini=param_ini+runif(2*d)
    m_QW_5=try(optim(param_ini,neg_log_lik_diff_cov_FFBS_fixed_A, fixed_A=A_est_QW_5,
                     upper=c(rep(log(10^10),d),rep(log(10^10),d)),lower=c(rep(log(10^{-10}),d),rep(log(10^{-10}),d)),control = list(maxit = 100), method="L-BFGS-B"),silent=T)
  } 
  
  #m=try(optim(param_ini,neg_log_lik_diff_cov_FFBS,A_ini=A_ini,lower=c(rep(log(10^{-10}),d),rep(log(10^{-10}),d))),silent=T)
  
  #sink()
  
  
  ##my method
  beta_hat=exp(m$par[1:d])
  tau_hat=exp(m$par[(d+1):(2*d)])
  
  
  
  G_log_det_cov_hat=Get_G_log_det_cov(m$par, output, delta_x,d,kernel_type = 'matern_5_2')
  
  for(i in 1:d){
    G_hat[[i]]=output_2-G_log_det_cov_hat[[i]]
  }
  
  A_hat=Optimization_Stiefel_Manifold(A_ini, G=G_hat,max_iter=100)
  
  
  sigma_2_0_hat=(trace_output_2+F_Funct(A_hat,G_hat))/(n*k)
  
  ###pc 
  
  G_log_det_cov_pc_hat=Get_G_log_det_cov( m_pc$par, output, delta_x,d,kernel_type = 'matern_5_2')
  
  for(i in 1:d){
    G_hat_pc[[i]]=output_2-G_log_det_cov_pc_hat[[i]]
  }
  
  sigma_2_0_pc_hat=(trace_output_2+F_Funct(A_est_pc,G_hat_pc) )/(n*k)
  
  ##QW_1
  G_log_det_cov_QW_1_hat=Get_G_log_det_cov( m_QW_1$par, output, delta_x,d,kernel_type = 'matern_5_2')
  
  for(i in 1:d){
    G_hat_QW_1[[i]]=output_2-G_log_det_cov_QW_1_hat[[i]]
  }
  
  sigma_2_0_QW_1_hat=(trace_output_2+F_Funct(A_est_QW_1,G_hat_QW_1) )/(n*k)
  
  ##QW_5
  G_log_det_cov_QW_5_hat=Get_G_log_det_cov( m_QW_5$par, output, delta_x,d,kernel_type = 'matern_5_2')
  
  for(i in 1:d){
    G_hat_QW_5[[i]]=output_2-G_log_det_cov_QW_5_hat[[i]]
  }
  
  sigma_2_0_QW_5_hat=(trace_output_2+F_Funct(A_est_QW_5,G_hat_QW_5) )/(n*k)
  
  
  record_angle_2[i_angle,1]=rospca::angle(A,A_est_pc)
  
  record_angle_2[i_angle,2]=rospca::angle(A,A_hat)
  
  
  record_angle_2[i_angle,3]=rospca::angle(A,A_est_QW_1)
  
  record_angle_2[i_angle,4]=rospca::angle(A,A_est_QW_5)
  
  
  
  record_param_2[i_angle,1:d]=beta_hat
  record_param_2[i_angle,(d+1):(2*d)]=tau_hat*(sigma_2_0_hat)
  record_param_2[i_angle,2*d+1]=sigma_2_0_hat
  
  ##pc
  record_param_pc_2[i_angle,1:d]=exp(m_pc$par[1:d])
  record_param_pc_2[i_angle,(d+1):(2*d)]=exp(m_pc$par[(d+1):(2*d) ])*(sigma_2_0_pc_hat)
  record_param_pc_2[i_angle,2*d+1]=sigma_2_0_pc_hat
  
  ##QW_1
  record_param_QW_1_2[i_angle,1:d]=exp(m_QW_1$par[1:d])
  record_param_QW_1_2[i_angle,(d+1):(2*d)]=exp(m_QW_1$par[(d+1):(2*d)])*as.numeric(sigma_2_0_QW_1_hat)
  record_param_QW_1_2[i_angle,2*d+1]=sigma_2_0_QW_1_hat
  
  ##QW_5
  record_param_QW_5_2[i_angle,1:d]=exp(m_QW_5$par[1:d])
  record_param_QW_5_2[i_angle,(d+1):(2*d)]=exp(m_QW_5$par[(d+1):(2*d)])*as.numeric(sigma_2_0_QW_5_hat)
  record_param_QW_5_2[i_angle,2*d+1]=sigma_2_0_QW_5_hat
  
  
  
  print(record_angle_2[i_angle,]) 
  
  ##prediction
  Y_hat_pc=A_est_pc%*%t(A_est_pc)%*%output
  Y_hat_QW_1=A_est_QW_1%*%t(A_est_QW_1)%*%output
  Y_hat_QW_5=A_est_QW_5%*%t(A_est_QW_5)%*%output
  
  pred_GPPCA=pred_FFBS_FastGaSP(param=record_param_1[i_angle,],A_hat=A_hat,input=input,testing_input=input, output_here=output,d=d,var_data=F)
  
  Y_hat=pred_GPPCA[[1]]
  Y_var=pred_GPPCA[[2]]
  Y_95_lower=Y_hat+sqrt(Y_var)*qnorm(0.025)
  Y_95_upper=Y_hat+sqrt(Y_var)*qnorm(0.975)
  
  
  Y_mean=A%*%Factor
  
  record_MSE_2[i_angle,1]=mean( (Y_hat_pc-Y_mean)^2)
  record_MSE_2[i_angle,2]=mean((Y_hat-Y_mean)^2)
  record_MSE_2[i_angle,3]=mean((Y_hat_QW_1-Y_mean)^2)
  record_MSE_2[i_angle,4]=mean((Y_hat_QW_5-Y_mean)^2)
  
  
  record_prop_diff_2[i_angle,1]=sum(Y_95_lower<=Y_mean & Y_95_upper>=Y_mean)/(n*k)
  record_prop_diff_2[i_angle,2]=mean(abs(Y_95_upper-Y_95_lower))
  
  Y_var_data=Y_var+rep(sigma_2_0_hat,k)
  Y_95_data_lower=Y_hat+sqrt(Y_var_data)*qnorm(0.025)
  Y_95_data_upper=Y_hat+sqrt(Y_var_data)*qnorm(0.975)
  
  record_prop_diff_data_2[i_angle,1]=sum(Y_95_data_lower<=output & Y_95_data_upper>=output)/(n*k)
  record_prop_diff_data_2[i_angle,2]=mean(abs(Y_95_data_upper-Y_95_data_lower))
  
  
  print(record_MSE_2[i_angle,])
  
}


record_angle=cbind(record_angle_1,record_angle_2)*pi/2
  
#pdf("simulation_diff_cov_1_sigma_0_0_5.pdf",height=5,width=4)
boxplot(record_angle,las=2,names=c('PCA','GPPCA','LY1','LY5','PCA','GPPCA','LY1','LY5'),
        col=c('red','blue','yellow','green','red','blue','yellow','green'),
        ylab='Largest Principal Angle',ylim=c(0,pi/2),at=c(1,2,3,4,6,7,8,9),main=expression(k==8~and~d==4))
#dev.off()


colMeans(record_angle_1)
colMeans(record_angle_2)
colMeans(record_MSE_1)
colMeans(record_MSE_2)




# ##k=8, d=4
# > colMeans(record_angle_1)
# [1] 0.2740322 0.1770248 0.4165510 0.3582556
# > colMeans(record_angle_2)
# [1] 0.10599092 0.08008695 0.20194743 0.14524094
# > colMeans(record_MSE_1)
# [1] 0.13360965 0.01427039 0.15861892 0.15342771
# > colMeans(record_MSE_2)
# [1] 0.12795342 0.03961855 0.13958358 0.13359466
# ##k=40, d=4
# > colMeans(record_angle_1)
# [1] 0.6228991 0.4077857 0.7333904 0.6094326
# > colMeans(record_angle_2)
# [1] 0.3060351 0.2120718 0.4449162 0.3502915
# > colMeans(record_MSE_1)
# [1] 0.037834892 0.007074085 0.048719902 0.043699084
# > colMeans(record_MSE_2)
# [1] 0.03044574 0.01108526 0.03410951 0.03176962
# ##k=16,d=8
# > colMeans(record_angle_1)
# [1] 0.2976958 0.1665977 0.4190692 0.3391856
# > colMeans(record_angle_2)
# [1] 0.09897629 0.06844120 0.14181685 0.09671753
# > colMeans(record_MSE_1)
# [1] 0.13116487 0.01311457 0.14065226 0.13698913
# > colMeans(record_MSE_2)
# [1] 0.12740170 0.03312051 0.12800191 0.12686855
# ##k=80, d=8
# > colMeans(record_angle_1)
# [1] 0.6011695 0.3909604 0.6958772 0.5541663
# > colMeans(record_angle_2)
# [1] 0.2571782 0.1745981 0.3140106 0.2423956
# > colMeans(record_MSE_1)
# [1] 0.034534719 0.005988451 0.036726032 0.034019493
# > colMeans(record_MSE_2)
# [1] 0.028936833 0.008034399 0.028806492 0.027905387
