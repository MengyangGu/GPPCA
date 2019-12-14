#library(FastGaSP)
library(Rcpp)
library(RcppEigen)
library(RobustGaSP)
library(rstiefel)
library(maps)#Install maps package if not done before
library(pcaMethods)
library(nloptr)
sourceCpp(file='../src/functions.cpp') 
# 
neg_log_lik_with_trend_shared_cov<-function(param,kernel_type){

  beta=exp(param[1])
  tau=exp(param[2])
  
  if(kernel_type=='matern_5_2'){
    Sigma=tau*matern_5_2_funct(R0_00,beta)
  }else if(kernel_type=='exp'){
    Sigma=tau*pow_exp_funct(R0_00,beta,alpha_i=c(1))
  }
  
  Sigma_L_tau=t(L_M)%*%Sigma%*%L_M
  
  Sigma_L_tau_tilde=Sigma_L_tau+diag(num_obs-q)
  Sigma_L_tau_tilde_inv=solve(Sigma_L_tau_tilde)
  matrix1=Sigma_L_tau-Sigma_L_tau%*%Sigma_L_tau_tilde_inv%*%Sigma_L_tau
  
  svd_matrix1=svd(matrix1)
  L_M_tau=svd_matrix1$u%*%sqrt(diag(svd_matrix1$d))
  

  YL=output_L_M%*%L_M_tau
  
  svd_YL=svd(YL)
  
  
  Sigma_tilde=Sigma+diag(num_obs)
  L=t(chol(Sigma_tilde))
  ##one should avoid solve here
  L_H=t(chol(t(H)%*%solve(Sigma_tilde)%*%H))
  sum_det=sum(log(diag(L)))+sum(log(diag(L_H)))
  
  #fixed measurement error
  neg_lik=-(-d*sum_det-(tr_output_M_output-sum(svd_YL$d[1:d]^2))/(2*est_measurement_var) )
  ##est measurement error
  #neg_lik=-(-d*sum_det-(k*(num_obs-q))/2*log(tr_output_M_output-sum(svd_YL$d[1:d]^2)) )
  
  return(neg_lik)
}

##these two functions were not used but can be used if one wants to use the posterior mode instead of the MLE
approx_ref_matern_5_2<-function(a,b,C,beta_i,eta_i ){
  t=C*beta_i+eta_i  ###JR prior
  a*log(t)-b*t
}

neg_log_post_shared_cov<-function(param,kernel_type){
  
  neg_log_lik_val=neg_log_lik_with_trend_shared_cov(param,kernel_type)
  
  beta=exp(param[1])
  tau=exp(param[2])
  eta=1/tau
  a=.2
  b=1
  neglogprior=-approx_ref_matern_5_2(a,b,C,beta,  eta)
  neg_log_lik_val+ neglogprior
}

Get_A_est_sigma_0_2_with_trend_large_k<-function(param,kernel_type){
  beta=exp(param[1])
  tau=exp(param[2])
  
  if(kernel_type=='matern_5_2'){
    Sigma=tau*matern_5_2_funct(R0_00,beta)
  }else if(kernel_type=='exp'){
    Sigma=tau*pow_exp_funct(R0_00,beta,alpha_i=c(1))
  }
  
  Sigma_L_tau=t(L_M)%*%Sigma%*%L_M
  
  Sigma_L_tau_tilde=Sigma_L_tau+diag(num_obs-q)
  Sigma_L_tau_tilde_inv=solve(Sigma_L_tau_tilde)
  matrix1=Sigma_L_tau-Sigma_L_tau%*%Sigma_L_tau_tilde_inv%*%Sigma_L_tau
  
  svd_matrix1=svd(matrix1)
  L_M_tau=svd_matrix1$u%*%sqrt(diag(svd_matrix1$d))
  
  YL=output_L_M%*%L_M_tau
  
  svd_YL=svd(YL)
  
  return_list=as.list(1:3)
  
  return_list[[1]]=svd_YL$u[,1:d]
  return_list[[2]]=(tr_output_M_output-sum(svd_YL$d[1:d]^2))/( (num_obs-q)*k)
  return_list[[3]]=svd_YL$d^2
  return_list
  
}


pred_mean_share_cov_with_trend<-function(param,sigma_2_0_here,A_hat,output_here,kernel_type){
  beta=exp(param[1])
  tau=exp(param[2])
  
  Z_hat=matrix(0,d,num_obs)
  z_star=matrix(0,d,num_obs_all)
  D_matrix=matrix(0,d,num_obs_all)
  
  sigma_2_0_here=as.numeric(sigma_2_0_here)
  sigma_2=(tau*sigma_2_0_here)
  if(kernel_type=='matern_5_2'){
    Sigma=sigma_2*matern_5_2_funct(R0_00,beta)
    Sigma_star=sigma_2*matern_5_2_funct(r0, beta)
  }else if(kernel_type=='exp'){
    Sigma=sigma_2*pow_exp_funct(R0_00,beta,alpha_i=c(1))
    Sigma_star=sigma_2*pow_exp_funct(r0,beta,alpha_i=c(1))
  }
  
  middle_part=t(A_hat)%*%output_here%*%M%*%solve(Sigma%*%M+sigma_2_0_here*diag(num_obs))
  Z_hat=middle_part%*%Sigma
  
  z_star=t(Sigma_star)%*%t(middle_part)
  
  ##calculating the variance
  Sigma_tilde=Sigma+sigma_2_0_here*diag(num_obs)
  Sigma_tilde_inv=solve(Sigma_tilde)
  H_t_inv_Sigma_tilde_H_inv=solve((t(H)%*%solve(Sigma_tilde)%*%H))
  
  Sigma_star_Sigma_tilde_inv_Sigma_star=t(Sigma_star)%*%Sigma_tilde_inv%*%Sigma_star
  h_t_star_minus_H_t_Sigma_inv_Sigma=t(H_testing)-t(H)%*%Sigma_tilde_inv%*%Sigma_star
  D_star=rep(0,(num_obs_all)) 
  for(i_test in 1:(num_obs_all)){
    D_star[i_test]=-t(Sigma_star[,i_test])%*%Sigma_tilde_inv%*%Sigma_star[,i_test]+t(h_t_star_minus_H_t_Sigma_inv_Sigma[,i_test])%*%(H_t_inv_Sigma_tilde_H_inv)%*%h_t_star_minus_H_t_Sigma_inv_Sigma[,i_test]
  }
  
  D_star=D_star+sigma_2+sigma_2_0_here
  
  pred_mean=(output_here-A_hat%*% Z_hat)%*%H%*%solve(t(H)%*%H)%*%t(H_testing)+A_hat%*%t(z_star)
  
  pred_list=as.list(1:2)
  pred_list[[1]]=pred_mean
  pred_list[[2]]=D_star
  
  return(pred_list)
}

    
gridded_data=readRDS("NOAA_gridded_temperature/NOAA_gridded_data_1999_2018.rds")
 output_all=gridded_data$temperature

grid_with_full_obs=gridded_data$loc_full_obs
time=gridded_data$months
Lat= seq(-87.5, 87.5, length=36)
Lon=seq(2.5, 357.5, length=72)


set.seed(1)  
num_obs_all=dim(output_all)[2]
k=dim(output_all)[1]

input_all=as.numeric(1:(num_obs_all))
index_testing_month=sample(input_all, 20)
index_testing_loc=sample(1:dim(output_all)[1], 1200)



input=input_all[-index_testing_month]
testing_input=input_all[index_testing_month]
index_training_loc=(1:k)[-index_testing_loc]


output=output_all[,-index_testing_month]
testing_output=output_all[index_testing_loc,index_testing_month]
output_partial_obs=output_all[-index_testing_loc,index_testing_month]


num_obs=length(input)
k=dim(output)[1]

C=(max(input)-min(input))/num_obs

q=2
H=matrix(0,num_obs,q)
H[,1]=rep(1,num_obs)
H[,2]=input
HH_inv=solve(t(H)%*%H)
HH=H%*%HH_inv%*%t(H)
M=diag(num_obs)-HH
eigen_M=eigen(M)
#q=dim(H)[2]
L_M=eigen_M$vector%*%diag(sqrt(abs(eigen_M$values) ))[,1:(num_obs-q)]
output_L_M=output%*%L_M

tr_output_M_output=sum(diag(t(output_L_M)%*%output_L_M))

###this value is the measurement var fixed or estimated from the station
est_measurement_var=0.1

d=100


R0_00=(abs(outer(input,input,'-')))

param_ini=c(0,10)



kernel_type='matern_5_2'

##maximum posterior mode using the JR prior. This one is more robust in optimization.
#neg_log_post_shared_cov(param_ini,kernel_type)

#m_trend_post=try(optim(param_ini,neg_log_post_shared_cov,kernel_type=kernel_type,lower=rep(log(10^{-10}),2),
#                       upper=rep(log(10^{10}),2),control = list(maxit = 100), method="L-BFGS-B"),silent=T)

## maximum marginal likelihood estimation
m_lik=try(optim(param_ini,neg_log_lik_with_trend_shared_cov,kernel_type=kernel_type,lower=rep(log(10^{-10}),2),
                upper=rep(log(10^{10}),2),control = list(maxit = 100), method="L-BFGS-B"),silent=T)


Get_A_est_all=Get_A_est_sigma_0_2_with_trend_large_k(m_lik$par,kernel_type=kernel_type)

A_with_mean_est=Get_A_est_all[[1]]
##estimated variance
#sigma_2_0_est=Get_A_est_all[[2]]
##fixed variance
sigma_2_0_est=est_measurement_var


r0=abs(outer(input, input_all,'-'))

H_testing=matrix(1,num_obs_all,q)
H_testing[,2]=input_all
H_testing=as.matrix(H_testing)


pred_GPPCA_with_mean_all=pred_mean_share_cov_with_trend(param=m_lik$par,sigma_2_0_here=sigma_2_0_est, A_hat=A_with_mean_est,output_here=output,kernel_type=kernel_type)


sqrt(mean( (pred_GPPCA_with_mean_all[[1]][index_testing_loc,index_testing_month]-testing_output)^2))

sd(as.numeric(testing_output))


AA_t=A_with_mean_est%*%t(A_with_mean_est)
I_minus_AA_t=diag(k)-AA_t
D_star=pred_GPPCA_with_mean_all[[2]]
pred_GPPCA_conditional_mean=matrix(0, length(index_testing_loc),num_obs_all-num_obs )
pred_GPPCA_conditional_var=matrix(0, length(index_testing_loc),num_obs_all-num_obs)

diag_Sigma_est_matrix=matrix(0,k,  length(index_testing_month))
for(i_testing in 1:length(index_testing_month) ){
  print(i_testing)
  Sigma_est=D_star[index_testing_month[i_testing] ]*(AA_t)+(sigma_2_0_est*as.numeric(1+t(H_testing[index_testing_month[i_testing],])%*%HH_inv%*%(H_testing[index_testing_month[i_testing],])  )*I_minus_AA_t)
  #Sigma_est=D_star[index_testing_year[i_testing] ]*(AA_t)
  #diag(Sigma_est)=diag(Sigma_est)+(sigma_2_0_est*as.numeric(1+t(H_testing[index_testing_year[i_testing],])%*%HH_inv%*%(H_testing[index_testing_year[i_testing],])  ))
  #Sigma_est=D_star[index_testing_year[i_testing] ]*(AA_t)-(sigma_2_0_est*as.numeric(1+t(H_testing[index_testing_year[i_testing],])%*%HH_inv%*%(H_testing[index_testing_year[i_testing],])  )*AA_t)+sigma_2_0_est*diag(k)
  diag_Sigma_est_matrix[,i_testing]=diag(Sigma_est)
  
  Sigma_est_00_inv=solve(Sigma_est[-index_testing_loc,-index_testing_loc])
  Useful_block=Sigma_est[index_testing_loc,-index_testing_loc]%*%Sigma_est_00_inv
  pred_GPPCA_conditional_mean[,i_testing]=pred_GPPCA_with_mean_all[[1]][index_testing_loc,index_testing_month[i_testing]]+Useful_block%*%( output_partial_obs[,i_testing]-pred_GPPCA_with_mean_all[[1]][-index_testing_loc,index_testing_month[i_testing]])
  
  pred_GPPCA_conditional_var[,i_testing]=diag(Sigma_est[index_testing_loc,index_testing_loc])-diag(Useful_block%*%(Sigma_est[-index_testing_loc,index_testing_loc]))
}

LB_95_GPPCA=pred_GPPCA_conditional_mean+sqrt(pred_GPPCA_conditional_var)*qnorm(0.025,0,1)
UB_95_GPPCA=pred_GPPCA_conditional_mean+sqrt(pred_GPPCA_conditional_var)*qnorm(0.975,0,1)

sqrt(mean( (pred_GPPCA_conditional_mean-testing_output)^2))
mean(abs(UB_95_GPPCA> testing_output & testing_output>LB_95_GPPCA))
mean(abs(UB_95_GPPCA-LB_95_GPPCA))


# #fixed measurement error, d=50, matern
# 
#> sqrt(mean( (pred_GPPCA_conditional_mean-testing_output)^2))
#[1] 0.3851825
#> mean(abs(UB_95_GPPCA> testing_output & testing_output>LB_95_GPPCA))
#[1] 0.9332083
#> mean(abs(UB_95_GPPCA-LB_95_GPPCA))
#[1] 1.332343
#fixed measurement error, d=100, matern
# > sqrt(mean( (pred_GPPCA_conditional_mean-testing_output)^2))
# [1] 0.3148709
# > mean(abs(UB_95_GPPCA> testing_output & testing_output>LB_95_GPPCA))
# [1] 0.97725
# > mean(abs(UB_95_GPPCA-LB_95_GPPCA))
# [1] 1.444499

# estimated measurement error, d=50, matern
#> sqrt(mean( (pred_GPPCA_conditional_mean-testing_output)^2))
#[1] 0.3859712
#> mean(abs(UB_95_GPPCA> testing_output & testing_output>LB_95_GPPCA))
#[1] 0.8698333
#> mean(abs(UB_95_GPPCA-LB_95_GPPCA))
#[1] 1.017599
# # estimated measurement error, d=100, matern
#> sqrt(mean( (pred_GPPCA_conditional_mean-testing_output)^2))
#[1] 0.3202735
#> mean(abs(UB_95_GPPCA> testing_output & testing_output>LB_95_GPPCA))
#[1] 0.77175
#> mean(abs(UB_95_GPPCA-LB_95_GPPCA))
#[1] 0.5632038



##1. GP to interpolate the grid for each missing year
library(RobustGaSP)
##here the grid is 72*36
input_spatial=matrix(0,dim(output_partial_obs)[1],2)
testing_input_spatial=matrix(0,k- dim(output_partial_obs)[1],2)

input_spatial[,1]=index_training_loc%%72
testing_input_spatial[,1]=index_testing_loc%%72

input_spatial[,2]=ceiling(index_training_loc/72)
testing_input_spatial[,2]=ceiling(index_testing_loc/72)

predict_rgasp_mean=matrix(0,dim(testing_output)[1],dim(testing_output)[2])
predict_rgasp_sd=matrix(0,dim(testing_output)[1],dim(testing_output)[2])
predict_rgasp_lower_95=matrix(0,dim(testing_output)[1],dim(testing_output)[2])
predict_rgasp_upper_95=matrix(0,dim(testing_output)[1],dim(testing_output)[2])

for(i in 1: length(index_testing_month) ){
  print(i)
  model.rgasp=rgasp(design=input_spatial,response=as.matrix(output_partial_obs[,i]),nugget.est=T)
  pred.rgasp.all=predict(model.rgasp,testing_input_spatial)
  predict_rgasp_mean[,i]=pred.rgasp.all$mean
  predict_rgasp_lower_95[,i]=pred.rgasp.all$lower95
  predict_rgasp_upper_95[,i]=pred.rgasp.all$upper95
  predict_rgasp_sd[,i]=pred.rgasp.all$sd
  print(sqrt(mean((pred.rgasp.all[[1]]-testing_output[,i])^2)))
}


(sqrt(mean((predict_rgasp_mean-testing_output)^2)))
mean(abs(predict_rgasp_upper_95> testing_output & testing_output>predict_rgasp_lower_95))
mean(abs(predict_rgasp_upper_95-predict_rgasp_lower_95))


###spatial grid
# > (sqrt(mean((predict_rgasp_mean-testing_output)^2)))
# [1] 0.5595441
# > mean(abs(predict_rgasp_upper_95> testing_output & testing_output>predict_rgasp_lower_95))
# [1] 0.9418333
# > mean(abs(predict_rgasp_upper_95-predict_rgasp_lower_95))
# [1] 2.230837


#i_plot=500
#plot(input_all,output_all[index_testing_loc[i_plot],])
#lines(testing_input,predict_rgasp_mean[i_plot,],col='red',type='p')
#lines(testing_input,pred_GPPCA_conditional_mean[i_plot,],col='blue',type='p')
#sqrt(mean( (pred_GPPCA_conditional_mean[i_plot,]-output_all[index_testing_loc[i_plot],index_testing_month])^2))
#sqrt(mean( (predict_rgasp_mean[i_plot,]-output_all[index_testing_loc[i_plot],index_testing_month])^2))

##2. GP to interpolate the grid for each missing year

predict_rgasp_mean_month=matrix(0,dim(testing_output)[1],dim(testing_output)[2])
predict_rgasp_sd_month=matrix(0,dim(testing_output)[1],dim(testing_output)[2])
predict_rgasp_lower_95_month=matrix(0,dim(testing_output)[1],dim(testing_output)[2])
predict_rgasp_upper_95_month=matrix(0,dim(testing_output)[1],dim(testing_output)[2])

for(i_test in 1: length(index_testing_loc) ){
  print(i_test)
  model.rgasp=rgasp(design=input,response=as.matrix(output[index_testing_loc[i_test],]),trend=H,nugget.est=T)
  
  pred.rgasp.all=predict(model.rgasp,testing_input=as.matrix(testing_input), testing_trend=H_testing[index_testing_month,])
  
  predict_rgasp_mean_month[i_test,]=pred.rgasp.all$mean
  predict_rgasp_lower_95_month[i_test,]=pred.rgasp.all$lower95
  predict_rgasp_upper_95_month[i_test,]=pred.rgasp.all$upper95
  predict_rgasp_sd_month[i_test,]=pred.rgasp.all$sd
  
  
}



sqrt(mean((predict_rgasp_mean_month-testing_output)^2))
mean(abs(predict_rgasp_upper_95_month> testing_output & testing_output>predict_rgasp_lower_95_month))
mean(abs(predict_rgasp_upper_95_month-predict_rgasp_lower_95_month))


#####by month
# > sqrt(mean((predict_rgasp_mean_month-testing_output)^2))
# [1] 0.9367765
# > mean(abs(predict_rgasp_upper_95_month> testing_output & testing_output>predict_rgasp_lower_95_month))
# [1] 0.944
# > mean(abs(predict_rgasp_upper_95_month-predict_rgasp_lower_95_month))
# [1] 2.27591


# i_plot=300
# plot(input_all,output_all[index_testing_loc[i_plot],])
# lines(testing_input,predict_rgasp_mean_month[i_plot,],col='red',type='p')
# lines(testing_input,pred_GPPCA_conditional_mean[i_plot,],col='blue',type='p')
# 
# sqrt(mean((predict_rgasp_mean_month[i_plot,]-testing_output[i_plot,])^2))
# sqrt(mean((pred_GPPCA_conditional_mean[i_plot,]-testing_output[i_plot,])^2))
# 
# 
# max(pred_GPPCA_conditional_mean-testing_output)
#lines(testing_input,predict_rgasp_mean[i_plot,],col='red',type='p')
#lines(testing_input,pred_GPPCA_conditional_mean[i_plot,],col='blue',type='p')
#sqrt(mean( (pred_GPPCA_conditional_mean[i_plot,]-output_all[index_testing_loc[i_plot],index_testing_month])^2))
#sqrt(mean( (predict_rgasp_mean[i_plot,]-output_all[index_testing_loc[i_plot],index_testing_month])^2))

##ppca

#plot(svd_output_sub_mean_div_sqrt_n$d)

d=100

row_mean_output=rowMeans(output)

output_sub_mean=(output-row_mean_output)


svd_output_sub_mean_div_sqrt_n=svd(output_sub_mean/sqrt(num_obs))

#S=output_sub_mean%*%t(output_sub_mean)/num_obs
#eigen_S=eigen(S)
#svd_output_sub_mean_div_sqrt_n$u[,1]-eigen_S$vectors[,1]
#svd_output_sub_mean_div_sqrt_n$u[1:5,1]
#eigen_S$vectors[1:5,1]
#eigen_S$values[1:5]
#svd_output_sub_mean_div_sqrt_n$d[1:5]^2
###fixed measurement error
est_error_var=0.1
##estimated measurement error
#est_error_var=(sum(svd_output_sub_mean_div_sqrt_n$d[ (d+1): length(svd_output_sub_mean_div_sqrt_n$d)]))/(k-d-1)
A_hat_ppca=svd_output_sub_mean_div_sqrt_n$u[,1:d]%*%sqrt( diag( (svd_output_sub_mean_div_sqrt_n$d^2-est_error_var)[1:d] ) )

Var_ppca=A_hat_ppca%*%t(A_hat_ppca)
diag(Var_ppca)=diag(Var_ppca)+est_error_var


Var_ppca_00_inv=solve(Var_ppca[-index_testing_loc,-index_testing_loc])
Useful_block=Var_ppca[index_testing_loc,-index_testing_loc]%*%Var_ppca_00_inv

pred_ppca_conditional_mean=Useful_block%*%output_partial_obs+matrix(row_mean_output[index_testing_loc],dim(testing_output)[1],dim(testing_output)[2])

conditional_var_ppca=diag(Var_ppca[index_testing_loc,index_testing_loc])-diag(Useful_block%*%(Var_ppca[-index_testing_loc,index_testing_loc]))

LB_95_ppca=pred_ppca_conditional_mean+matrix(sqrt(conditional_var_ppca)*qnorm(0.025,0,1),length(index_testing_loc),length(index_testing_month))
UB_95_ppca=pred_ppca_conditional_mean+matrix(sqrt(conditional_var_ppca)*qnorm(0.975,0,1),length(index_testing_loc),length(index_testing_month))


sqrt(mean( (pred_ppca_conditional_mean-testing_output)^2))
mean(abs(UB_95_ppca> testing_output & testing_output>LB_95_ppca))
mean(abs(UB_95_ppca-LB_95_ppca))

# #rospca::angle(Get_A_est_all[[1]],A_hat_ppca)
# plot(Get_A_est_all[[3]]/sum(Get_A_est_all[[3]]))
# 
# ppca_eigenvalue=svd_output_sub_mean_div_sqrt_n$d^2-est_error_var
# plot( ppca_eigenvalue/sum(ppca_eigenvalue))
# 
# sum( (Get_A_est_all[[3]]/sum(Get_A_est_all[[3]]))[1:100])
# sum( (ppca_eigenvalue/sum(ppca_eigenvalue))[1:100])

# ##ppca estimated measurement error, d=50
# > sqrt(mean( (pred_ppca_conditional_mean-testing_output)^2))
# [1] 0.6203872
# > mean(abs(UB_95_ppca> testing_output & testing_output>LB_95_ppca))
# [1] 0.6772083
# > mean(abs(UB_95_ppca-LB_95_ppca))
# [1] 1.082893
# ##ppca fixed measurement error, d=100
# 
# > sqrt(mean( (pred_ppca_conditional_mean-testing_output)^2))
# [1] 0.6024701
# > mean(abs(UB_95_ppca> testing_output & testing_output>LB_95_ppca))
# [1] 0.525125
# > mean(abs(UB_95_ppca-LB_95_ppca))
# [1] 0.8029068
# 
# ##ppca fixed measurement error, d=50
# 
# > sqrt(mean( (pred_ppca_conditional_mean-testing_output)^2))
# [1] 0.6173095
# > mean(abs(UB_95_ppca> testing_output & testing_output>LB_95_ppca))
# [1] 0.7654583
# > mean(abs(UB_95_ppca-LB_95_ppca))
# [1] 1.321161
# 
# ##ppca fixed measurement error, d=100
# 
# > sqrt(mean( (pred_ppca_conditional_mean-testing_output)^2))
# [1] 0.5852149
# > mean(abs(UB_95_ppca> testing_output & testing_output>LB_95_ppca))
# [1] 0.819
# > mean(abs(UB_95_ppca-LB_95_ppca))
# [1] 1.396627
# 
###PCA  
# row_mean_output=rowMeans(output)
# 
# output_sub_mean=(output-row_mean_output)
# 
# svd_output_sub_mean_div_sqrt_n=svd(output_sub_mean)
# A_hat_pca=svd_output_sub_mean_div_sqrt_n$u[,1:d]
# 
# Var_pca=A_hat_pca%*%t(A_hat_pca)
# diag(Var_pca)=diag(Var_pca)+est_measurement_var
# 
# 
# Var_pca_00_inv=solve(Var_pca[-index_testing_loc,-index_testing_loc])
# Useful_block_pca=Var_pca[index_testing_loc,-index_testing_loc]%*%Var_pca_00_inv
# 
# pred_pca_conditional_mean=Useful_block_pca%*%output_partial_obs+row_mean_output[index_testing_loc]
# 
# sqrt(mean( (pred_pca_conditional_mean-testing_output)^2))

###random forest by time
dim(output)
dim(output_partial_obs)
library(randomForest)

input_loc=output[-index_testing_loc,]
pred_rf_time_record=matrix(0,length(index_testing_loc),length(index_testing_month) )
for(i_test in 1:length(index_testing_month)){
  print(i_test)
  output_loc=output_partial_obs[,i_test]
  input_testing_loc=output[index_testing_loc,]
  rf_loc=randomForest(x=(input_loc),y=as.matrix(output_loc))
  rf_loc_pred=predict(rf_loc,input_testing_loc)
  pred_rf_time_record[,i_test]=rf_loc_pred
}

sqrt( mean( (pred_rf_time_record-testing_output)^2))

###random forest by time 
# > sqrt( mean( (pred_rf_time_record-testing_output)^2))
# [1] 0.4412568

dim(output)
dim(output_partial_obs)
library(randomForest)

input_time=output[-index_testing_loc,]
pred_rf_loc_record=matrix(0,length(index_testing_loc),length(index_testing_month) )
for(i_test in 1:length(index_testing_loc)){
  print(i_test)
  output_time=output[index_testing_loc[i_test],]
  input_testing_loc=output_partial_obs
  rf_loc=randomForest(x=t(input_time),y=as.matrix(output_time))
  rf_loc_pred=predict(rf_loc,t(input_testing_loc))
  pred_rf_loc_record[i_test,]=rf_loc_pred
}

sqrt( mean( (pred_rf_loc_record-testing_output)^2))

# > sqrt( mean( (pred_rf_loc_record-testing_output)^2))
# [1] 0.3909926



# Spatial temporal
R0_spatial=as.list(1:2)
R0_spatial[[1]]=abs(outer(input_spatial[,1],input_spatial[,1],'-'))
R0_spatial[[2]]=abs(outer(input_spatial[,2],input_spatial[,2],'-'))

R0_temporal=abs(outer(input_all,input_all,'-'))
k_obs=dim(input_spatial)[1]

##spatial temporal output for estimation parameter
output_st=output_all[-index_testing_loc,]

##only use the baseline  mean, otherwise the computation is hard
X_spatial=rep(1,k_obs)

neg_log_lik_spatial_temporal<-function(param){
  
  beta_spatial=exp(param[1:2])
  beta_temporal=exp(param[3])
  tau_spatial=exp(param[4])
  tau_temporal=exp(param[5])
  
  R_spatial=separable_kernel(R0_spatial,beta_spatial,'matern_5_2',rep(1,2))
  
  R_temporal=matern_5_2_funct(R0_temporal,beta_temporal)
  
  
  R_spatial_tilde=R_spatial+tau_spatial*diag(k_obs)
  R_temporal_tilde=R_temporal+tau_temporal*diag(num_obs_all)
  
  
  ##we should use cholesky to make it more robust it the next round
  #L_spatial=t(chol(R_spatial_tilde))
  #L_temporal=t(chol(R_temporal_tilde))
  R_spatial_tilde_inv=solve(R_spatial_tilde)
  R_temporal_tilde_inv=solve(R_temporal_tilde)
  
  R_inv_X_spatial=R_spatial_tilde_inv%*%X_spatial
  

  X_spatial_t_R_inv_X_spatial_inv=solve(t(X_spatial)%*%R_inv_X_spatial)
  beta_hat=X_spatial_t_R_inv_X_spatial_inv%*%t(R_inv_X_spatial)%*%(output_st)
  output_st_tilde= output_st-(X_spatial%*%beta_hat)
  S_2= sum(diag(R_spatial_tilde_inv%*%output_st_tilde%*%R_temporal_tilde_inv%*%t(output_st_tilde)))
  -(-num_obs_all/2*determinant(R_spatial_tilde)$modulus[1]-k_obs/2*determinant(R_temporal_tilde)$modulus[1] -k_obs*num_obs_all/2*log(S_2))
}

#neg_log_lik_spatial_temporal(c(-1,-1,-1,-3,-3))

##run  five times to find the mle
##seems to give me the same result
for(i in 1:5){
  set.seed(i)
  par_ini=-5*runif(5)
  
  par_ini=c(-5*runif(3),-10*runif(2))
  
  m0=optim(par_ini,neg_log_lik_spatial_temporal,method="L-BFGS-B")
  print(m0$par)
  print(m0$value)
  if(i==1){
    m=m0
  }else if(m0$value<m$value){
    m=m0
  }
}

input_temporal=input
testing_input_temporal=testing_input
#num_testing_temporal=length(testing_input_temporal)

k_testing=dim(testing_input_spatial)[1]
X_testing_spatial=rep(1,k_testing)




beta_spatial=exp(m$par[1:2])
beta_temporal=exp(m$par[3])
tau_spatial=exp(m$par[4])
tau_temporal=exp(m$par[5])

R_spatial=separable_kernel(R0_spatial,beta_spatial,'matern_5_2',rep(1,2))

R_temporal=matern_5_2_funct(R0_temporal,beta_temporal)


R_spatial_tilde=R_spatial+tau_spatial*diag(k_obs)
R_temporal_tilde=R_temporal+tau_temporal*diag(num_obs_all)

##we should use cholesky to make it more robust it the next round
R_spatial_tilde_inv=solve(R_spatial_tilde)
R_temporal_tilde_inv=solve(R_temporal_tilde)


R_inv_X_spatial=R_spatial_tilde_inv%*%X_spatial


X_spatial_t_R_inv_X_spatial_inv=solve(t(X_spatial)%*%R_inv_X_spatial)
beta_hat=X_spatial_t_R_inv_X_spatial_inv%*%t(R_inv_X_spatial)%*%(output_st)
output_st_tilde= output_st-(X_spatial%*%beta_hat)
S_2= sum(diag(R_spatial_tilde_inv%*%output_st_tilde%*%R_temporal_tilde_inv%*%t(output_st_tilde)))

sigma_hat_spatial_temporal=S_2/(k_obs*num_obs_all)

r0_spatial=as.list(1:2)
r0_spatial[[1]]=abs(outer(input_spatial[,1],testing_input_spatial[,1],'-'))
r0_spatial[[2]]=abs(outer(input_spatial[,2],testing_input_spatial[,2],'-'))

r0_temporal=abs(outer(input,testing_input_temporal,'-'))


r_spatial=separable_kernel(r0_spatial,beta_spatial,'matern_5_2',rep(1,2))

r_temporal=matern_5_2_funct(r0_temporal,beta_temporal)

rr0_spatial=as.list(1:2)
rr0_spatial[[1]]=abs(outer(testing_input_spatial[,1],testing_input_spatial[,1],'-'))
rr0_spatial[[2]]=abs(outer(testing_input_spatial[,2],testing_input_spatial[,2],'-'))
rr_spatial=separable_kernel(rr0_spatial,beta_spatial,'matern_5_2',rep(1,2))



pred_mean_spatial=X_testing_spatial%*%beta_hat+t(r_spatial)%*%R_spatial_tilde_inv%*% (output_st-(X_spatial%*%beta_hat))

#sqrt(mean((pred_mean_spatial[,index_testing_month]-testing_output)^2))


R_temporal_tilde_part=R_temporal_tilde[-index_testing_month,-index_testing_month]

R_temporal_tilde_part_inv=solve(R_temporal_tilde_part)
pred_mean_spatial_temporal=pred_mean_spatial[,index_testing_month]+t(t(r_temporal)%*%R_temporal_tilde_part_inv%*%t(output_all[index_testing_loc,-index_testing_month]-pred_mean_spatial[,-index_testing_month]))

R_star_tilde_temporal=R_temporal_tilde[index_testing_month,index_testing_month]-t(r_temporal)%*%R_temporal_tilde_part_inv%*%r_temporal
R_star_tilde_spatial=rr_spatial+tau_spatial*diag(k_testing)-t(r_spatial)%*%R_spatial_tilde_inv%*%r_spatial

var_est=sigma_hat_spatial_temporal*kronecker(diag(R_star_tilde_temporal),diag(R_star_tilde_spatial))

var_est_matrix=matrix(var_est,k_testing,num_testing_temporal)

LB_spatial_temporal=pred_mean_spatial_temporal+sqrt(var_est_matrix)*qnorm(0.025)
UB_spatial_temporal=pred_mean_spatial_temporal+sqrt(var_est_matrix)*qnorm(0.975)

sqrt(mean((pred_mean_spatial_temporal-testing_output)^2))


mean(abs(UB_spatial_temporal> testing_output & testing_output>LB_spatial_temporal))
mean(abs(UB_spatial_temporal-LB_spatial_temporal))


#> sqrt(mean((pred_mean_spatial_temporal-testing_output)^2))
#[1] 0.4920091
#> mean(abs(UB_spatial_temporal> testing_output & testing_output>LB_spatial_temporal))
#[1] 0.956625
#> mean(abs(UB_spatial_temporal-LB_spatial_temporal))
#[1] 2.103004

##make some plots 
library(maps)#Install maps package if not done before

#par(mfrow=c(1,3))
#i_plot=5
i_plot=5
mapmat=matrix(NA,72,36)
mapmat[grid_with_full_obs ]=output_all[,index_testing_month[i_plot]]
max(output_all[,index_testing_month[i_plot]])
min(output_all[,index_testing_month[i_plot]])
int=seq(-9.05,9.05,length.out=81)
max(mapmat, na.rm = T)
min(mapmat, na.rm = T)

pdf(file='real_temp_2016_Nov.pdf',height=6.2,width=6.2)
rgb.palette=colorRampPalette(c('black','blue', 'darkgreen','green', 'yellow','pink','red','maroon'),interpolate='spline')
filled.contour(Lon, Lat, mapmat, color.palette=rgb.palette, levels=int,
               plot.title=title(main="Observated temperature anomalies",
                                xlab="Longitude",ylab="Latitude", cex.lab=1.5),
               plot.axes={axis(1, cex.axis=1.5);axis(2, cex.axis=1.5);map('world2', add=TRUE);grid()},
               key.title=title(main=~degree*C),
               key.axes={axis(4, cex.axis=1.5)})
dev.off()

mapmat=matrix(NA,72,36)
mapmat[grid_with_full_obs[-index_testing_loc] ]=output_all[-index_testing_loc,index_testing_month[i_plot]]

mapmat[grid_with_full_obs[index_testing_loc] ]=pred_GPPCA_conditional_mean[,i_plot]
int=seq(-9.05,9.05,length.out=81)
max(mapmat, na.rm = T)
min(mapmat, na.rm = T)

pdf(file='GPPCA_est_temp_2016_Nov.pdf',height=6.2,width=6.2)
rgb.palette=colorRampPalette(c('black','blue', 'darkgreen','green', 'yellow','pink','red','maroon'),interpolate='spline')
filled.contour(Lon, Lat, mapmat, color.palette=rgb.palette, levels=int,
               plot.title=title(main="Interpolation by the GPPCA",
                                xlab="Longitude",ylab="Latitude", cex.lab=1.5),
               plot.axes={axis(1, cex.axis=1.5);axis(2, cex.axis=1.5);map('world2', add=TRUE);grid()},
               key.title=title(main=~degree*C),
               key.axes={axis(4, cex.axis=1.5)})
dev.off()


mapmat=matrix(NA,72,36)
mapmat[grid_with_full_obs[-index_testing_loc] ]=output_all[-index_testing_loc,index_testing_month[i_plot]]
mapmat[grid_with_full_obs[index_testing_loc] ]=pred_mean_spatial_temporal[,i_plot]
int=seq(-9.05,9.05,length.out=81)

max(mapmat, na.rm = T)
min(mapmat, na.rm = T)

pdf(file='spatio_temporal_smoothing_temp_2016_Nov.pdf',height=6.2,width=6.2)
rgb.palette=colorRampPalette(c('black','blue', 'darkgreen','green', 'yellow','pink','red','maroon'),interpolate='spline')
filled.contour(Lon, Lat, mapmat, color.palette=rgb.palette, levels=int,
               plot.title=title(main="Interpolation by the spatio-temporal model",
                                xlab="Longitude",ylab="Latitude", cex.lab=1.5),
               plot.axes={axis(1, cex.axis=1.5);axis(2, cex.axis=1.5);map('world2', add=TRUE);grid()},
               key.title=title(main=~degree*C),
               key.axes={axis(4, cex.axis=1.5)})
dev.off()

#i_plot=5
index_testing_month[i_plot]
sqrt(mean( (pred_GPPCA_conditional_mean[,i_plot]-testing_output[,i_plot])^2))
sqrt(mean( (predict_rgasp_mean[,i_plot]-testing_output[,i_plot])^2))
sqrt(mean( (pred_ppca_conditional_mean[,i_plot]-testing_output[,i_plot])^2))



mapmat=matrix(NA,72,36)
mapmat[grid_with_full_obs[-index_testing_loc] ]=output_all[-index_testing_loc,index_testing_month[i_plot]]
mapmat[grid_with_full_obs[index_testing_loc] ]=pred_ppca_conditional_mean[,i_plot]
int=seq(-9.05,9.05,length.out=81)

max(mapmat, na.rm = T)
min(mapmat, na.rm = T)

rgb.palette=colorRampPalette(c('black','blue', 'darkgreen','green', 'yellow','pink','red','maroon'),interpolate='spline')
filled.contour(Lon, Lat, mapmat, color.palette=rgb.palette, levels=int,
               plot.title=title(main="Interpolation by the PPCA",
                                xlab="Longitude",ylab="Latitude", cex.lab=1.5),
               plot.axes={axis(1, cex.axis=1.5);axis(2, cex.axis=1.5);map('world2', add=TRUE);grid()},
               key.title=title(main=~degree*C),
               key.axes={axis(4, cex.axis=1.5)})


# 
# 
# #linear model
# input_loc=output[-index_testing_loc,]
# pred_lm_loc_record=matrix(0,length(index_testing_loc),length(index_testing_month) )
# for(i_test in 1:length(index_testing_month)){
#   print(i_test)
#   output_loc=output_partial_obs[,i_test]
#   #input_testing_loc=output[index_testing_loc,]
# 
#   input_data_frame=data.frame(input_loc)
#   testing_input_data_frame=data.frame(output[index_testing_loc,])
#   colnames(testing_input_data_frame)=colnames(input_data_frame)
# 
# 
#   #rf_loc=randomForest(x=(input_loc),y=as.matrix(output_loc))
#   #rf_loc_pred=predict(rf_loc,input_testing_loc)
# 
#   y.lm=lm(as.matrix(output_loc)~.,data=input_data_frame) ##need to make sure the name are the same
#   pred_all=predict(y.lm, testing_input_data_frame,interval="prediction", level = 0.95)  ###okay this looks to work
# 
#   pred_lm_loc_record[,i_test]=pred_all[,1]
# }
# 
# sqrt(mean((pred_lm_loc_record-testing_output)^2))



