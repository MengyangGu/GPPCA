library(Rcpp)
library(RcppEigen)
library(RobustGaSP)
library(rstiefel)
library(maps)#Install maps package if not done before
#library(pcaMethods)
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
  L_H=t(chol(t(H)%*%solve(Sigma_tilde)%*%H))
  sum_det=sum(log(diag(L)))+sum(log(diag(L_H)))
  

  neg_lik=-(-d*sum_det-(tr_output_M_output-sum(svd_YL$d[1:d]^2))/(2*est_measurement_var) )
  
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
  
  #neg_lik=-(-sum_det -(k*(num_obs-q))/2*log(tr_output_M_output+F_Funct_Large_k(A_est_here,UD)))
}


pred_mean_share_cov_with_trend<-function(param,sigma_2_0_here,A_hat,output_here,kernel_type){
  beta=exp(param[1])
  tau=exp(param[2])
  
  Z_hat=matrix(0,d,num_obs)
  z_star=matrix(0,d,num_obs)
  D_matrix=matrix(0,d,num_obs)
  
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
  D_star=rep(0,(num_obs)) 
  for(i_test in 1:(num_obs)){
    D_star[i_test]=-t(Sigma_star[,i_test])%*%Sigma_tilde_inv%*%Sigma_star[,i_test]+t(h_t_star_minus_H_t_Sigma_inv_Sigma[,i_test])%*%(H_t_inv_Sigma_tilde_H_inv)%*%h_t_star_minus_H_t_Sigma_inv_Sigma[,i_test]
  }
  
  D_star=D_star+sigma_2+sigma_2_0_here
  
  pred_mean=(output_here-A_hat%*% Z_hat)%*%H%*%solve(t(H)%*%H)%*%t(H_testing)+A_hat%*%t(z_star)
  
  Theta_est=(output_here-A_hat%*% Z_hat)%*%H%*%solve(t(H)%*%H)
  pred_list=as.list(1:3)
  pred_list[[1]]=pred_mean
  pred_list[[2]]=D_star
  pred_list[[3]]=Theta_est
  
  return(pred_list)
}



gridded_data=readRDS("NOAA_gridded_temperature/NOAA_gridded_data_1999_2018.rds")
output_all=gridded_data$temperature

grid_with_full_obs=gridded_data$loc_full_obs
time=gridded_data$months
Lat= seq(-87.5, 87.5, length=36)
Lon=seq(2.5, 357.5, length=72)


set.seed(1)  
#num_obs_all=dim(output_all)[2]
k=dim(output_all)[1]
input=as.numeric(1:(dim(output_all)[2]))
output=output_all

  
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
L_M=eigen_M$vector%*%diag(sqrt(abs(eigen_M$values) ))[,1:(num_obs-q)]
output_L_M=output%*%L_M

tr_output_M_output=sum(diag(t(output_L_M)%*%output_L_M))

###this value is the measurement var estimated from the station
est_measurement_var=0.1
d=100


R0_00=(abs(outer(input,input,'-')))

param_ini=c(0,10)



kernel_type='matern_5_2'


#####maximum posterior mode
#neg_log_post_shared_cov(param_ini,kernel_type)

#m_trend_post=try(optim(param_ini,neg_log_post_shared_cov,kernel_type=kernel_type,lower=rep(log(10^{-10}),2),
#                       upper=rep(log(10^{10}),2),control = list(maxit = 100), method="L-BFGS-B"),silent=T)

##marximum likelihood
m_lik=try(optim(param_ini,neg_log_lik_with_trend_shared_cov,kernel_type=kernel_type,lower=rep(log(10^{-10}),2),
                       upper=rep(log(10^{10}),2),control = list(maxit = 100), method="L-BFGS-B"),silent=T)

Get_A_est_all=Get_A_est_sigma_0_2_with_trend_large_k(m_lik$par,kernel_type=kernel_type)
A_with_mean_est=Get_A_est_all[[1]]

##estimated variance of the noise
#sigma_2_0_est=Get_A_est_all[[2]]
##fixed variance of the noise
sigma_2_0_est=est_measurement_var

r0=abs(outer(input, input,'-'))

H_testing=matrix(1,num_obs,q)
H_testing[,2]=input
H_testing=as.matrix(H_testing)

pred_GPPCA_with_mean_all=pred_mean_share_cov_with_trend(param=m_lik$par,sigma_2_0_here=sigma_2_0_est, A_hat=A_with_mean_est,output_here=output,kernel_type=kernel_type)


##theta est
theta_1_hat=pred_GPPCA_with_mean_all[[3]][,1]
theta_2_hat=pred_GPPCA_with_mean_all[[3]][,2]






mapmat=matrix(NA,72,36)
mapmat[grid_with_full_obs ]=theta_1_hat
max(mapmat[grid_with_full_obs])
min(mapmat[grid_with_full_obs])
int=seq(-1.6,1.6,length.out=81)
#pdf(file='estimated_theta_1.pdf',height=6.2,width=6.2)
rgb.palette=colorRampPalette(c('black','blue', 'darkgreen','green', 'yellow','pink','red','maroon'),interpolate='spline')
filled.contour(Lon, Lat, mapmat, color.palette=rgb.palette, levels=int,
               plot.title=title(main="Estimated intercept",
                                xlab="Longitude",ylab="Latitude", cex.lab=1.5),
               plot.axes={axis(1, cex.axis=1.5);axis(2, cex.axis=1.5);map('world2', add=TRUE);grid()},
               key.title=title(main=~degree*C),
               key.axes={axis(4, cex.axis=1)})
#dev.off()


mapmat=matrix(NA,72,36)
mapmat[grid_with_full_obs ]=theta_2_hat
max(mapmat[grid_with_full_obs])
min(mapmat[grid_with_full_obs])
int=seq(-0.015,0.015,length.out=81)
#pdf(file='estimated_theta_2.pdf',height=6.2,width=6.2)
rgb.palette=colorRampPalette(c('black','blue', 'darkgreen','green', 'yellow','pink','red','maroon'),interpolate='spline')
filled.contour(Lon, Lat, mapmat, color.palette=rgb.palette, levels=int,
               plot.title=title(main="Estimated monthly temperature change rate",
                                xlab="Longitude",ylab="Latitude", cex.lab=1.5),
               plot.axes={axis(1, cex.axis=1.5);axis(2, cex.axis=1.5);map('world2', add=TRUE);grid()},
               key.title=title(main=~degree*C),
               key.axes={axis(4, cex.axis=1)})
#dev.off()

##plot some of them
i_plot=which(grid_with_full_obs==1857)

plot(output[i_plot,],ylim=c(-3,5))
lines(theta_1_hat[i_plot]+theta_2_hat[i_plot]*(1:240))

i_plot=which(grid_with_full_obs==1776)
plot(output[i_plot,],ylim=c(-3,5))
lines(theta_1_hat[i_plot]+theta_2_hat[i_plot]*(1:240))

which(grid_with_full_obs==1857)
##

((Lon)+180)%%360-180

Lat[26]
Lon[57]

##this seems roughly New york
25* 72+57

#this is roughly Los Angeles
Lat[25]
Lat[48]

24*72+48


##new york
i_plot=which(grid_with_full_obs==1857)
plot(output[i_plot,],ylim=c(-3,5))
lines(pred_GPPCA_with_mean_all[[1]][i_plot,],type='l',col='blue')
lines(theta_1_hat[i_plot]+theta_2_hat[i_plot]*(1:240))



##los angeles
i_plot=which(grid_with_full_obs==1776)
plot(output[i_plot,],ylim=c(-3,5))
lines(pred_GPPCA_with_mean_all[[1]][i_plot,],type='l',col='blue')
lines(theta_1_hat[i_plot]+theta_2_hat[i_plot]*(1:240))
