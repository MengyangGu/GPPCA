library(RobustGaSP)
library(rospca)
library(Rcpp)
library(RcppEigen)
library(nloptr)
library(FastGP)
library(rstiefel)

sourceCpp(file='../src/functions.cpp') 

####


neg_log_lik_shared_cov_FFBS<-function(param,kernel_type){

  G_log_det_cov=Get_G_log_det_cov(param, output, delta_x,d=1,kernel_type = kernel_type)
  G=output_2-G_log_det_cov[[1]]
  
  eigen_G=eigen(G)
  
  -(-(sum(G_log_det_cov[[2]]))*d/2-(n*k)/2*log(trace_output_2-sum(eigen_G$values[1:d]) ))
  
  
}

posterior_AZ_funct<-function(param,A_hat){
  
  beta=param[1]
  sigma_2=param[2]
  sigma_2_0=param[3]
  
  R0_00=abs(outer(input,input,'-'))
  R=matern_5_2_funct(R0_00,beta)
  Sigma=sigma_2*R
  L_tilde=t(chol(Sigma+sigma_2_0*diag(n)))
  

  weighted_output_t=(t(output)%*%A_hat)

  
  post_mean=A_hat%*%t(Sigma%*%backsolve(t(L_tilde),forwardsolve(L_tilde,weighted_output_t)) ) 
  
  D=Sigma-Sigma%*%backsolve(t(L_tilde),forwardsolve(L_tilde,Sigma))  
  
  post_var=A_hat^2%*%t(as.matrix(diag(D),n,1)) ##this only for d=1

  return_list=as.list(1:2)
  return_list[[1]]=post_mean
  return_list[[2]]=post_var
  return_list
}

 


n=100
k=2
d=1

beta_real=.01

sigma_real=1
sigma_0_real=sqrt(1)

input=as.numeric(seq(1,n,(n-1)/(n-1)))

R0=abs(outer(input,input,'-'))
R=matern_5_2_funct(R0,beta_real)
L_sample=t(chol(R))



set.seed(1)

input=sort(input)
delta_x=input[2:length(input)]-input[1:(length(input)-1)]

A=rustiefel(k, d)  ##sample from Stiefel manifold

Factor=matrix(0,d,n)

for(i in 1: d){
  Factor[i,]=sigma_real^2*L_sample%*%rnorm(n)
}

output=A%*%Factor+matrix(rnorm(n*k,mean=0,sd=sigma_0_real),k,n)


plot(output[2,])
output_2=output%*%t(output)

trace_output_2=sum(output^2)

svd_output=svd(output)

A_est_pc=svd_output$u[,1:d]

####
A_ini=svd_output$u[,1:d] 

param_ini=c(log(.1),log(10))

kernel_type='matern_5_2'
m=try(optim(param_ini,neg_log_lik_shared_cov_FFBS, kernel_type=kernel_type, method="L-BFGS-B"),silent=T)
m

beta_hat=exp(m$par[1])
tau_hat=exp(m$par[2])

G_log_det_cov_hat=Get_G_log_det_cov( m$par, output, delta_x,d=1,kernel_type)

G_hat=output_2- G_log_det_cov_hat[[1]]

eigen_G_hat=eigen(G_hat)

sigma_2_0_hat=(trace_output_2-sum(eigen_G_hat$values[1:d]) )/(n*k)

plot(output[1,],output[2,],col='black',pch=20,xlab=expression(y[1]),ylab=expression(y[2]))
abline(a=0,b=eigen_G_hat$vectors[2]/eigen_G_hat$vectors[1],col='blue',lty=2)
abline(a=0,b=A_est_pc[2]/A_est_pc[1],col='red',lty=3)
abline(a=0,b=A[2]/A[1],col='black',lty=1)




rospca::angle(A,A_est_pc)
rospca::angle(A,eigen_G_hat$vectors[,1])


R0_00=abs(outer(input,input,'-'))
R=matern_5_2_funct(R0_00,beta_hat)

##cholesky
#tilde_Sigma_hat=tau_hat*R
#LL=t(chol(diag(n)-solve(tilde_Sigma_hat+diag(n))))
#middle_term=solve(solve(tilde_Sigma_hat)+diag(n))
#output_trans=output%*%as.matrix(t(chol(middle_term)))

eigen_middle_matrix=eigen(diag(n)-solve(tau_hat*R+diag(n)))

output_transformation=output%*%eigen_middle_matrix$vectors%*%diag(eigen_middle_matrix$values)




cor(output[1,],output[2,])
cor(output_transformation[1,],output_transformation[2,])

#pdf(file='demonstration_sigma_0_1.pdf',height=4,width=10)
par(mfrow=c(1,3))
plot(output[1,],pch=1,ylim=c(min(output)-.5,max(output)+.5),cex=1.5,xlab='x',ylab='y',cex.lab=1.7,cex.axis=1.7, mgp=c(2.5,1,0))
lines(output[2,],cex=1.5,pch=20,type='p')
plot(output_transformation[1,],ylim=c(min(output_transformation),max(output_transformation)),pch=1,cex=1.5,xlab='x',ylab='',cex.lab=1.7,cex.axis=1.7)
title(ylab=expression(tilde(y)), mgp=c(2.2,1,0),cex.lab=2)
lines(output_transformation[2,],cex=1.5,pch=20,type='p')
plot(output[1,],output[2,],ylim=c(min(output),max(output)),col='black',pch=17,cex=1.5,xlab=expression(y[1]),ylab=expression(y[2]),cex.lab=1.7, cex.axis=1.7,mgp=c(2.5,1,0))
abline(a=0,b=eigen_G_hat$vectors[2]/eigen_G_hat$vectors[1],cex=1.5,col='blue',lty=2)
abline(a=0,b=A_est_pc[2]/A_est_pc[1],col='red',cex=1.5,lty=3)
abline(a=0,b=A[2]/A[1],col='black',cex=1.5,lty=1)
#dev.off()






AZ_posterior=posterior_AZ_funct(param=c(beta_hat,sigma_2_0_hat*tau_hat,sigma_2_0_hat),A_hat=eigen_G_hat$vectors[,1:d])

Y_hat=AZ_posterior[[1]]
Y_var=AZ_posterior[[2]]
Y_95_lower=Y_hat+sqrt(Y_var)*qnorm(0.025)
Y_95_upper=Y_hat+sqrt(Y_var)*qnorm(0.975)

pdf(file='demonstration_est_AZ_sigma_0_1.pdf',height=4,width=8)
par(mfrow=c(1,2))
index_d=1
plot(input, (A%*%Factor)[index_d,],type='l',ylim=c(min(output[index_d,]),max(output[index_d,])),cex=1.5,cex.lab=1,cex.axis=1,xlab='x',ylab=expression(hat(Y)[1]), mgp=c(2,1,0))
polygon(c(input,rev(input)),c(Y_95_lower[index_d,],rev(Y_95_upper[index_d,])),col = "grey80", border = F)
lines(input, (A%*%Factor)[index_d,],type='l',ylim=c(min(output[index_d,]),max(output[index_d,])),cex=1.5)
lines(input,Y_hat[index_d,],col='blue',type='l',lty=2,cex=1.5)
lines(input,(A_est_pc%*%t(A_est_pc)%*%output)[index_d,],type='l',col='red',lty=3,cex=1.5)
legend('bottomleft',legend=c('PCA', 'GPPCA','Truth'),lty=c(3,2,1),col=c('red','blue','black'),cex=.6)
lines(input, output[index_d,], type='p',pch=20,cex=.8)
#sum((A%*%Factor>Y_95_lower)& (A%*%Factor<Y_95_upper))/(n*k)
index_d=2
plot(input, (A%*%Factor)[index_d,],type='l',ylim=c(min(output[index_d,]),max(output[index_d,])),cex=1.5,cex.lab=1,cex.axis=1,xlab='x',ylab=expression(hat(Y)[2]), mgp=c(2,1,0))
polygon(c(input,rev(input)),c(Y_95_lower[index_d,],rev(Y_95_upper[index_d,])),col = "grey80", border = F)
lines(input, (A%*%Factor)[index_d,],type='l',ylim=c(min(output[index_d,]),max(output[index_d,])),cex=1.5)
lines(input,Y_hat[index_d,],col='blue',type='l',lty=2,cex=1.5)
lines(input,(A_est_pc%*%t(A_est_pc)%*%output)[index_d,],type='l',col='red',lty=3,cex=1.5)
legend('bottomright',legend=c('PCA', 'GPPCA','Truth'),lty=c(3,2,1),col=c('red','blue','black'),cex=.6)
lines(input, output[index_d,], type='p',pch=20,cex=.8)
#lines(output[index_d,],type='p',pch=20,cex=0.8)
dev.off()



##sigma_0=0.1
#> cor(output[1,],output[2,])
#[1] -0.8313914
#> cor(output_transformation[1,],output_transformation[2,])
#[1] -0.9993769


##sigma_0=1

#> cor(output[1,],output[2,])
#[1] -0.1841232
#> cor(output_transformation[1,],output_transformation[2,])
#[1] -0.9921484

