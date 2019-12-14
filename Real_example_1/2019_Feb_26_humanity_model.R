library(RobustGaSP)
library(tseries)
library(Rcpp)
library(RcppEigen)
sourceCpp(file='../src/functions.cpp')

Y<-read.matrix(file="data_humanity_model/Y.txt",header=FALSE)
X<-read.matrix(file="data_humanity_model/X.txt",header=FALSE)
Yt<-read.matrix(file="data_humanity_model/Yt.txt",header=FALSE)
Xt<-read.matrix(file="data_humanity_model/Xt.txt",header=FALSE)



k=dim(Y)[2]
## Rows of output matrix
n<-dim(Y)[1]
## Number of rows for each combination of categorical input variables
n0<-30

## Row numbers of output matrix for each combination of categorical input variables
B11<-1:n0
B12<-(n0+1):(2*n0)
B21<-(2*n0+1):(3*n0)
B22<-(3*n0+1):(4*n0)
## Define dummy variables for each categorical input variables
z2<-rep(0,n)
z2[B12]<-1
z2[B22]<-1
z1<-rep(0,n)
z1[B21]<-1
z1[B22]<-1



## Create dataframes for input variables for training and test data.
X<-data.frame(weight=X[,1],plan=X[,2],helsp=X[,3],capacity=X[,4],engsp=X[,5],
              hospG=X[,6],shelG=X[,7],foodG=X[,8],hospC=X[,9],shelC=X[,10],foodC=X[,11],aid=z1,loc=z2)
Xt<-data.frame(weight=Xt[,1],plan=Xt[,2],helsp=Xt[,3],capacity=Xt[,4],engsp=Xt[,5],
               hospG=Xt[,6],shelG=Xt[,7],foodG=Xt[,8],hospC=Xt[,9],shelC=Xt[,10],foodC=Xt[,11],aid=z1,loc=z2)



##plot the data
#par(mfrow=c(3,4))
#for(i_test in 1:11){
# plot(X[,i_test],Y[,4])
#}
###end of the code from the paper

##set up the Ind GP
sum_squares_error=rep(0,k)
predict_mean_ind_GP=matrix(0,dim(Yt)[1],dim(Yt)[2])
predict_ind_GP_lower_95=matrix(0,dim(Yt)[1],dim(Yt)[2])
predict_ind_GP_upper_95=matrix(0,dim(Yt)[1],dim(Yt)[2])

kernel_type_here="matern_5_2"

#H=matrix(1,n,1)
#H_testing=matrix(1,num_testing,1)

#H=cbind(matrix(1,n,1),X[,c(2,8,11,12,13)],X[,2]^2,X[,c(8,11)]*X[,12],X[,11]*X[,13],X[,12]*X[,13])
#H_testing=cbind(matrix(1,num_testing,1),Xt[,c(2,8,11,12,13)],Xt[,2]^2,Xt[,c(8,11)]*Xt[,12],Xt[,11]*Xt[,13],Xt[,12]*Xt[,13])

num_testing=dim(Xt)[1]

##here the mean is the intercept and the 11th one
H=cbind(matrix(1,n,1),X[,11])
H_testing=cbind(matrix(1,num_testing,1),Xt[,11])

H=as.matrix(H)
H_testing=as.matrix(H_testing)
for(i_k in 1:k){
  #m=rgasp(design=X,response=Y[,i_k],nugget.est= TRUE,kernel_type='pow_exp',alpha=rep(2,dim(X)[2]))
  m=rgasp(design=X,response=Y[,i_k],trend=H,nugget.est= TRUE,kernel_type=kernel_type_here,alpha=rep(2,dim(X)[2]))
  
  m_pred=predict(m,Xt,testing_trend=H_testing)
  predict_mean_ind_GP[,i_k]=m_pred$mean
  predict_ind_GP_lower_95[,i_k]=m_pred$lower95
  predict_ind_GP_upper_95[,i_k]=m_pred$upper95
  
  #  sum_squares_error[i_k]= sum( (m_pred$mean-Yt[,i_k])^2)
  
}


sqrt(mean((predict_mean_ind_GP-Yt)^2))

length(which(predict_ind_GP_lower_95<Yt & predict_ind_GP_upper_95>Yt ))/length(Yt)
mean(abs(predict_ind_GP_upper_95-predict_ind_GP_lower_95))



###############################


###simulation code to explore if one subtract and add back the mean structure 
##and deal with it in a fully probabilistic way (marginalizing out, etc.)

###the optimization function and detivative 
F_funct<-function(X,G){##p \times k matrix matrix X, G is a list for k items
  return_val=0
  for(i in 1: dim(X)[2]){
    return_val=return_val+t(X[,i])%*%G[[i]]%*%X[,i]
  }
  -return_val
}


##this is the Y'(0)
F_funct_dev<-function(X,G){##p \times k matrix matrix X, G is a list for k items
  return_matrix=matrix(0,dim(X)[1],dim(X)[2])
  for(i in 1: dim(X)[2]){
    return_matrix[,i]=G[[i]]%*%X[,i]
  }
  -2*return_matrix
}


neg_log_lik_diff_cov_with_trend<-function(param,A_ini){
  
  #print(param)
  
  #beta=exp(param[1:d])
  #tau=exp(param[(d+1):(2*d) ])
  G=as.list(1:d)
  sum_det=0
  
  for(i_d in 1:d){
    #beta=exp(param[ ((i_d-1)*p+1):(i_d*p)])
    #tau=exp(param[d*p+i_d])
    
    ##same range parameter but different variance
    beta=exp(param[ 1:p])
    tau=exp(param[p+i_d])
    
    #tau=exp(param[p+1])
    
    
    Sigma=tau*separable_kernel(R0_list, beta, kernel_type=kernel_type_here, alpha=rep(2,p) )

    Sigma_tilde=Sigma+diag(n)
    
    
    L=t(chol(Sigma_tilde))
    ##one should avoid solve here
    L_H=t(chol(t(H)%*%solve(Sigma_tilde)%*%H))
    
    ##one may use cholesky here, but I just inversion for now
    G[[i_d]]=output%*%M%*%solve(Sigma%*%M+diag(n))%*%Sigma%*%M%*%t(output)
    sum_det=sum_det+(sum(log(diag(L)))+sum(log(diag(L_H))))
    
  }
  A_est_here=Optimization_Stiefel_Manifold(A_ini, G=G,max_iter=100)
  
  #print((-sum_det -(k*(n-q))/2*log(tr_output_M_output+F_funct(A_est_here,G))))
  -(-sum_det -(k*(n-q))/2*log(tr_output_M_output+F_funct(A_est_here,G)))
  
  
}

Get_A_est_sigma_0_2_with_trend<-function(param,A_ini){
  G=as.list(1:d)
  sum_det=0
  
  for(i_d in 1:d){
    #beta=exp(param[ ((i_d-1)*p+1):(i_d*p)])
    beta=exp(param[ 1:p])
    
    tau=exp(param[p+i_d])
    
    ##same range parameter but different variance
    #beta=exp(param[ 1:p])
    #tau=exp(param[p+i_d])
    
    
    
    Sigma=tau*separable_kernel(R0_list, beta, kernel_type=kernel_type_here, alpha=rep(2,p) )
    
    Sigma_tilde=Sigma+diag(n)
    
    
    L=t(chol(Sigma_tilde))
    ##one should avoid solve here
    L_H=t(chol(t(H)%*%solve(Sigma_tilde)%*%H))
    
    ##one may use cholesky here, but I just inversion for now
    G[[i_d]]=output%*%M%*%solve(Sigma%*%M+diag(n))%*%Sigma%*%M%*%t(output)
    sum_det=sum_det+(sum(log(diag(L)))+sum(log(diag(L_H))))
    
  }
  A_est_here=Optimization_Stiefel_Manifold(A_ini, G=G,max_iter=100)
  #A_est_here
  return_list=list(1:2)
  return_list[[1]]=A_est_here
  return_list[[2]]=(tr_output_M_output+F_funct(A_est_here,G))/(k*(n-q))
  return_list
}


pred_mean_diff_cov_with_trend<-function(param,sigma_2_0_here,A_hat){
  #beta=exp(param[1:d])
  #tau=exp(param[(d+1):(2*d)])
  
  Z_hat=matrix(0,d,n)
  z_star=matrix(0,d,num_testing)
  D_matrix=matrix(0,d,num_testing)
  for(i_d in 1:d){
    beta=exp(param[ 1:p])
    tau=exp(param[p+i_d])
    sigma_2_0_here=as.numeric(sigma_2_0_here)
    sigma_2=(tau*sigma_2_0_here)

    
    Sigma=sigma_2*separable_kernel(R0_list, beta, kernel_type=kernel_type_here, alpha=rep(2,p) )
    Sigma_star=sigma_2*separable_kernel(r0_list, beta, kernel_type=kernel_type_here, alpha=rep(2,p) )

    middle_part=t(A_hat[,i_d])%*%output%*%M%*%solve(Sigma%*%M+sigma_2_0_here*diag(n))
    Z_hat[i_d,]=middle_part%*%Sigma
    z_star[i_d,]=t(Sigma_star)%*%t(middle_part)

    ##calculating the variance
    Sigma_tilde=Sigma+sigma_2_0_here*diag(n)
    L=t(chol(Sigma_tilde))
    Sigma_tilde_inv_star=backsolve(t(L),forwardsolve(L,Sigma_star))            
    #R.inv_X=backsolve(t(L),forwardsolve(L,X))            
    
    ## In the future, I will use Cholesky here
    H_t_inv_Sigma_tilde_H_inv=solve((t(H)%*%solve(Sigma_tilde)%*%H))
    
    
    Sigma_inv_H=backsolve(t(L),forwardsolve(L,H))            
    D_star=rep(0,num_testing)
    for(i_testing in 1:num_testing){
      left_term=H_testing[i_testing,]-t(Sigma_star[,i_testing])%*%Sigma_inv_H
      diff=left_term%*%H_t_inv_Sigma_tilde_H_inv%*%t(left_term)
      D_star[i_testing]=-t(Sigma_star[,i_testing])%*%Sigma_tilde_inv_star[,i_testing]+diff
    }
    
    D_star=D_star+sigma_2+sigma_2_0_here
    
    D_matrix[i_d,]=D_star
  }
  var_hat=matrix(0,k,num_testing)
  
  for(i_testing in 1:num_testing){
    var_hat[,i_testing]=diag(A_hat%*%diag(D_matrix[,i_testing])%*%t(A_hat))+diag(sigma_2_0_here*as.numeric(1+t(H_testing[i_testing,])%*%HH_inv%*%(H_testing[i_testing,]))*(diag(k)- A_hat%*%t(A_hat)))
  }
  #pred_mean=(output-A_hat%*% Z_hat)%*%H%*%solve(t(H)%*%H)%*%H_testing+A_hat%*%Z_hat
  pred_mean=(output-A_hat%*% Z_hat)%*%H%*%solve(t(H)%*%H)%*%t(H_testing)+A_hat%*%z_star
  
  pred_list=as.list(1:2)
  pred_list[[1]]=pred_mean
  pred_list[[2]]=var_hat
  return(pred_list)
}
##shared cov
get_chol<-function(x,beta){
  R0_00=abs(outer(x,x,'-'))
  R=matern_5_2_funct(R0_00,beta)
  #R=pow_exp_funct(R0_00,beta,alpha=1)
  
  rcppeigen_get_chol(R)
}





###GPPCA
input=as.matrix(X)
p=dim(X)[2]

HH_inv=solve(t(H)%*%H)
HH=H%*%HH_inv%*%t(H)
M=diag(n)-HH
R0_list=list(1:p)
for(i_p in 1:p){
  R0_list[[i_p]]=abs(outer(input[,i_p],input[,i_p],'-'))
}
q=1

output=t(Y)
output_sub_mean=matrix(0,k,n)
d=k

for(i_k in 1:k){
  output_sub_mean[i_k,]=output[i_k,]-mean(output[i_k,])
}


##
svd_output_sub_mean=svd(output_sub_mean)
A_ini=svd_output_sub_mean$u[,1:d] 
tr_output_M_output=sum(diag(output%*%M%*%t(output)))

##initial parameters
param_ini=c(rep(-3,p),rep(8,d))

##takes 7 min or so to optimize
system.time(
  for(i in 1:1){
m_trend=try(optim(param_ini,neg_log_lik_diff_cov_with_trend,A_ini=A_ini,control = list(maxit = 100), 
                    method="L-BFGS-B"),silent=T)
  }
)


A_sigma_2_0=Get_A_est_sigma_0_2_with_trend(m_trend$par,A_ini=A_ini)

A_with_mean_est=A_sigma_2_0[[1]]
sigma_2_0_est=A_sigma_2_0[[2]]


testing_input=as.matrix(Xt)
r0_list=list(1:p)
for(i_p in 1:p){
  r0_list[[i_p]]=abs(outer(input[,i_p], testing_input[,i_p],'-'))
}


  
pred_GPPCA_with_mean_all=pred_mean_diff_cov_with_trend(param=m_trend$par,sigma_2_0_here=sigma_2_0_est, A_hat=A_with_mean_est)


LB95_GPPCA=pred_GPPCA_with_mean_all[[1]]+sqrt(pred_GPPCA_with_mean_all[[2]])*qnorm(0.025)
UB95_GPPCA=pred_GPPCA_with_mean_all[[1]]+sqrt(pred_GPPCA_with_mean_all[[2]])*qnorm(0.975)


sqrt(mean((t(pred_GPPCA_with_mean_all[[1]])-Yt)^2))

length(which( t(LB95_GPPCA)<Yt& t(UB95_GPPCA)>Yt))/length(Yt)
mean(abs(UB95_GPPCA-LB95_GPPCA))



#sqrt(mean((mean(Y)-Yt)^2))

sigma_hat=as.numeric(sigma_2_0_est)*exp(m_trend$par[(p+1):(p+d)])

Cov_est_data=A_with_mean_est%*%diag(sigma_hat)%*%t(A_with_mean_est)+as.numeric(sigma_2_0_est)*diag(k)


library(ggplot2)
library(reshape2)
rng=c(min(Cov_est_data),max(Cov_est_data))

#pdf(file='Cov_est_data.pdf',height=4,width=8)

qplot(x=Var1+1, y=Var2+1, data=melt(Cov_est_data), fill=value, xlab='Day',ylab='Day',geom="tile")+
  scale_fill_gradient2(limits=c((rng[1]), (rng[2])),low = "royalblue4", high = "royalblue1") +
  labs(fill = "Covariance")+theme(text = element_text(size=15),axis.text = element_text(hjust=1.5))
#dev.off()

#pdf(file='pred_two_days.pdf',height=5,width=7)

plot_i=4
plot(Yt[,plot_i],type='p',col='black',pch=1,cex=1.2,xlab='Held out runs',ylab='Output',ylim=c(0,30000))
lines((predict_mean_ind_GP)[,plot_i],col='red',pch=20,type='p',cex=1.2)
lines(t(pred_GPPCA_with_mean_all[[1]])[,plot_i],col='blue',pch=20,type='p',cex=1)

plot_j=5
lines(Yt[,plot_j],type='p',col='black',pch=0,cex=1.2,xlab='Held out runs',ylab='Output')
lines((predict_mean_ind_GP)[,plot_j],col='red',pch=15,type='p',cex=1.2)
lines(t(pred_GPPCA_with_mean_all[[1]])[,plot_j],col='blue',pch=15,type='p',cex=1)

legend("topright", legend=c("Casualties on day 5", "GPPCA prediction for day 5", "Ind GP prediction for day 5",
                           "Casualties on day 6", "GPPCA prediction for day 6", "Ind GP prediction for day 6"),
       col=rep(c('black','red','blue'),2), pch=c(c(1,20,20),c(0,15,15)),ncol = 2,cex=0.9)
#dev.off()
##Gaussian kernel with constant mean
# ##ind GP
# > sqrt(mean((predict_mean_ind_GP-Yt)^2))
# [1] 364.304
# > length(which(predict_ind_GP_lower_95<Yt & predict_ind_GP_upper_95>Yt ))/length(Yt)
# [1] 0.9183333
# > mean(abs(predict_ind_GP_upper_95-predict_ind_GP_lower_95))
# [1] 1176.64
# ##GPPCA
# > sqrt(mean((t(pred_GPPCA_with_mean_all[[1]])-Yt)^2))
# [1] 333.8055
# > length(which( t(LB95_GPPCA)<Yt& t(UB95_GPPCA)>Yt))/length(Yt)
# [1] 0.9483333
# > mean(abs(UB95_GPPCA-LB95_GPPCA))
# [1] 1516.124
# 
# ##Gaussian kernel with the 11th parameter as the mean
# ##ind
# > sqrt(mean((predict_mean_ind_GP-Yt)^2))
# [1] 404.0626
# > length(which(predict_ind_GP_lower_95<Yt & predict_ind_GP_upper_95>Yt ))/length(Yt)
# [1] 0.9183333
# > mean(abs(predict_ind_GP_upper_95-predict_ind_GP_lower_95))
# [1] 1165.86
# ##GPPCA
# > sqrt(mean((t(pred_GPPCA_with_mean_all[[1]])-Yt)^2))
# [1] 318.4018
# > length(which( t(LB95_GPPCA)<Yt& t(UB95_GPPCA)>Yt))/length(Yt)
# [1] 0.9566667
# > mean(abs(UB95_GPPCA-LB95_GPPCA))
# [1] 1309.283
# 
# 
# ##matern 2.5 with constant mean
# ##ind GP
# > sqrt(mean((predict_mean_ind_GP-Yt)^2))
# [1] 339.5839
# > length(which(predict_ind_GP_lower_95<Yt & predict_ind_GP_upper_95>Yt ))/length(Yt)
# [1] 0.93
# > mean(abs(predict_ind_GP_upper_95-predict_ind_GP_lower_95))
# [1] 984.2133
# #GPPCA
# > sqrt(mean((t(pred_GPPCA_with_mean_all[[1]])-Yt)^2))
# [1] 282.3827
# > length(which( t(LB95_GPPCA)<Yt& t(UB95_GPPCA)>Yt))/length(Yt)
# [1] 0.9616667
# > mean(abs(UB95_GPPCA-LB95_GPPCA))
# [1] 1220.965
# 
# ##matern 2.5 with selected mean
# ##ind
# > sqrt(mean((predict_mean_ind_GP-Yt)^2))
# [1] 331.5376
# > length(which(predict_ind_GP_lower_95<Yt & predict_ind_GP_upper_95>Yt ))/length(Yt)
# [1] 0.9266667
# > mean(abs(predict_ind_GP_upper_95-predict_ind_GP_lower_95))
# [1] 966.724
# ##GPPCA
# > sqrt(mean((t(pred_GPPCA_with_mean_all[[1]])-Yt)^2))
# [1] 273.6336
# > length(which( t(LB95_GPPCA)<Yt& t(UB95_GPPCA)>Yt))/length(Yt)
# [1] 0.9566667
# > mean(abs(UB95_GPPCA-LB95_GPPCA))
# [1] 1175.54




