## Install and load required packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("limma")


library(tseries)
library(limma)
library(BMS)

## Read-in computer model output (Y) and inputs (X) for both training and test data.
## The "t" indicates test data.
Y<-read.matrix(file="data_humanity_model/Y.txt",header=FALSE)
X<-read.matrix(file="data_humanity_model/X.txt",header=FALSE)
Yt<-read.matrix(file="data_humanity_model/Yt.txt",header=FALSE)
Xt<-read.matrix(file="data_humanity_model/Xt.txt",header=FALSE)

## Pre-compute inner products of output matrices
YY<-t(Y)%*%Y
YYt<-t(Yt)%*%Yt

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

## Formulas for different models. Comment out the ones you don't want.
## intercept model
#formula<-~1 
## linear model
#formula<-~weight+plan+helsp+capacity+engsp+hospG+shelG+foodG+hospC+shelC+foodC+aid+loc	
## maximal model
#formula<-~(weight+plan+helsp+capacity+engsp+hospG+shelG+foodG+hospC+shelC+foodC+aid+loc)^2+I(weight^2)+I(plan^2)+I(helsp^2)+I(capacity^2)+I(engsp^2)+I(hospG^2)+I(shelG^2)+I(foodG^2)+I(hospC^2)+I(shelC^2)+I(foodC^2)
## modal model
formula<-~plan+foodG+foodC+loc+I(plan^2)+foodC:loc

## Computer model matrix and inner product for training data.
H<-model.matrix(formula,data=X)
HH<-t(H)%*%H

## Computer model matrix and inner product for test data.
Ht<-model.matrix(formula,data=Xt)
HHt<-t(Ht)%*%Ht

## Define number of columns of output and model matrices
k<-dim(Y)[2]
m<-dim(H)[2]

## Create an array of squared differences between the elements of the training inputs
X.array<-array(data=rep(0,n*dim(X)[2]*n),dim=c(n,dim(X)[2],n))
for(i in 1:n){
for(j in 1:dim(X)[2]){
X.array[,j,i]<-(X[,j]-X[i,j])^2}}

## Create an array of squared differences between the elements of the test inputs
Xt.array<-array(data=rep(0,n*dim(Xt)[2]*n),dim=c(n,dim(Xt)[2],n))
for(i in 1:n){
for(j in 1:dim(Xt)[2]){
Xt.array[,j,i]<-(Xt[,j]-Xt[i,j])^2}}

## Create an array of differences between the elements of the test inputs multiplied by differences between the elements of the training inputs
T.array<-array(data=rep(0,n*dim(Xt)[2]*n),dim=c(n,dim(Xt)[2],n))
for(i in 1:n){
for(j in 1:n){
T.array[i,,j]<-(X[i,]-Xt[j,])^2
}
  }

## Set up values of prior hyperparameters
invOmega<-0*diag(m)
M<-matrix(0,nrow=m,ncol=k)
S<-0*diag(k)
delta<--k+1

## Function giving the log marginal posterior density of the correlation and nugget parameters (on log scale)
logpi<-function(z){
r<-exp(z)
A<-0*diag(n)
for(j in 1:dim(X)[2]){
A<-A-r[j]*X.array[,j,]}
A<-exp(A)+r[14]*diag(n)
iA<-solve(A)
dA<-determinant(A)$modulus[1]
tHA<-t(H)%*%iA
iOmega<-tHA%*%H
Omega<-solve(iOmega)
M<-Omega%*%tHA%*%Y
S<-t(Y)%*%iA%*%Y-t(M)%*%iOmega%*%M
sum(z)-sum(exp(z))-0.5*k*dA+0.5*k*determinant(Omega)$modulus[1]-0.5*(delta+n+k-1)*determinant(S)$modulus[1]}

## Maximise the log marginal posterior density of the correlation and nugget parameters to find plug-in values
## (Minimum values are specified to avoid numerical problems).
opt<-optim(fn=logpi,par=rep(-5,14),method="L-BFGS-B",control=list(fnscale=-1),lower=rep(-10,14))

## Find plug-in values of correlation parameters and corresponding values of correlation matrix and its inverse
r<-exp(opt$par)
A<-0*diag(n)
for(j in 1:dim(X)[2]){
A<-A-r[j]*X.array[,j,]}
A<-exp(A)+r[14]*diag(n)
iA<-solve(A)

## Calculate updated values of hyperparameters
invOmega.hat<-t(H)%*%iA%*%H+invOmega
Omega.hat<-solve(invOmega.hat)
M.hat<-Omega.hat%*%(t(H)%*%iA%*%Y+invOmega%*%M)
S.hat<-t(Y)%*%iA%*%Y+t(M)%*%invOmega%*%M+S-t(M.hat)%*%invOmega.hat%*%M.hat
delta.hat<-delta+n

## Calculate correlation matrix for test data
At<-0*diag(n)
for(j in 1:dim(Xt)[2]){
At<-At-r[j]*Xt.array[,j,]}
At<-exp(At)+r[14]*diag(n)

## Calculate cross-correlation matrix between test and training data
T<-0*diag(n)
for(j in 1:dim(Xt)[2]){
T<-T-r[j]*T.array[,j,]}
T<-exp(T)

## Compute values of matrix t-distribution
mu<-Ht%*%M.hat+t(T)%*%iA%*%(Y-H%*%M.hat)
R<-At-t(T)%*%iA%*%T+(Ht-t(T)%*%iA%*%H)%*%Omega.hat%*%t(Ht-t(T)%*%iA%*%H)
Q<-S.hat
mmm<-n+delta.hat+k-1

## Calculate RMSE from Table 5
sqrt(mean(((Yt-mu)^2)))

## Calculate RRMSE 
set.seed(1)
rmu<-0*mu
for(i in 1:n){
for(j in 1:k){
sam<-mu[i,j]+sqrt(R[i,i]*S.hat[j,j]/delta.hat)*rt(n=10000,df=delta.hat)
rmu[i,j]<-mean(1/sam)/mean(1/(sam^2))}}

sqrt(mean(((Yt-rmu)^2)/(Yt^2)))

chol_R<-chol(R)
chol_Q<-chol(Q)

B<-t(solve(chol_Q))
A<-solve(chol_R)

## Calculate Cholesky of matrices R and Q (and their inverses)
chol_R<-chol(R)
chol_Q<-chol(Q)

B<-t(solve(chol_Q))
A<-solve(chol_R)

## Calculate matrix of standardised errors E
E<-A%*%(Yt-mu)%*%B

## Draw QQ-plot 
#qqt(sqrt(delta.hat)*as.vector(E),df=delta.hat,pch=16)
#abline(a=0,b=1,lty=2)

## Compute value of omnibus diagnostic
U<-1/det(diag(k)+t(E)%*%E)
U

## Compute coverage matrix and its mean (coverage and average length in Table 5)
cov<-matrix(0,nrow=n,ncol=k)
sum_length=0
for(i in 1:n){
for(j in 1:k){
diff<-qt(0.975,df=delta.hat)*sqrt(R[i,i]*S.hat[j,j]/delta.hat)
low<-mu[i,j]-diff
upp<-mu[i,j]+diff
sum_length=sum_length+upp-low
cov[i,j]<-ifelse(Yt[i,j]>low & Yt[i,j]<upp,1,0)}}

mean(cov)
sum_length/length(Yt)


