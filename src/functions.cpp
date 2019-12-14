
// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*- 
 
// we only include RcppEigen.h which pulls Rcpp.h in for us 
#include <iostream> 
#include <RcppEigen.h> 
#include <Rcpp.h> 
#include <cmath> 
// [[Rcpp::depends(RcppEigen)]] 

using namespace Rcpp;
using namespace std;
using namespace Eigen; 
 
 

//July 16, 2018
////Construct_W_matern_5_20_matern_5_2 
// [[Rcpp::export]] 
MatrixXd Construct_W0_matern_5_2(const double sigma2, const double lambda){ 
  //int num_dim=sigma2.size(); 
  
  Eigen::MatrixXd W0= Eigen::MatrixXd::Zero(3,3); 
  //Eigen::MatrixXd d= Eigen::MatrixXd::Zero(3,3); //the first row has all zeros  
  
  W0(0,0)=sigma2; 
  W0(0,2)=W0(2,0)=-sigma2*pow(lambda,2.0)/3.0; 
  W0(1,1)=sigma2*pow(lambda,2.0)/3.0; 
  W0(2,2)=sigma2*pow(lambda,4.0); 

  return W0; 
} 
////Construct_W0_exp 
// [[Rcpp::export]] 
MatrixXd Construct_W0_exp(const double sigma2, const double lambda){ 
  //int num_dim=sigma2.size(); 
  
  Eigen::MatrixXd W0= Eigen::MatrixXd::Zero(1,1); 
  
  W0(0,0)=sigma2; 
  
  return W0; 
} 

 
 ////Construct_G_matern_5_2 
 // [[Rcpp::export]] 
 List Construct_G_matern_5_2(Eigen::VectorXd delta_x, double lambda){  //be careful about delta_x,lambda if one only wants to sample one 
   int num_obs=delta_x.size()+1; 
   //int num_dim=lambda.size();  
   List GG(num_obs);  
   GG[0]=Eigen::MatrixXd::Zero(3,3); 
   
   Eigen::MatrixXd d= Eigen::MatrixXd::Zero(3,3); //the first row has all zeros  
   
   // num_dim list, each is 3(num_obs)\times 3 list 
  // for(int i_GG=0;i_GG<num_dim;i_GG++){ 
  //   Eigen::MatrixXd d= Eigen::MatrixXd::Zero(num_obs,9);  //the first row has all zeros  
     for(int j_GG=0;j_GG<(num_obs-1);j_GG++){ 
       int j_GG_1=j_GG+1;    
       d(0,0)=pow(lambda,2.0)*pow(delta_x[j_GG],2.0)+2*lambda*delta_x[j_GG]+2; 
       d(1,0)=-pow(lambda,3.0)*pow(delta_x[j_GG],2.0); 
       d(2,0)=pow(lambda,4.0)*pow(delta_x[j_GG],2.0)-2*pow(lambda,3.0)*delta_x[j_GG]; 
       d(0,1)=2*(lambda*pow(delta_x[j_GG],2.0)+delta_x[j_GG]); 
       d(1,1)=-2*(pow(lambda,2.0)*pow(delta_x[j_GG],2.0)-lambda*delta_x[j_GG]-1); 
       d(2,1)=2*(pow(lambda,3.0)*pow(delta_x[j_GG],2.0)-3*pow(lambda,2.0)*delta_x[j_GG]); 
       d(0,2)=pow(delta_x[j_GG],2); 
       d(1,2)=2*delta_x[j_GG]-lambda*pow(delta_x[j_GG],2.0); 
       d(2,2)=pow(lambda,2.0)*pow(delta_x[j_GG],2.0)-4*lambda*delta_x[j_GG]+2;     
       d=exp(-lambda*delta_x[j_GG])/2.0*d;
       GG[j_GG_1]=d; 
     } 
    // GG[i_GG]=d; 
   //} 
   return GG; 
 } 



////Construct_G_exp
// [[Rcpp::export]] 
List Construct_G_exp(Eigen::VectorXd delta_x, double lambda){  //be careful about delta_x,lambda if one only wants to sample one 
  int num_obs=delta_x.size()+1; 
  //int num_dim=lambda.size();  
  List GG(num_obs);  
  GG[0]=Eigen::MatrixXd::Zero(1,1); 
  Eigen::MatrixXd d= Eigen::MatrixXd::Zero(1,1); 
  
  for(int j_GG=0;j_GG<(num_obs-1);j_GG++){ 
    d(0,0)=exp(-delta_x[j_GG]*lambda);
    GG[j_GG+1]=d; 
  }
  
  return GG;
}


////Construct_W_matern_5_2  
// [[Rcpp::export]] 
List Construct_W_matern_5_2(double sigma2,Eigen::VectorXd delta_x, double lambda, MatrixXd W0){  //be careful about delta_x,lambda if one only wants to sample one 
  int num_obs=delta_x.size()+1; 
  //int num_dim=sigma2.size();  
  List Wi(num_obs);  
  Wi[0]=W0; 
  Eigen::MatrixXd d= Eigen::MatrixXd::Zero(3,3); //the first row has all zeros  
  
  
 // List Wi(num_obs); 
  //for(int i_Wi=0;i_Wi<num_dim;i_Wi++){ 
    //Eigen::MatrixXd d= Eigen::MatrixXd::Zero(num_obs,9);   
    double  lambda_delta_x;
    double exp_neg_2_lambda_delta_x;
    int  j_Wi_1;
    for(int j_Wi=0;j_Wi<(num_obs-1);j_Wi++){ 
      j_Wi_1= j_Wi+1; 
      lambda_delta_x=lambda*delta_x[j_Wi];  //close and jump then it is... 
      exp_neg_2_lambda_delta_x=exp(-2*lambda_delta_x); 
      
      d(0,0)=(exp_neg_2_lambda_delta_x*(3+6*lambda_delta_x+6*pow(lambda_delta_x,2.0)+4*pow(lambda_delta_x,3.0)+2*pow(lambda_delta_x,4.0))-3 )/(-4*pow(lambda,5.0)); 
      d(1, 0)=  d(0, 1)=exp_neg_2_lambda_delta_x*pow(delta_x[j_Wi],4.0)/2.0; 
      d(2, 0)=  d(0, 2)=(exp_neg_2_lambda_delta_x*(1+2*lambda_delta_x+2*pow(lambda_delta_x,2.0)+4*pow(lambda_delta_x,3.0)-2*pow(lambda_delta_x,4.0))-1 )/(4*pow(lambda,3.0)); 
      d(1, 1)= (exp_neg_2_lambda_delta_x*(1+2*lambda_delta_x+2*pow(lambda_delta_x,2.0)-4*pow(lambda_delta_x,3.0)+2*pow(lambda_delta_x,4.0))-1 )/(-4*pow(lambda,3.0)); 
      d(2, 1)=  d(1, 2)=exp_neg_2_lambda_delta_x*pow(delta_x[j_Wi],2.0)*(4-4*lambda_delta_x+pow(lambda_delta_x,2.0) )/2.0; 
      d(2, 2)=(exp_neg_2_lambda_delta_x*(-3+10*lambda_delta_x-22*pow(lambda_delta_x,2.0)+12*pow(lambda_delta_x,3.0)-2*pow(lambda_delta_x,4.0))+3 )/(4*lambda)     ;  
      d=d*(4*sigma2*pow(lambda,5.0)/3.0); 
      Wi[j_Wi_1]=d; 
      
    //} 
  } 
  return Wi; 
} 

////Construct_W_matern_5_2  
// [[Rcpp::export]] 
List Construct_W_exp(double sigma2, Eigen::VectorXd delta_x, double lambda, MatrixXd W0){  //be careful about delta_x,lambda if one only wants to sample one 
  int num_obs=delta_x.size()+1; 
  //int num_dim=sigma2.size();  
  List Wi(num_obs);  
  Wi[0]=W0; 
  Eigen::MatrixXd d= Eigen::MatrixXd::Zero(1,1); 
  
  for(int j_Wi=0;j_Wi<(num_obs-1);j_Wi++){ 
    d(0,0)=1-exp(-2*delta_x[j_Wi]*lambda);
    Wi[j_Wi+1]=d;
  }
  
  return Wi;
}

////Get_Q_K  
// [[Rcpp::export]] 
List Get_Q_K(const List GG,const List  W,const Eigen::MatrixXd C0,const double VV){ 

   int n=GG.size();
   int k=C0.rows();
   
   Eigen::VectorXd Q=Eigen::VectorXd::Zero(n);
   Eigen::MatrixXd K=Eigen::MatrixXd::Zero(n,k);
   Eigen::MatrixXd C=C0;
    
   Eigen::MatrixXd GG_matrix;
   Eigen::MatrixXd W_matrix;
   
   Eigen::MatrixXd RR;
   
      
   // num_dim list, each is 3(num_obs)\times 3 list 
   for(int t=0;t<n;t++){ 
     GG_matrix=GG[t];
     W_matrix=W[t];
     RR=GG_matrix*C*GG_matrix.transpose()+W_matrix;
     //Q[t]=RR(0,0);
     Q[t]=RR(0,0)+VV;
     K.row(t)=RR.col(0).transpose()/Q[t];
     C=RR-RR.col(0)*RR.row(0)/Q[t];
   }

   List return_list;
   return_list.push_back(Q);
   return_list.push_back(K);
   
   return return_list;
}



     
     
////Get_Y_minus_a_1_scaled_matrix_2d  
// [[Rcpp::export]] 
Eigen::MatrixXd Get_Y_minus_a_1_scaled_matrix_2d(const Eigen::MatrixXd output_KF,const List GG,const Eigen::VectorXd Q,const Eigen::MatrixXd K){
  
     int n1=output_KF.rows();
     int n2=output_KF.cols();
     int k=K.cols();
     
     Eigen::MatrixXd m=Eigen::MatrixXd::Zero(k,n2);

     Eigen::MatrixXd Y_minus_a_1_scaled_matrix=Eigen::MatrixXd::Zero(n1,n2); 
     Eigen::VectorXd sqrt_Q=Q.array().sqrt();

     Eigen::MatrixXd GG_matrix;
     
     Eigen::MatrixXd a;
     for(int t=0;t<n1;t++){
       GG_matrix=GG[t];
       a=GG_matrix*m;
       Y_minus_a_1_scaled_matrix.row(t)=(output_KF.row(t)-a.row(0))/sqrt_Q[t];
       m=a+K.row(t).transpose()*(output_KF.row(t)-a.row(0));
     }
     
     
     return Y_minus_a_1_scaled_matrix;
}




  
  
// [[Rcpp::export]] 
List Get_G_log_det_cov(const Eigen::VectorXd param,const Eigen::MatrixXd output,const Eigen::VectorXd delta_x,int d,
                       const String kernel_type){
        Eigen::VectorXd gamma=(1.0/param.head(d).array().exp()).matrix();
        Eigen::VectorXd tau=param.tail(d).array().exp().matrix();
        

        Eigen::MatrixXd    output_t=output.transpose();
        
        Eigen::VectorXd Q;


        Eigen::MatrixXd W0;
        List GG;
        List W;
        List Q_K;
        
        Eigen::VectorXd VV=1.0/tau.array();
        Eigen::VectorXd log_det_cov=Eigen::VectorXd::Zero(d);
          
        //int n=output.cols();
          
        Eigen::MatrixXd    z;
          
        List return_list;
        
        Eigen::VectorXd lambda;
        if(kernel_type=="matern_5_2"){
        // lambda=sqrt(5.0)/gamma;
          lambda=(sqrt(5.0)/gamma.array()).matrix();
          
          for(int i=0;i<d;i++){  
             W0=Construct_W0_matern_5_2(1.0,lambda[i]);  
            
             GG=Construct_G_matern_5_2(delta_x,lambda[i]);  
            
             W=Construct_W_matern_5_2(1.0,delta_x,lambda[i],W0);
            
             Q_K=Get_Q_K(GG,W,W0,VV[i]);
             
             Q=Q_K[0];
             
             log_det_cov[i]=(tau[i]*Q.array()).log().sum();
            
             z=Get_Y_minus_a_1_scaled_matrix_2d(output_t,GG,Q_K[0],Q_K[1]);
             
              return_list.push_back( (1.0/tau[i]*(z.transpose()*z).array()).matrix());
            // return_list.push_back( ((z.transpose()*z).array()).matrix());
             
          }
        }else if(kernel_type=="exp"){
          lambda=(1.0/gamma.array()).matrix();
          
          for(int i=0;i<d;i++){  
            W0=Construct_W0_exp(1.0,lambda[i]);  
            GG=Construct_G_exp(delta_x,lambda[i]);  
            W=Construct_W_exp(1.0,delta_x,lambda[i],W0);
            
            

            Q_K=Get_Q_K(GG,W,W0,VV[i]);
            
            Q=Q_K[0];
            
            log_det_cov[i]=(tau[i]*Q.array()).log().sum();
            
            z=Get_Y_minus_a_1_scaled_matrix_2d(output_t,GG,Q_K[0],Q_K[1]);
            
            return_list.push_back( (1.0/tau[i]*(z.transpose()*z).array()).matrix());
            // return_list.push_back( ((z.transpose()*z).array()).matrix());
            
          }
          
          }
        return_list.push_back(log_det_cov);
        
        return return_list;
}


/*
//[[Rcpp::export]] 
List  Posterior_AZ_funct(const Eigen::VectorXd param,const Eigen::VectorXd delta_x, const Eigen::MatrixXd A_hat,
                           const Eigen::MatrixXd output,int d){
  
  //double beta_hat=param[0];
  //double sigma_2_hat=param[1];
  //double sigma_2_0_hat=param[2];
  
  
  Eigen::VectorXd gamma=(1.0/param.head(d).array()).matrix();
  Eigen::VectorXd sigma_2=param.segment(d,d);
  double sigma_2_0=param[2*d];
  
  Eigen::VectorXd tau=(sigma_2.array()/sigma_2_0).matrix();
  Eigen::VectorXd    lambda=(sqrt(5.0)/gamma.array()).matrix();
  
  Eigen::VectorXd VV=1.0/tau.array();
  
  Eigen::MatrixXd    output_t=output.transpose();
  
  Eigen::MatrixXd    z;
  
  
  List return_list;
  
  Eigen::VectorXd Q;
  
  
  Eigen::MatrixXd W0;
  List GG;
  List W;
  List Q_K;
  
  Eigen::VectorXd output_t_a_i;
  
  for(int i=0;i<d;i++){  
    W0=Construct_W_matern_5_20_matern_5_2(1.0,lambda[i]);  
    
    GG=Construct_G_matern_5_2(delta_x,lambda[i]);  
    
    W=Construct_W_matern_5_2(1.0,delta_x,lambda[i],W0);
    
    Q_K=Get_Q_K(GG,W,W0,VV[i]);
    
    
    //Q=Q_K[0];
    
    //log_det_cov[i]=(tau[i]*Q.array()).log().sum();
    
    output_t_a_i=output_t*A_hat.col(i);
    
    z=Get_Y_minus_a_1_scaled_matrix_2d(output_t_a_i,GG,Q_K[0],Q_K[1]);
    
    return_list.push_back( (1.0/tau[i]*(z.transpose()*z).array()).matrix());
    // return_list.push_back( ((z.transpose()*z).array()).matrix());
    
  }
  

  
}
*/


// [[Rcpp::export]] 
double F_Funct(const Eigen::MatrixXd A_cur,const List G){
  double return_val=0;
  int d=A_cur.cols();
  Eigen::MatrixXd G_matrix;
  for(int i=0;i<d;i++){
    G_matrix=G[i];
    return_val=return_val+A_cur.col(i).transpose()*G_matrix*A_cur.col(i);
  }
  return -return_val;
}
// [[Rcpp::export]] 
Eigen::MatrixXd F_Funct_Dev(const Eigen::MatrixXd A_cur,const List G){
  int k=A_cur.rows();
  int d=A_cur.cols();
  Eigen::MatrixXd return_matrix=Eigen::MatrixXd::Zero(k,d); 
  Eigen::MatrixXd G_matrix;
  
  for(int i=0;i<d;i++){
    G_matrix=G[i];
    return_matrix.col(i)=G_matrix*A_cur.col(i);
  }
  return -2*return_matrix;
}






//[[Rcpp::export]] 
List Get_B_U_V(const Eigen::MatrixXd A_cur,const List G){
  Eigen::MatrixXd B= F_Funct_Dev(A_cur,G);
  int k=A_cur.rows();
  int d=A_cur.cols();
  Eigen::MatrixXd U=Eigen::MatrixXd::Zero(k,2*d); 
  Eigen::MatrixXd V=Eigen::MatrixXd::Zero(k,2*d); 
  
  U.leftCols(d)=B;
  U.rightCols(d)=A_cur;
  
  V.leftCols(d)=A_cur;
  V.rightCols(d)=-B;
  
  List return_list;
  
  return_list.push_back(B);
  return_list.push_back(U);
  return_list.push_back(V);
  
  return return_list;
  
}


//[[Rcpp::export]] 
List Y_Funct(const Eigen::MatrixXd A_cur, const List B_U_V,  double tau){
  int d=A_cur.cols();
  int k=A_cur.rows();
  Eigen::MatrixXd B=B_U_V[0];
  Eigen::MatrixXd U=B_U_V[1];
  Eigen::MatrixXd V=B_U_V[2];
  
  
  //I may need to consider what if this is singular
  
  Eigen::MatrixXd middle_middle_term=Eigen::MatrixXd::Identity(2*d,2*d)+tau/2.0*V.transpose()*U;
    
  JacobiSVD<MatrixXd> svd(middle_middle_term);
  double cond = svd.singularValues()(0)/svd.singularValues()(svd.singularValues().size()-1);
  
  //cout << tau;    
  
  while(cond>pow(10.0,16.0)){
   // tau=tau/(2*log(cond/pow(10.0,15.0)+1));
    tau=tau/2;
    
    middle_middle_term=Eigen::MatrixXd::Identity(2*d,2*d)+tau/2.0*V.transpose()*U;
    
    JacobiSVD<MatrixXd>  svd(middle_middle_term);
    
    cond = svd.singularValues()(0)/svd.singularValues()(svd.singularValues().size()-1);
    
  }
  //cout << tau;    
  

  Eigen::MatrixXd middle_term=U*((Eigen::MatrixXd::Identity(2*d,2*d)+tau/2.0*V.transpose()*U).lu().solve(V.transpose()));
 //Eigen::MatrixXd middle_term=U*((middle_middle_term).lu().solve(V.transpose()));
  
  Eigen::MatrixXd Y_tau= (Eigen::MatrixXd::Identity(k,k)-tau*middle_term)*A_cur;
  
  // Eigen::MatrixXd Y_tau= A_cur -tau*middle_term*A_cur;
    
  Eigen::MatrixXd Y_dev_tau=-middle_term*(A_cur+Y_tau)/2;
  
  List return_list;
  return_list.push_back(Y_tau);
  return_list.push_back(Y_dev_tau);
  return_list.push_back(tau);
  
  return return_list;
    
}

    
//[[Rcpp::export]] 
Eigen::MatrixXd Optimization_Stiefel_Manifold(const Eigen::MatrixXd A_ini, const List G, int max_iter){
  int k=A_ini.cols();
  int d=A_ini.rows();
  
  double rho_1=pow(10.0,-4.0);
  double delta=0.2;
  double eta=0.85;
  double epsilon_1=pow(10.0,-5.0);
  double epsilon_2=pow(10.0,-10.0);
  
  double C_cur=F_Funct(A_ini,G);
  Eigen::MatrixXd  A_cur=A_ini;
  double  Q_cur=1.0;
  //double tau_cur=0.001;
  
  //List B_U_V;
  
  List B_U_V=Get_B_U_V(A_cur, G);
  Eigen::MatrixXd B=B_U_V[0];
  Eigen::MatrixXd U=B_U_V[1];
  Eigen::MatrixXd V=B_U_V[2];
  
  //double tau_cur=0.01;
  
  
  double tau_cur=1.0/((V.transpose()*U).diagonal().array().abs().sum());
    
 //   1/sum(abs(diag(t(V)%*%U)));
    
  
  Eigen::MatrixXd gradient_F_Y_tau=B-A_cur*B.transpose()*A_cur;
  
  double F_Y_tau=F_Funct(A_cur,G);
  double F_Y_0;
  
  bool find_tau;
  Eigen::MatrixXd gradient_F_A_cur;
  double norm_gradient_cur;
  
  double F_cur_val=F_Y_tau;
  double F_last_val=F_cur_val-1;
  
  List Y_tau_dev_tau_all;
  
  Eigen::MatrixXd Y_tau;
  Eigen::MatrixXd Y_tau_dev;

  Eigen::MatrixXd Y_0_dev;
  
  Eigen::MatrixXd B_A;
  
  Eigen::MatrixXd diff_graident;
  double F_dev_Y_0;
  
  Eigen::MatrixXd S;
  
  double SS;
  double SY_abs;
  double YY;
  
  double tau_next;  
  
  double Q_next;
  
  int count;
  for(int i_iter=0; i_iter<max_iter;i_iter++){
    find_tau=true;
    //B_U_V=B_U_V_next;
    
    gradient_F_A_cur=gradient_F_Y_tau;
    norm_gradient_cur=pow(gradient_F_A_cur.array().pow(2.0).sum(),0.5);
    
    F_cur_val=F_Y_tau;
      
    if(i_iter>1){
      if((norm_gradient_cur/(k*d))<epsilon_1 ||  (abs(F_cur_val- F_last_val))<epsilon_2 ){
        break;
      }
    }
    
    F_last_val=F_cur_val;
    
    count=0;
    while(find_tau || count>50 ){
      count++;
      Y_tau_dev_tau_all=Y_Funct(A_cur,B_U_V,tau_cur);
      
      Y_tau=Y_tau_dev_tau_all[0];
      Y_tau_dev=Y_tau_dev_tau_all[1];
      tau_cur=Y_tau_dev_tau_all[2];
      
      
      Y_0_dev=-U*V.transpose()*A_cur;
      
      F_Y_tau=F_Funct(Y_tau,G);
      
      F_Y_0=F_Funct(A_cur,G);
      
      B_A=F_Funct_Dev(A_cur,G);
        
      F_dev_Y_0=(B_A.array()*Y_0_dev.array()).sum();
      
        
     if( (F_Y_tau- (C_cur+rho_1*tau_cur*F_dev_Y_0)) <=0){
       find_tau=false;
     }else{
       tau_cur=delta*tau_cur;
     }
    }
        
      
     //t(Y_tau)%*%Y_tau
     B_U_V=Get_B_U_V(Y_tau, G);
        
     B=B_U_V[0];
     U=B_U_V[1];
     V=B_U_V[2];
        
     //compute that trace thing for tau
     S=Y_tau-A_cur;
          
     gradient_F_Y_tau=B-Y_tau*B.transpose()*Y_tau;
     
     diff_graident=gradient_F_Y_tau-gradient_F_A_cur;
       
       SS=(S.array()*S.array()).sum();
       SY_abs=abs((S.array()*diff_graident.array()).sum());
       YY=(diff_graident.array()*diff_graident.array()).sum();

       if(i_iter%2==0){
         tau_next=SS/SY_abs;
       }else{
         tau_next=SY_abs/YY;
       }
       
       A_cur=Y_tau;
       Q_next=eta*Q_cur+1;
       C_cur=(eta*Q_cur*C_cur+F_Y_tau)/Q_next;
       Q_cur=Q_next;
         
       //this is controversial   
       if(!isnan(tau_next)){
          tau_cur=max(min(tau_next,pow(10.0,20)),pow(10.0,-20));
       }
         
       //tau_cur=max(min(tau_next,pow(10.0,20)),pow(10.0,-20));
      
  }
     
     return A_cur; 
  
}
  
  
  
  
  


////update_FRFt 
// [[Rcpp::export]] 
Eigen::MatrixXd Chol_rcppeigen(const Eigen::MatrixXd A){
  
  LLT<MatrixXd> lltOfR(A);    // compute the cholesky decomposition of R called lltofR
  MatrixXd L = lltOfR.matrixL();   //retrieve factor L  in the decomposition
  return L;
}

  
  
  
  
//March 20, 2019
//The following code are for the function with a mean structure

//This function get and all the matrices for Kalman filtering  
//no value is required


// [[Rcpp::export]] 
List Prep_Kalman_filtering_GPPCA(const Eigen::VectorXd param,const Eigen::VectorXd delta_x,int d, const String kernel_type){
  
  Eigen::VectorXd gamma=(1.0/param.head(d).array().exp()).matrix();
  Eigen::VectorXd tau=param.tail(d).array().exp().matrix();
  
  Eigen::MatrixXd W0;
  List GG;
  List W;
  List Q_K;
  
  List W0_list;
  List GG_list;
  List W_list;
  List Q_K_list;
  
  Eigen::VectorXd VV=1.0/tau.array();

  List return_list;
  
  Eigen::VectorXd lambda;
  if(kernel_type=="matern_5_2"){
    // lambda=sqrt(5.0)/gamma;
    lambda=(sqrt(5.0)/gamma.array()).matrix();
    
    for(int i=0;i<d;i++){  
      W0=Construct_W0_matern_5_2(1.0,lambda[i]);  
      
      GG=Construct_G_matern_5_2(delta_x,lambda[i]);  
      
      W=Construct_W_matern_5_2(1.0,delta_x,lambda[i],W0);
      
      Q_K=Get_Q_K(GG,W,W0,VV[i]);
      
      W0_list.push_back(W0);
      GG_list.push_back(GG);
      W_list.push_back(W);
      Q_K_list.push_back(Q_K);
    }
  }else if(kernel_type=="exp"){
    lambda=(1.0/gamma.array()).matrix();
    
    for(int i=0;i<d;i++){  
      W0=Construct_W0_exp(1.0,lambda[i]);  
      GG=Construct_G_exp(delta_x,lambda[i]);  
      W=Construct_W_exp(1.0,delta_x,lambda[i],W0);
      Q_K=Get_Q_K(GG,W,W0,VV[i]);
      
      W0_list.push_back(W0);
      GG_list.push_back(GG);
      W_list.push_back(W);
      Q_K_list.push_back(Q_K);
      
    }
  }
  return_list.push_back(W0_list);
  return_list.push_back(GG_list);
  return_list.push_back(W_list);
  return_list.push_back(Q_K_list);
  
  return return_list;
}

// [[Rcpp::export]] 
double F_Funct_with_mean(const Eigen::MatrixXd A_cur, const Eigen::MatrixXd output, const Eigen::MatrixXd H,const Eigen::MatrixXd M, 
                         const Eigen::VectorXd tau,const List W0_list, const List GG_list,
                         const List W_list,const List Q_K_list){
  //int k=A_cur.rows();
  int d=A_cur.cols();
  
  Eigen::MatrixXd H_half;
  Eigen::MatrixXd y_a_half;
  Eigen::MatrixXd H_T_tau_R_tilde_inv_H;
  Eigen::VectorXd H_T_tau_R_tilde_inv_y_a;
  
  Eigen::VectorXd theta_a;
  Eigen::VectorXd y_a_tilde;
  Eigen::VectorXd z;
  Eigen::VectorXd a_t_G_a=VectorXd::Zero(d);
  
  Eigen::MatrixXd A_t_Y=A_cur.transpose()*output;
    
  List Q_K;  
  for(int i_d=0;i_d<d;i_d++){
     Q_K=Q_K_list[i_d];
      
     H_half=Get_Y_minus_a_1_scaled_matrix_2d(H,GG_list[i_d],Q_K[0],Q_K[1]);
     y_a_half=Get_Y_minus_a_1_scaled_matrix_2d(A_t_Y.row(i_d).transpose(),GG_list[i_d],Q_K[0],Q_K[1]);
     H_T_tau_R_tilde_inv_H=1.0/tau[i_d]*H_half.transpose()*H_half;
     H_T_tau_R_tilde_inv_y_a=1.0/tau[i_d]*H_half.transpose()*y_a_half;
     
     theta_a=H_T_tau_R_tilde_inv_H.lu().solve(H_T_tau_R_tilde_inv_y_a);
     
     
     
     //H_T_tau_R_tilde_inv_H_inv=solve(H_T_tau_R_tilde_inv_H)
     //theta_a=H_T_tau_R_tilde_inv_H_inv%*%H_T_tau_R_tilde_inv_y_a
      
     y_a_tilde=A_t_Y.row(i_d).transpose()-H*theta_a;
     //  y_a-H*theta_a;
     z=Get_Y_minus_a_1_scaled_matrix_2d(y_a_tilde,GG_list[i_d],Q_K[0],Q_K[1]);
     a_t_G_a[i_d]=(A_t_Y.row(i_d)*M*A_t_Y.row(i_d).transpose()- 1.0/tau[i_d]*z.transpose()*z)(0);
      
  }
  
  return -a_t_G_a.sum();
}


// [[Rcpp::export]] 
MatrixXd F_Funct_with_mean_dev(const Eigen::MatrixXd A_cur, const Eigen::MatrixXd output, const Eigen::MatrixXd H,const Eigen::MatrixXd M, 
                        const Eigen::VectorXd tau,const List W0_list, const List GG_list,
                         const List W_list,const List Q_K_list){
  int k=A_cur.rows();
  int d=A_cur.cols();
  
  
  Eigen::MatrixXd H_half;
  Eigen::MatrixXd y_a_half;
  Eigen::MatrixXd H_T_tau_R_tilde_inv_H;
  Eigen::VectorXd H_T_tau_R_tilde_inv_y_a;
  
  Eigen::VectorXd theta_a;
  Eigen::VectorXd y_a_tilde;
  Eigen::VectorXd z1;
  Eigen::MatrixXd z2;
  
  Eigen::MatrixXd return_matrix=Eigen::MatrixXd::Zero(k,d); 
  
  //Eigen::VectorXd dev_terms=VectorXd::Zero(d);
  
  Eigen::MatrixXd A_t_Y=A_cur.transpose()*output;
  
  List Q_K;  
  for(int i_d=0;i_d<d;i_d++){
    Q_K=Q_K_list[i_d];
    
    H_half=Get_Y_minus_a_1_scaled_matrix_2d(H,GG_list[i_d],Q_K[0],Q_K[1]);
    y_a_half=Get_Y_minus_a_1_scaled_matrix_2d(A_t_Y.row(i_d).transpose(),GG_list[i_d],Q_K[0],Q_K[1]);
    H_T_tau_R_tilde_inv_H=1.0/tau[i_d]*H_half.transpose()*H_half;
    H_T_tau_R_tilde_inv_y_a=1.0/tau[i_d]*H_half.transpose()*y_a_half;
    
    theta_a=H_T_tau_R_tilde_inv_H.lu().solve(H_T_tau_R_tilde_inv_y_a);
    
    
    //H_T_tau_R_tilde_inv_H_inv=solve(H_T_tau_R_tilde_inv_H)
    //theta_a=H_T_tau_R_tilde_inv_H_inv%*%H_T_tau_R_tilde_inv_y_a
    
    y_a_tilde=A_t_Y.row(i_d).transpose()-H*theta_a;
    //  y_a-H*theta_a;
    z1=Get_Y_minus_a_1_scaled_matrix_2d(y_a_tilde,GG_list[i_d],Q_K[0],Q_K[1]);
    z2=Get_Y_minus_a_1_scaled_matrix_2d(output.transpose(),GG_list[i_d],Q_K[0],Q_K[1]);
    

    return_matrix.col(i_d)=output*M*A_t_Y.row(i_d).transpose()-1.0/tau[i_d]*z2.transpose()*z1;
  }
  return -2*return_matrix;
  
}
    //  (A_t_Y.row(i_d)*M*A_t_Y.row(i_d).transpose()- 1.0/tau[i_d]*z.transpose()*z)(0);
    
    
    
    
//March 25, code to explore the large k small n strategy
    
// [[Rcpp::export]] 
double F_Funct_Large_k(const Eigen::MatrixXd A_cur,const List UD){
      double return_val=0;
      int d=A_cur.cols();
      Eigen::VectorXd z;
      Eigen::MatrixXd UD_matrix;
      for(int i=0;i<d;i++){
         UD_matrix=UD[i];
        z=(A_cur.col(i).transpose()*UD_matrix).transpose();
        //t(X[,i])%*%UD[[i]]
        return_val=return_val+(z.transpose()*z);
      // A_cur.col(i).transpose()*G_matrix*A_cur.col(i);
      }
      return -return_val;
}

// [[Rcpp::export]] 
Eigen::MatrixXd F_Funct_Dev_Large_k(const Eigen::MatrixXd A_cur,const List UD){
  int k=A_cur.rows();
  int d=A_cur.cols();
  Eigen::MatrixXd return_matrix=Eigen::MatrixXd::Zero(k,d); 
  Eigen::MatrixXd UD_matrix;
  Eigen::VectorXd z;
  
  for(int i=0;i<d;i++){
    UD_matrix=UD[i];
    z=(A_cur.col(i).transpose()*UD_matrix).transpose();
    //return_matrix[,i]=UD[[i]]%*%t(z_t)
      
    return_matrix.col(i)=UD_matrix*z;
  }
  return -2*return_matrix;
}


//[[Rcpp::export]] 
List Get_B_U_V_Large_k(const Eigen::MatrixXd A_cur,const List UD){
  Eigen::MatrixXd B= F_Funct_Dev_Large_k(A_cur,UD);
  int k=A_cur.rows();
  int d=A_cur.cols();
  Eigen::MatrixXd U=Eigen::MatrixXd::Zero(k,2*d); 
  Eigen::MatrixXd V=Eigen::MatrixXd::Zero(k,2*d); 
  
  U.leftCols(d)=B;
  U.rightCols(d)=A_cur;
  
  V.leftCols(d)=A_cur;
  V.rightCols(d)=-B;
  
  List return_list;
  
  return_list.push_back(B);
  return_list.push_back(U);
  return_list.push_back(V);
  
  return return_list;
  
}



//[[Rcpp::export]] 
List Y_Funct_Large_k(const Eigen::MatrixXd A_cur, const List B_U_V,  double tau){
  int d=A_cur.cols();
  //int k=A_cur.rows();
  Eigen::MatrixXd B=B_U_V[0];
  Eigen::MatrixXd U=B_U_V[1];
  Eigen::MatrixXd V=B_U_V[2];
  
  
  //I may need to consider what if this is singular
  
  Eigen::MatrixXd middle_middle_term=Eigen::MatrixXd::Identity(2*d,2*d)+tau/2.0*V.transpose()*U;
  
  JacobiSVD<MatrixXd> svd(middle_middle_term);
  double cond = svd.singularValues()(0)/svd.singularValues()(svd.singularValues().size()-1);
  
  //cout << tau;    
  
  while(cond>pow(10.0,16.0)){
    // tau=tau/(2*log(cond/pow(10.0,15.0)+1));
    tau=tau/2;
    
    middle_middle_term=Eigen::MatrixXd::Identity(2*d,2*d)+tau/2.0*V.transpose()*U;
    
    JacobiSVD<MatrixXd>  svd(middle_middle_term);
    
    cond = svd.singularValues()(0)/svd.singularValues()(svd.singularValues().size()-1);
    
  }
  //cout << tau;    
  
  //this term is k by k...
  //Eigen::MatrixXd middle_term=U*((Eigen::MatrixXd::Identity(2*d,2*d)+tau/2.0*V.transpose()*U).lu().solve(V.transpose()));

  //Eigen::MatrixXd Y_tau= (Eigen::MatrixXd::Identity(k,k)-tau*middle_term)*A_cur;
  
  Eigen::MatrixXd middle_term_no_U=((Eigen::MatrixXd::Identity(2*d,2*d)+tau/2.0*V.transpose()*U).lu().solve(V.transpose()));
  
  
   Eigen::MatrixXd Y_tau= A_cur -tau*U*(middle_term_no_U*A_cur);
  
  Eigen::MatrixXd Y_dev_tau=- U*(middle_term_no_U*(A_cur+Y_tau))/2;
  
  List return_list;
  return_list.push_back(Y_tau);
  return_list.push_back(Y_dev_tau);
  return_list.push_back(tau);
  
  return return_list;
  
}

//[[Rcpp::export]] 
Eigen::MatrixXd Optimization_Stiefel_Manifold_Large_k(const Eigen::MatrixXd A_ini, const List UD, int max_iter){
  int k=A_ini.cols();
  int d=A_ini.rows();
  
  double rho_1=pow(10.0,-4.0);
  double delta=0.2;
  double eta=0.85;
  double epsilon_1=pow(10.0,-5.0);
  double epsilon_2=pow(10.0,-10.0);
  
  double C_cur=F_Funct_Large_k(A_ini,UD);
  Eigen::MatrixXd  A_cur=A_ini;
  double  Q_cur=1.0;
  //double tau_cur=0.001;
  
  //List B_U_V;
  
  List B_U_V=Get_B_U_V_Large_k(A_cur, UD);
  Eigen::MatrixXd B=B_U_V[0];
  Eigen::MatrixXd U=B_U_V[1];
  Eigen::MatrixXd V=B_U_V[2];
  
  //double tau_cur=0.01;
  
  
  double tau_cur=1.0/((V.transpose()*U).diagonal().array().abs().sum());
  
  //   1/sum(abs(diag(t(V)%*%U)));
  
  
  Eigen::MatrixXd gradient_F_Y_tau=B-A_cur*(B.transpose()*A_cur);
  
  double F_Y_tau=F_Funct_Large_k(A_cur,UD);
  double F_Y_0;
  
  bool find_tau;
  Eigen::MatrixXd gradient_F_A_cur;
  double norm_gradient_cur;
  
  double F_cur_val=F_Y_tau;
  double F_last_val=F_cur_val-1;
  
  List Y_tau_dev_tau_all;
  
  Eigen::MatrixXd Y_tau;
  Eigen::MatrixXd Y_tau_dev;
  
  Eigen::MatrixXd Y_0_dev;
  
  Eigen::MatrixXd B_A;
  
  Eigen::MatrixXd diff_graident;
  double F_dev_Y_0;
  
  Eigen::MatrixXd S;
  
  double SS;
  double SY_abs;
  double YY;
  
  double tau_next;  
  
  double Q_next;
  
  int count;
  for(int i_iter=0; i_iter<max_iter;i_iter++){
    find_tau=true;
    //B_U_V=B_U_V_next;
    
    gradient_F_A_cur=gradient_F_Y_tau;
    norm_gradient_cur=pow(gradient_F_A_cur.array().pow(2.0).sum(),0.5);
    
    F_cur_val=F_Y_tau;
    
    if(i_iter>1){
      if((norm_gradient_cur/(k*d))<epsilon_1 ||  (abs(F_cur_val- F_last_val))<epsilon_2 ){
        break;
      }
    }
    
    F_last_val=F_cur_val;
    
    count=0;
    //cout << count;   
    
    while(find_tau || count>50 ){
      count++;
      Y_tau_dev_tau_all=Y_Funct_Large_k(A_cur,B_U_V,tau_cur);
      
      Y_tau=Y_tau_dev_tau_all[0];
      Y_tau_dev=Y_tau_dev_tau_all[1];
      tau_cur=Y_tau_dev_tau_all[2];
      
      
      Y_0_dev=-U*(V.transpose()*A_cur);
      
      F_Y_tau=F_Funct_Large_k(Y_tau,UD);
      
      F_Y_0=F_Funct_Large_k(A_cur,UD);
      
      B_A=F_Funct_Dev_Large_k(A_cur,UD);
      
      F_dev_Y_0=(B_A.array()*Y_0_dev.array()).sum();
      
      
      if( (F_Y_tau- (C_cur+rho_1*tau_cur*F_dev_Y_0)) <=0){
        find_tau=false;
      }else{
        tau_cur=delta*tau_cur;
      }
    }
    
    //cout << count<<endl;    
    
    //t(Y_tau)%*%Y_tau
    B_U_V=Get_B_U_V_Large_k(Y_tau, UD);
    //cout << 2;    
    
    B=B_U_V[0];
    U=B_U_V[1];
    V=B_U_V[2];
    
    //compute that trace thing for tau
    S=Y_tau-A_cur;
    
    gradient_F_Y_tau=B-Y_tau*(B.transpose()*Y_tau);
    
    diff_graident=gradient_F_Y_tau-gradient_F_A_cur;
    
    SS=(S.array()*S.array()).sum();
    SY_abs=abs((S.array()*diff_graident.array()).sum());
    YY=(diff_graident.array()*diff_graident.array()).sum();
    
    if(i_iter%2==0){
      tau_next=SS/SY_abs;
    }else{
      tau_next=SY_abs/YY;
    }
    
    A_cur=Y_tau;
    Q_next=eta*Q_cur+1;
    C_cur=(eta*Q_cur*C_cur+F_Y_tau)/Q_next;
    Q_cur=Q_next;
    
    //this is controversial   
    if(!isnan(tau_next)){
      tau_cur=max(min(tau_next,pow(10.0,20)),pow(10.0,-20));
    }
    
    //tau_cur=max(min(tau_next,pow(10.0,20)),pow(10.0,-20));
    
  }
  
  return A_cur; 
  
}









///old codes///////////////////////////////////////////////////////

////construct_Q0 
// [[Rcpp::export]] 
List Construct_Q0(const Eigen::VectorXd sigma2, const Eigen::VectorXd lambda){ 
  int num_dim=sigma2.size(); 
  List  Q0(num_dim); 
  for(int i_Q0=0;i_Q0<num_dim;i_Q0++){  // Is there a way to parallel 
    Eigen::MatrixXd d= Eigen::MatrixXd::Zero(3,3); 
    d(0,0)=sigma2[i_Q0]; 
    d(0,2)=d(2,0)=-sigma2[i_Q0]*pow(lambda[i_Q0],2.0)/3.0; 
    d(1,1)=sigma2[i_Q0]*pow(lambda[i_Q0],2.0)/3.0; 
    d(2,2)=sigma2[i_Q0]*pow(lambda[i_Q0],4.0); 
    Q0[i_Q0]=d; 
    
  } 
  return Q0; 
} 


////construct_Qi 
// [[Rcpp::export]] 
List Construct_Qi(Eigen::VectorXd sigma2,Eigen::VectorXd delta_x, Eigen::VectorXd lambda, List Q0){  //be careful about delta_x,lambda if one only wants to sample one 
  int num_obs=delta_x.size()+1; 
  int num_dim=sigma2.size();  
  
  List Qi(num_dim); 
  for(int i_Qi=0;i_Qi<num_dim;i_Qi++){ 
    Eigen::MatrixXd d= Eigen::MatrixXd::Zero(num_obs,9);   
    
    for(int j_Qi=0;j_Qi<(num_obs-1);j_Qi++){ 
      int  j_Qi_1= j_Qi+1; 
      double  lambda_delta_x=lambda[i_Qi]*delta_x[j_Qi];  //close and jump then it is... 
      double exp_neg_2_lambda_delta_x=exp(-2*lambda_delta_x); 
      
      d(j_Qi_1, 0)=(exp_neg_2_lambda_delta_x*(3+6*lambda_delta_x+6*pow(lambda_delta_x,2.0)+4*pow(lambda_delta_x,3.0)+2*pow(lambda_delta_x,4.0))-3 )/(-4*pow(lambda[i_Qi],5.0)); 
      d(j_Qi_1, 1)=  d(j_Qi_1, 3)=exp_neg_2_lambda_delta_x*pow(delta_x[j_Qi],4.0)/2.0; 
      d(j_Qi_1, 2)=  d(j_Qi_1, 6)=(exp_neg_2_lambda_delta_x*(1+2*lambda_delta_x+2*pow(lambda_delta_x,2.0)+4*pow(lambda_delta_x,3.0)-2*pow(lambda_delta_x,4.0))-1 )/(4*pow(lambda[i_Qi],3.0)); 
      d(j_Qi_1, 4)= (exp_neg_2_lambda_delta_x*(1+2*lambda_delta_x+2*pow(lambda_delta_x,2.0)-4*pow(lambda_delta_x,3.0)+2*pow(lambda_delta_x,4.0))-1 )/(-4*pow(lambda[i_Qi],3.0)); 
      d(j_Qi_1, 5)=  d(j_Qi_1, 7)=exp_neg_2_lambda_delta_x*pow(delta_x[j_Qi],2.0)*(4-4*lambda_delta_x+pow(lambda_delta_x,2.0) )/2.0; 
      d(j_Qi_1, 8)=(exp_neg_2_lambda_delta_x*(-3+10*lambda_delta_x-22*pow(lambda_delta_x,2.0)+12*pow(lambda_delta_x,3.0)-2*pow(lambda_delta_x,4.0))+3 )/(4*lambda[i_Qi])     ;  
    } 
    d=d*(4*sigma2[i_Qi]*pow(lambda[i_Qi],5.0)/3.0); 
    Eigen::MatrixXd Q0_i_Qi=Q0[i_Qi]; 
    Q0_i_Qi.resize(1,9); 
    d.row(0)=Q0_i_Qi;     // the first row is Q0[i_Qi]; 
    
    Qi[i_Qi]=d; 
  } 
  return Qi; 
} 



////Construct_GG
// [[Rcpp::export]] 
List Construct_GG(Eigen::VectorXd delta_x, Eigen::VectorXd lambda){  //be careful about delta_x,lambda if one only wants to sample one 
  int num_obs=delta_x.size()+1; 
  int num_dim=lambda.size();  
  List GG(num_dim);                                       // num_dim list, each is 3(num_obs)\times 3 list 
  for(int i_GG=0;i_GG<num_dim;i_GG++){ 
    Eigen::MatrixXd d= Eigen::MatrixXd::Zero(num_obs,9);  //the first row has all zeros  
    for(int j_GG=0;j_GG<(num_obs-1);j_GG++){ 
      int j_GG_1=j_GG+1;    
      d(j_GG_1,0)=pow(lambda[i_GG],2.0)*pow(delta_x[j_GG],2.0)+2*lambda[i_GG]*delta_x[j_GG]+2; 
      d(j_GG_1,1)=-pow(lambda[i_GG],3.0)*pow(delta_x[j_GG],2.0); 
      d(j_GG_1,2)=pow(lambda[i_GG],4.0)*pow(delta_x[j_GG],2.0)-2*pow(lambda[i_GG],3.0)*delta_x[j_GG]; 
      d(j_GG_1,3)=2*(lambda[i_GG]*pow(delta_x[j_GG],2.0)+delta_x[j_GG]); 
      d(j_GG_1,4)=-2*(pow(lambda[i_GG],2.0)*pow(delta_x[j_GG],2.0)-lambda[i_GG]*delta_x[j_GG]-1); 
      d(j_GG_1,5)=2*(pow(lambda[i_GG],3.0)*pow(delta_x[j_GG],2.0)-3*pow(lambda[i_GG],2.0)*delta_x[j_GG]); 
      d(j_GG_1,6)=pow(delta_x[j_GG],2); 
      d(j_GG_1,7)=2*delta_x[j_GG]-lambda[i_GG]*pow(delta_x[j_GG],2.0); 
      d(j_GG_1,8)=pow(lambda[i_GG],2.0)*pow(delta_x[j_GG],2.0)-4*lambda[i_GG]*delta_x[j_GG]+2;     
      d.row(j_GG_1)=exp(-lambda[i_GG]*delta_x[j_GG])/2.0*d.row(j_GG_1); 
    } 
    GG[i_GG]=d; 
  } 
  return GG; 
} 


////update_FRFt 
// [[Rcpp::export]] 
Eigen::VectorXd Update_FRFt(const List dlm_filter_matern_5_2_U_R,const Eigen::MatrixXd dlm_filter_matern_5_2_D_R, const int num_obs){ 
  Eigen::VectorXd FRFt(num_obs); 
  for( int i_t=0;i_t<num_obs;i_t++){ 
    Eigen::MatrixXd dlm_filter_matern_5_2_U_R_i_t=dlm_filter_matern_5_2_U_R[i_t]; 
    FRFt[i_t]=(dlm_filter_matern_5_2_U_R_i_t.row(0).array()*(dlm_filter_matern_5_2_D_R.row(i_t).array().pow(2.0))*dlm_filter_matern_5_2_U_R_i_t.row(0).array()).matrix().sum(); 
  } 
  return  FRFt; 
} 



