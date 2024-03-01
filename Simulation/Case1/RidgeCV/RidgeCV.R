##Tikhonov method 
##tuning parameters selected by five-fold CV
##'GetResult' is the mean function
## 'sigma' is the tunning parameter
GetResult<-function(r){
  library(fdapace)
  library(TruncatedNormal)
  n<-c(50,100,200,500);sigma=0.5
  M1<-matrix(0,nrow = length(n),ncol =1)
  M2<-matrix(0,nrow = length(n),ncol =1)
  E<-matrix(0,nrow = length(n),ncol =1)
  epsilon<-matrix(0,nrow = length(n),ncol =1)
  F_main<- function(n,r) {
    #n<-50;r=1;sigma=0.1;fold=1
    set.seed(r)
    K<-20;M<-20
    Get_F<-function(x,a,b){
      return(approx(a,b,x,rule=2)$y)
    }
    Get_Fq<-function(Data_T,F_q,t){
      temp<-list()
      for (i in 1:length(Data_T)) {
        temp[[i]]<-F_q+t^{-1}*Data_T[[i]]
      }
      return(temp)
    }
    #####Generate X#####
    f_xd<-function(t,x){
      n<-length(x)
      out<-c()
      for (i in 1:n) {
        out[i]<-dbeta(x[i],(2+t[i]),(3-t[i]^2/2-t[i]/2))
      }
      return(out)
    }
    f_xq<-function(t,x){
      n<-length(x)
      out<-c()
      for (i in 1:n) {
        out[i]<-qbeta(x[i],(2+t[i]),(3-t[i]^2/2-t[i]/2))
      }
      return(out)
    }
    f_xp<-function(t,x){
      n<-length(x)
      out<-c()
      for (i in 1:n) {
        out[i]<-pbeta(x[i],(2+t[i]),(3-t[i]^2/2-t[i]/2))
      }
      return(out)
    }
    t<-seq(0,1,length=101)
    x<-t
    F_xd<-outer(t,x,f_xd)
    F_xq<-outer(t,x,f_xq)
    F_xp<-outer(t,x,f_xp)
    Basis_sin<-list()
    for (i in 1:K) {
      temp<-matrix(0,nrow = length(t),ncol = length(x))
      for (j in  1:length(t)) {
        temp[j,]<-CreateBasis(K,t,type = 'sin')[,i]
      }
      Basis_sin[[i]]<-temp
      remove(temp)
    }
    Phi_x<-list()
    for (i in 1:K) {
      temp<-matrix(0,nrow = length(t),ncol = length(x))
      for (j in 1:length(t)) {
        temp[j,]<-CreateBasis(K,F_xp[j,],type = 'sin')[,i]
      }
      Phi_x[[i]]<-temp
      remove(temp)
    }
    Xi_x<-matrix(0,nrow = n,ncol = K)
    for (i in 1:n) {
      for (k in 1:K) {
        Xi_x[i,k]<-(2^-k/(1.78*sqrt(2)*pi*k ))*rtnorm(n = 1, mu = 0, lb = -1, ub = 1, method = "fast",sd=1)
      }
    }
    Data_Tx<-list()
    for (i in 1:n) {
      Data_Tx[[i]]<-matrix(0,nrow =length(t),ncol = length(x) )
      for (k in 1:K) {
        Data_Tx[[i]]<-Data_Tx[[i]]+Xi_x[i,k]*Basis_sin[[k]]
      }
    }
    Data_Fqx<-Get_Fq(Data_Tx,F_xq,1)
    #####Generate Y#####
    f_yd<-function(t,x){
      n<-length(x)
      out<-c()
      for (i in 1:n) {
        out[i]<-dbeta(x[i],(3-t[i]),(t[i]^2/2+t[i]/2+2 ))
      }
      return(out)
    }
    f_yq<-function(t,x){
      n<-length(x)
      out<-c()
      for (i in 1:n) {
        out[i]<-qbeta(x[i],(3-t[i]),(t[i]^2/2+t[i]/2+2 ))
      }
      return(out)
    }
    f_yp<-function(t,x){
      n<-length(x)
      out<-c()
      for (i in 1:n) {
        out[i]<-pbeta(x[i],(3-t[i]),(t[i]^2/2+t[i]/2+2 ))
      }
      return(out)
    }
    F_yd<-outer(t,x,f_yd)
    F_yp<-outer(t,x,f_yp)
    F_yq<-outer(t,x,f_yq)
    Phi_y<-list()
    for (i in 1:K) {
      temp<-matrix(0,nrow = length(t),ncol = length(x))
      for (j in 1:length(t)) {
        temp[j,]<-CreateBasis(K,F_yp[j,],type = 'sin')[,i]
      }
      Phi_y[[i]]<-temp
      remove(temp)
    }
    Xi_y<-matrix(0,nrow = n,ncol = K)
    for (i in 1:n) {
      Xi_y[i,1]<-(2^-1/(1.78*sqrt(2)*pi ))*rtnorm(n = 1, mu = 0, lb = -1, ub = 1, method = "fast",sd=1)
      Xi_y[i,2]<-  0.5*(Xi_x[i,1]+Xi_x[i,2])+sigma*sqrt(17)/(1.78*sqrt(2)*pi*8)*rtnorm(n = 1, mu = 0, lb = -1, ub = 1, method = "fast",sd=1)
      for (k in 3:K) {
        Xi_y[i,k]<-(2^-k/(1.78*sqrt(2)*pi*k ))*rtnorm(n = 1, mu = 0, lb = -1, ub = 1, method = "fast",sd=1)
      }
    }
    Data_Ty<-list()
    for (i in 1:n) {
      Data_Ty[[i]]<-matrix(0,nrow =length(t),ncol = length(x) )
      for (k in 1:K) {
        Data_Ty[[i]]<-Data_Ty[[i]]+Xi_y[i,k]*Basis_sin[[k]]
      }
    }
    Data_Fqy<-Get_Fq(Data_Ty,F_yq,1)
    ######Mean Estimation#####
    ##X
    Fxq_hat<-matrix(0,nrow = length(t),ncol = length(x))
    for (i in 1:n) {
      Fxq_hat<-Fxq_hat+Data_Fqx[[i]]/n
    }
    Fxp_hat<-matrix(0,nrow = length(t),ncol = length(x))
    for (j in 1:length(t)) {
      for (k in 2:100) {
        Fxp_hat[j,k]<-uniroot((function(x) Get_F(x,a=seq(0,1,length=101),b=Fxq_hat[j,])-seq(0,1,length=101)[k]  ),lower = 0,upper = 1)[1]$root
      }
      Fxp_hat[j,101]<-1
    }
    dx_hat<-matrix(0,nrow = 101,ncol = 101)
    for (j in 1:101) {
      for (s in 1:100) {
        dx_hat[j,s]<-(Fxp_hat[j,s+1]-Fxp_hat[j,s])*100
      }
    }
    ##Y
    Fyq_hat<-matrix(0,nrow = length(t),ncol = length(x))
    for (i in 1:n) {
      Fyq_hat<-Fyq_hat+Data_Fqy[[i]]/n
    }
    Fyp_hat<-matrix(0,nrow = length(t),ncol = length(x))
    for (j in 1:length(t)) {
      for (k in 2:100) {
        Fyp_hat[j,k]<-uniroot((function(x) Get_F(x,a=seq(0,1,length=101),b=Fyq_hat[j,])-seq(0,1,length=101)[k]  ),lower = 0,upper = 1)[1]$root
      }
      Fyp_hat[j,101]<-1
    }
    dy_hat<-matrix(0,nrow = 101,ncol = 101)
    for (j in 1:101) {
      for (s in 1:100) {
        dy_hat[j,s]<-(Fyp_hat[j,s+1]-Fyp_hat[j,s])*100
      }
    }
    #####Log-mapping#####
    Data_Tx_hat<-list()
    for (i in 1:n) {
      temp<-matrix(0,nrow =length(t),ncol = 101 )
      for (j in 1:length(t)) {
        temp[j,]<-approx(seq(0,1,length=101),Data_Fqx[[i]][j,],Fxp_hat[j,],rule = 2)$y-x
      }
      Data_Tx_hat[[i]]<-temp
      rm(temp)
    }
    Data_Ty_hat<-list()
    for (i in 1:n) {
      temp<-matrix(0,nrow =length(t),ncol = 101 )
      for (j in 1:length(t)) {
        temp[j,]<-approx(seq(0,1,length=101),Data_Fqy[[i]][j,],Fyp_hat[j,],rule = 2)$y-x
      }
      Data_Ty_hat[[i]]<-temp
      rm(temp)
    }
    #####Projection#####
    ####X
    Phix_hat<-list()
    for (i in 1:M) {
      temp<-matrix(0,nrow = length(t),ncol = length(x))
      for (j in 1:length(t)) {
        temp[j,]<-CreateBasis(M,Fxp_hat[j,],type = 'sin')[,i]
      }
      Phix_hat[[i]]<-temp
      remove(temp)
    }
    Xix_hat<-matrix(0, nrow = n, ncol = M)
    for (i in 1:n) {
      for (k in 1:M) {
        temp<-0
        for (j in 1:length(t)) {
          temp<-temp+sum(Data_Tx_hat[[i]][j,]*Phix_hat[[k]][j,]*dx_hat[j,])
        }
        Xix_hat[i,k]<-temp*0.01*0.01
        remove(temp)
      }
    }
    ####y
    Phiy_hat<-list()
    for (i in 1:M) {
      temp<-matrix(0,nrow = length(t),ncol = length(x))
      for (j in 1:length(t)) {
        temp[j,]<-CreateBasis(M,Fyp_hat[j, ],type = 'sin')[,i]
      }
      Phiy_hat[[i]]<-temp
      remove(temp)
    }
    Xiy_hat<-matrix(0, nrow = n, ncol = M)
    for (i in 1:n) {
      for (k in 1:M) {
        temp<-0
        for (j in 1:length(t)) {
          temp<-temp+sum(Data_Ty_hat[[i]][j,]*Phiy_hat[[k]][j,]*dy_hat[j,])
        }
        Xiy_hat[i,k]<-temp*0.01*0.01
        remove(temp)
      }
    }
    ######
    SSplit<-function(){
      Train_ax<-list();Train_ay<-list();Train_axy<-list();Train_ayx<-list();
      Train_phix<-list();Train_phiy<-list()
      Train_dx<-list();Train_dy<-list()
      Test_Tx<-list();Test_Ty<-list()
      for (fold in 1:5) {
        Train_x<-list();Train_y<-list()
        Test_x<-list();Test_y<-list()
        Test_id<-seq((fold-1)*n/5+1,fold*n/5,1)
        Train_id<-setdiff(1:n,Test_id)
        for (i in 1:length(Train_id) ) {
          Train_x[[i]]<-Data_Fqx[[Train_id[i]]]
          Train_y[[i]]<-Data_Fqy[[Train_id[i]]]
        }
        for (i in 1:length(Test_id) ) {
          Test_x[[i]]<-Data_Fqx[[Test_id[i]]]
          Test_y[[i]]<-Data_Fqy[[Test_id[i]]]
        }
        ######Mean Estimation#####
        ##X
        Fxq_hat<-matrix(0,nrow = length(t),ncol = length(x))
        for (i in 1:length(Train_x)) {
          Fxq_hat<-Fxq_hat+Train_x[[i]]/length(Train_x)
        }
        Fxp_hat<-matrix(0,nrow = length(t),ncol = length(x))
        for (j in 1:length(t)) {
          for (k in 2:100) {
            Fxp_hat[j,k]<-uniroot((function(x) Get_F(x,a=seq(0,1,length=101),b=Fxq_hat[j,])-seq(0,1,length=101)[k]  ),lower = 0,upper = 1)[1]$root
          }
          Fxp_hat[j,101]<-1
        }
        dx_hat<-matrix(0,nrow = 101,ncol = 101)
        for (j in 1:101) {
          for (s in 1:100) {
            dx_hat[j,s]<-(Fxp_hat[j,s+1]-Fxp_hat[j,s])*100
          }
        }
        Train_dx[[fold]]<-dx_hat
        ##Y
        Fyq_hat<-matrix(0,nrow = length(t),ncol = length(x))
        for (i in 1:length(Train_y)) {
          Fyq_hat<-Fyq_hat+Train_y[[i]]/length(Train_y)
        }
        Fyp_hat<-matrix(0,nrow = length(t),ncol = length(x))
        for (j in 1:length(t)) {
          for (k in 2:100) {
            Fyp_hat[j,k]<-uniroot((function(x) Get_F(x,a=seq(0,1,length=101),b=Fyq_hat[j,])-seq(0,1,length=101)[k]  ),lower = 0,upper = 1)[1]$root
          }
          Fyp_hat[j,101]<-1
        }
        dy_hat<-matrix(0,nrow = 101,ncol = 101)
        for (j in 1:101) {
          for (s in 1:100) {
            dy_hat[j,s]<-(Fyp_hat[j,s+1]-Fyp_hat[j,s])*100
          }
        }
        Train_dy[[fold]]<-dy_hat
        #####Log-mapping#####
        Data_Tx_hat<-list()
        for (i in 1:length(Train_x)) {
          temp<-matrix(0,nrow =length(t),ncol = 101 )
          for (j in 1:length(t)) {
            temp[j,]<-approx(seq(0,1,length=101),Train_x[[i]][j,],Fxp_hat[j,],rule = 2)$y-x
          }
          Data_Tx_hat[[i]]<-temp
          rm(temp)
        }
        Data_Ty_hat<-list()
        for (i in 1:length(Train_y)) {
          temp<-matrix(0,nrow =length(t),ncol = 101 )
          for (j in 1:length(t)) {
            temp[j,]<-approx(seq(0,1,length=101),Train_y[[i]][j,],Fyp_hat[j,],rule = 2)$y-x
          }
          Data_Ty_hat[[i]]<-temp
          rm(temp)
        }
        Test_Tx_hat<-list()
        for (i in 1:length(Test_x)) {
          temp<-matrix(0,nrow =length(t),ncol = 101 )
          for (j in 1:length(t)) {
            temp[j,]<-approx(seq(0,1,length=101),Test_x[[i]][j,],Fxp_hat[j,],rule = 2)$y-x
          }
          Test_Tx_hat[[i]]<-temp
          rm(temp)
        }
        Test_Tx[[fold]]<-Test_Tx_hat
        Test_Ty_hat<-list()
        for (i in 1:length(Test_y)) {
          temp<-matrix(0,nrow =length(t),ncol = 101 )
          for (j in 1:length(t)) {
            temp[j,]<-approx(seq(0,1,length=101),Test_y[[i]][j,],Fyp_hat[j,],rule = 2)$y-x
          }
          Test_Ty_hat[[i]]<-temp
          rm(temp)
        }
        Test_Ty[[fold]]<-Test_Ty_hat
        #####Projection#####
        ####X
        Phix_hat<-list()
        for (i in 1:M) {
          temp<-matrix(0,nrow = length(t),ncol = length(x))
          for (j in 1:length(t)) {
            temp[j,]<-CreateBasis(M,Fxp_hat[j,],type = 'sin')[,i]
          }
          Phix_hat[[i]]<-temp
          remove(temp)
        }
        Train_phix[[fold]]<-Phix_hat
        Xix_hat<-matrix(0, nrow = length(Train_x), ncol = M)
        for (i in 1:length(Train_x)) {
          for (k in 1:M) {
            temp<-0
            for (j in 1:length(t)) {
              temp<-temp+sum(Data_Tx_hat[[i]][j,]*Phix_hat[[k]][j,]*dx_hat[j,])
            }
            Xix_hat[i,k]<-temp*0.01*0.01
            remove(temp)
          }
        }
        ####y
        Phiy_hat<-list()
        for (i in 1:M) {
          temp<-matrix(0,nrow = length(t),ncol = length(x))
          for (j in 1:length(t)) {
            temp[j,]<-CreateBasis(M,Fyp_hat[j, ],type = 'sin')[,i]
          }
          Phiy_hat[[i]]<-temp
          remove(temp)
        }
        Train_phiy[[fold]]<-Phiy_hat
        Xiy_hat<-matrix(0, nrow = length(Train_y), ncol = M)
        for (i in 1:length(Train_y)) {
          for (k in 1:M) {
            temp<-0
            for (j in 1:length(t)) {
              temp<-temp+sum(Data_Ty_hat[[i]][j,]*Phiy_hat[[k]][j,]*dy_hat[j,])
            }
            Xiy_hat[i,k]<-temp*0.01*0.01
            remove(temp)
          }
        }
        Ax<-matrix(0,nrow = M,ncol = M)
        Ay<-matrix(0,nrow = M,ncol = M)
        Axy<-matrix(0,nrow = M,ncol = M)
        Ayx<-matrix(0,nrow = M,ncol = M)
        for (i in 1:length(Train_id)) {
          for (j in 1:M) {
            for (k in 1:M) {
              Ax[j,k]<-Ax[j,k]+Xix_hat[i,j]*Xix_hat[i,k]/length(Train_id)
              Ay[j,k]<-Ay[j,k]+Xiy_hat[i,j]*Xiy_hat[i,k]/length(Train_id)
              Ayx[j,k]<-Ayx[j,k]+Xiy_hat[i,j]*Xix_hat[i,k]/length(Train_id)
              Axy[j,k]<-Axy[j,k]+Xix_hat[i,j]*Xiy_hat[i,k]/length(Train_id)
            }
          }
        }
        Train_ax[[fold]]<-Ax;Train_ay[[fold]]<-Ay
        Train_axy[[fold]]<-Axy;Train_ayx[[fold]]<-Ayx
      }
      return(list(Train_ax=Train_ax,Train_ay=Train_ay,Train_axy=Train_axy,Train_ayx=Train_ayx,Train_phix=Train_phix,Train_phiy=Train_phiy,Train_dx=Train_dx,Train_dy=Train_dy,Test_Tx=Test_Tx,Test_Ty=Test_Ty))
    }
    res<-SSplit()
    CV<-function(epsilon){
      lgep<-length(epsilon)
      S<-rep(0,lgep)
      for (e in 1:lgep) {
        for (fold in 1:5) {
          Test_x<-list();Test_y<-list()
          Test_id<-seq((fold-1)*n/5+1,fold*n/5,1)
          Train_id<-setdiff(1:n,Test_id)
          for (i in 1:length(Test_id) ) {
            Test_x[[i]]<-Data_Fqx[[Test_id[i]]]
            Test_y[[i]]<-Data_Fqy[[Test_id[i]]]
          }
          C<-solve(res$Train_ax[[fold]]+epsilon[e]*diag(M))%*%res$Train_axy[[fold]]%*%solve(res$Train_ay[[fold]]+epsilon[e]*diag(M))%*%res$Train_ayx[[fold]]
          R<-Re(eigen(C)$vectors[,1])
          if(sum((R-c(sqrt(2)/2,sqrt(2)/2,rep(0,K-2)))^2)>=sum((R+c(sqrt(2)/2,sqrt(2)/2,rep(0,K-2)))^2)){R<--R}
          G<-solve(res$Train_ay[[fold]]+epsilon[e]*diag(M))%*%res$Train_ayx[[fold]]%*%R/norm((solve(res$Train_ay[[fold]]+epsilon[e]*diag(M)))%*%res$Train_ayx[[fold]]%*%R,type = "2")
          f_hat<-matrix(0,nrow = 101,ncol = 101)
          g_hat<-matrix(0,nrow = 101,ncol = 101)
          for (k in 1:M) {
            f_hat<-f_hat+R[k]*res$Train_phix[[fold]][[k]]
            g_hat<-g_hat+G[k]*res$Train_phiy[[fold]][[k]]
          }
          A_x<-matrix(0,nrow =length(Test_x),ncol = 1 )
          A_y<-matrix(0,nrow =length(Test_x),ncol = 1 )
          for (i in 1:length(Test_x)) {
            temp<-0
            for (j in 1:length(t)) {
              temp<-temp+sum(res$Test_Tx[[fold]][[i]][j,]*f_hat[j,]*res$Train_dx[[fold]][j,])
            }
            A_x[i,1]<-temp*0.01*0.01
            temp<-0
            for (j in 1:length(t)) {
              temp<-temp+sum(res$Test_Ty[[fold]][[i]][j,]*g_hat[j,]*res$Train_dx[[fold]][j,])
            }
            A_y[i,1]<-temp*0.01*0.01
          }
          S[e]<-S[e]+cor(A_x,A_y)^2
        }
      }
      
      return(-S)
    }
    epsilon<-optimize(CV,interval = c(10^-8,10^-1),tol = .Machine$double.eps)$minimum
    ######Cov Estmation#####
    Ax<-matrix(0,nrow = M,ncol = M)
    Ay<-matrix(0,nrow = M,ncol = M)
    Axy<-matrix(0,nrow = M,ncol = M)
    Ayx<-matrix(0,nrow = M,ncol = M)
    for (i in 1:n) {
      for (j in 1:M) {
        for (k in 1:M) {
          Ax[j,k]<-Ax[j,k]+Xix_hat[i,j]*Xix_hat[i,k]/n
          Ay[j,k]<-Ay[j,k]+Xiy_hat[i,j]*Xiy_hat[i,k]/n
          Ayx[j,k]<-Ayx[j,k]+Xiy_hat[i,j]*Xix_hat[i,k]/n
          Axy[j,k]<-Axy[j,k]+Xix_hat[i,j]*Xiy_hat[i,k]/n
        }
      }
    }
    f_true<-sqrt(2)/2*Phi_x[[1]]+ sqrt(2)/2*Phi_x[[2]]
    g_true<-Phi_y[[2]]
    C<-solve(Ax+epsilon*diag(M))%*%Axy%*%solve(Ay+epsilon*diag(M))%*%Ayx
    R<-Re(eigen(C)$vectors[,1])
    g<-solve(Ay+epsilon*diag(M))%*%Ayx%*%R/norm((solve(Ay+epsilon*diag(M)))%*%Ayx%*%R,type = "2")
    out2<-abs(Re(eigen(C)$values[1])-0.5^2/(0.5^2+sigma^2))
    f_hat<-matrix(0,nrow = 101,ncol = 101)
    g_hat<-matrix(0,nrow = 101,ncol = 101)
    for (k in 1:M) {
      f_hat<-f_hat+R[k]*Phix_hat[[k]]
      g_hat<-g_hat+g[k]*Phiy_hat[[k]]
    }
    S1<-0;S2<-0
    for (s in 1:99) {
      S1<-S1+sum((f_true[s,]-f_hat[s,])^2*dx_hat[s,] )*0.01*0.01
      S2<-S2+sum((f_true[s,]+f_hat[s,])^2*dx_hat[s,] )*0.01*0.01
    }
    out1<-sqrt(min(S1,S2))
    S1<-0;S2<-0
    for (s in 1:99) {
      S1<-S1+sum((g_true[s,]-g_hat[s,])^2*dx_hat[s,] )*0.01*0.01
      S2<-S2+sum((g_true[s,]+g_hat[s,])^2*dx_hat[s,] )*0.01*0.01
    }
    out3<-sqrt(min(S1,S2))
    return(list(f=out1,g=out3,Corr=out2,Ep=epsilon))
  }
  for (j in 1:length(n)) {
    temp<-F_main(n[j],r)
    M1[j,]<-temp$f
    M2[j,]<-temp$g
    E[j,]<-temp$Corr
    epsilon[j,]<-temp$Ep
  }
  return(list(f=M1,g=M2,Corr=E,Ep=epsilon))
}
system.time({
  cl <- makeCluster(getOption("cl.cores", 100))
  Res<-parLapply(cl, 1:200,GetResult)
  stopCluster(cl)
  remove(cl)
})
########
IMSE_F<-matrix(0,4,1);IMSE_G<-matrix(0,4,1);Abs_Cor<-matrix(0,4,1);Ep<-matrix(0,4,1)
for (r in 1:200) {
  IMSE_F<-IMSE_F+Res[[r]]$f/200
  IMSE_G<-IMSE_G+Res[[r]]$g/200
  Abs_Cor<-Abs_Cor+Res[[r]]$Corr/200
  Ep<-Ep+Res[[r]]$Ep/200
}
IMSE_F
IMSE_G
Abs_Cor
Ep
