##FPCA method 
##tuning parameters selected by five-fold CV
##'GetResult' is the mean function
GetResult<-function(r){
  library(fdapace)
  library(TruncatedNormal)
  n<-c(50,100,200,500);sigma=0.1;K_s<-5
  M1<-matrix(0,nrow = length(n),ncol =1)
  M2<-matrix(0,nrow = length(n),ncol =1)
  E<-matrix(0,nrow = length(n),ncol =1)
  SelectM<-matrix(0,nrow = length(n),ncol =1)
  F_main<- function(n,r) {
    #n<-50;r=1;sigma=0.1;
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
        Xi_x[i,k]<-(2^-k/(1.78*sqrt(2)*pi*k ))*runif(1,-1,1)
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
      Xi_y[i,1]<-(2^-1/(1.78*sqrt(2)*pi ))*runif(1,-1,1)
      Xi_y[i,2]<-  0.5*(Xi_x[i,1]+Xi_x[i,2])+sigma*sqrt(17)/(1.78*sqrt(2)*pi*8)*runif(1,-1,1)
      for (k in 3:K) {
        Xi_y[i,k]<-(2^-k/(1.78*sqrt(2)*pi*k ))*runif(1,-1,1)
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
    f_true<-sqrt(2)/2*Phi_x[[1]]+ sqrt(2)/2*Phi_x[[2]]
    g_true<-Phi_y[[2]]
    ######CV
    CV<-function(){
      S<-rep(0,5)
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
        Test_Ty_hat<-list()
        for (i in 1:length(Test_y)) {
          temp<-matrix(0,nrow =length(t),ncol = 101 )
          for (j in 1:length(t)) {
            temp[j,]<-approx(seq(0,1,length=101),Test_y[[i]][j,],Fyp_hat[j,],rule = 2)$y-x
          }
          Test_Ty_hat[[i]]<-temp
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
        ###
        for (m in 1:K_s) {
          Ax<-matrix(0,nrow = m,ncol = m)
          Ay<-matrix(0,nrow = m,ncol = m)
          Axy<-matrix(0,nrow = m,ncol = m)
          Ayx<-matrix(0,nrow = m,ncol = m)
          for (i in 1:length(Train_id)) {
            for (j in 1:m) {
              for (k in 1:m) {
                Ax[j,k]<-Ax[j,k]+Xix_hat[i,j]*Xix_hat[i,k]/length(Train_id)
                Ay[j,k]<-Ay[j,k]+Xiy_hat[i,j]*Xiy_hat[i,k]/length(Train_id)
                Ayx[j,k]<-Ayx[j,k]+Xiy_hat[i,j]*Xix_hat[i,k]/length(Train_id)
                Axy[j,k]<-Axy[j,k]+Xix_hat[i,j]*Xiy_hat[i,k]/length(Train_id)
              }
            }
          }
          C<-solve(Ax)%*%Axy%*%solve(Ay)%*%Ayx
          R<-Re(eigen(C)$vectors[,1])
          g<-solve(Ay)%*%Ayx%*%R/norm((solve(Ay))%*%Ayx%*%R,type = "2")
          f_hat<-matrix(0,nrow = 101,ncol = 101)
          g_hat<-matrix(0,nrow = 101,ncol = 101)
          for (k in 1:m) {
            f_hat<-f_hat+R[k]*Phix_hat[[k]]
            g_hat<-g_hat+g[k]*Phiy_hat[[k]]
          }
          S1<-0;S2<-0
          for (s in 1:99) {
            S1<-S1+sum((f_true[s,]-f_hat[s,])^2*dx_hat[s,] )*0.01*0.01
            S2<-S2+sum((f_true[s,]+f_hat[s,])^2*dx_hat[s,] )*0.01*0.01
          }
          if(S1>=S2){f_hat<--f_hat}
          S1<-0;S2<-0
          for (s in 1:99) {
            S1<-S1+sum((g_true[s,]-g_hat[s,])^2*dx_hat[s,] )*0.01*0.01
            S2<-S2+sum((g_true[s,]+g_hat[s,])^2*dx_hat[s,] )*0.01*0.01
          }
          if(S1>=S2){g_hat<--g_hat}
          A_x<-matrix(0,nrow =length(Test_x),ncol = 1 )
          A_y<-matrix(0,nrow =length(Test_x),ncol = 1 )
          for (i in 1:length(Test_x)) {
            
            temp<-0
            for (j in 1:length(t)) {
              temp<-temp+sum(Test_Tx_hat[[i]][j,]*f_hat[j,]*dx_hat[j,])
            }
            A_x[i,1]<-temp*0.01*0.01
            temp<-0
            for (j in 1:length(t)) {
              temp<-temp+sum(Test_Ty_hat[[i]][j,]*g_hat[j,]*dy_hat[j,])
            }
            A_y[i,1]<-temp*0.01*0.01
          }
          S[m]<-S[m]+cor(A_x,A_y)^2
        }
      }
      return(which.max(S))
    }
    m<-CV()
    ######Cov Estmation#####
    Ax<-matrix(0,nrow = m,ncol = m)
    Ay<-matrix(0,nrow = m,ncol = m)
    Axy<-matrix(0,nrow = m,ncol = m)
    Ayx<-matrix(0,nrow = m,ncol = m)
    for (i in 1:n) {
      for (j in 1:m) {
        for (k in 1:m) {
          Ax[j,k]<-Ax[j,k]+Xix_hat[i,j]*Xix_hat[i,k]/n
          Ay[j,k]<-Ay[j,k]+Xiy_hat[i,j]*Xiy_hat[i,k]/n
          Ayx[j,k]<-Ayx[j,k]+Xiy_hat[i,j]*Xix_hat[i,k]/n
          Axy[j,k]<-Axy[j,k]+Xix_hat[i,j]*Xiy_hat[i,k]/n
        }
      }
    }
    C<-solve(Ax)%*%Axy%*%solve(Ay)%*%Ayx
    R<-Re(eigen(C)$vectors[,1])
    g<-solve(Ay)%*%Ayx%*%R/norm((solve(Ay))%*%Ayx%*%R,type = "2")
    out2<-abs(Re(eigen(C)$values[1])-0.5^2/(0.5^2+sigma^2))
    f_hat<-matrix(0,nrow = 101,ncol = 101)
    g_hat<-matrix(0,nrow = 101,ncol = 101)
    for (k in 1:m) {
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
    return(list(f=out1,g=out3,Corr=out2,Selm=m))
  }
  for (j in 1:length(n)) {
    temp<-F_main(n[j],r)
    M1[j,]<-temp$f
    M2[j,]<-temp$g
    E[j,]<-temp$Corr
    SelectM[j,]<-temp$Selm
  }
  return(list(f=M1,g=M2,Corr=E,SelectM=SelectM))
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
M