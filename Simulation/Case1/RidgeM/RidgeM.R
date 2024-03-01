##Tikhonov method 
##'GetResult' is the mean function
## 'sigma' is the tunning parameter
GetResult<-function(r){
  library(TruncatedNormal)
  library(fdapace)
  n<-c(10,50,100,200,500,1000);epsilon<-10^seq(-8,-1,length.out = 41);sigma=0.3
  M1<-matrix(0,nrow = length(n),ncol =length(epsilon))
  M2<-matrix(0,nrow = length(n),ncol =length(epsilon))
  M11<-matrix(0,nrow = length(n),ncol =length(epsilon))
  M22<-matrix(0,nrow = length(n),ncol =length(epsilon))
  E<-matrix(0,nrow = length(n),ncol =length(epsilon))
  F_main<- function(n,r,epsilon) {
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
    out2<-c();out1<-c();out3<-c();out11<-c();out33<-c()
    for (j in 1:length(epsilon)) {
      C<-solve(Ax+epsilon[j]*diag(M))%*%Axy%*%solve(Ay+epsilon[j]*diag(M))%*%Ayx
      R<-Re(eigen(C)$vectors[,1])
      out2[j]<-abs(Re(eigen(C)$values[1])- 0.5^2/(0.5^2+sigma^2))
      g<-solve(Ay+epsilon[j]*diag(M))%*%Ayx%*%R/norm((solve(Ay+epsilon[j]*diag(M)))%*%Ayx%*%R,type = "2")
      ####
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
      out11[j]<-sqrt(min(S1,S2))
      S1<-0;S2<-0
      for (s in 1:99) {
        S1<-S1+sum((g_true[s,]-g_hat[s,])^2*dx_hat[s,] )*0.01*0.01
        S2<-S2+sum((g_true[s,]+g_hat[s,])^2*dx_hat[s,] )*0.01*0.01
      }
      out33[j]<-sqrt(min(S1,S2))
    }
    return(list(f=out1,g=out3,Corr=out2,f_1=out11,g_1=out33 ))
  }
  for (j in 1:length(n)) {
    temp<-F_main(n[j],r,epsilon)
    E[j,]<-temp$Corr
    M11[j,]<-temp$f_1
    M22[j,]<-temp$g_1
  }
  return(list(Corr=E,f_1=M11,g_1=M22))
}
#######
library(parallel)
system.time({
  cl <- makeCluster(getOption("cl.cores", 100))
  Res<-parLapply(cl, 1:200,GetResult)
  stopCluster(cl)
  remove(cl)
})
E<-matrix(0,nrow =6,ncol = 41 )
M11<-matrix(0,nrow =6,ncol = 41 )
M22<-matrix(0,nrow =6,ncol = 41 )
for (r in 1:200) {
  E<-E+Res[[r]]$Corr/200
  M11<-M11+Res[[r]]$f_1/200
  M22<-M22+Res[[r]]$g_1/200
}
#######
library(ggplot2)
epsilon<-10^rep(seq(-8,-1,length=41),times=6)
Sample_size<-rep(c('n=10','n=50','n=100','n=200','n=500','n=1000'),each=41)
rho<-as.vector(t(E))
df_rho<-data.frame(epsilon=epsilon,Sample_size=Sample_size,rho=rho )
ggplot(data = df_rho, mapping = aes(x = epsilon, y = rho, colour = Sample_size)) + geom_line()+ 
  scale_colour_discrete(name="Sample Size",breaks = c('n=10','n=50','n=100','n=200','n=500','n=1000'))+xlab('Epsilon')+ylab('Absolute Error')+ labs(title = expression(sigma~'=0.5'))+theme(plot.title = element_text(hjust = 0.5)) +scale_y_log10()+scale_x_log10()
#####
F<-as.vector(t(M11))
df_rho<-data.frame(epsilon=epsilon,Sample_size=Sample_size,F=F )
ggplot(data = df_rho, mapping = aes(x = epsilon, y = F, colour = Sample_size)) + geom_line()+ 
  scale_colour_discrete(name="Sample Size",breaks = c('n=10','n=50','n=100','n=200','n=500','n=1000'))+xlab('Epsilon')+ylab('IMSE of U')+ labs(title = expression(sigma~'=0.5'))+theme(plot.title = element_text(hjust = 0.5)) +scale_y_log10()+scale_x_log10()
#####
G<-as.vector(t(M22))
df_rho<-data.frame(epsilon=epsilon,Sample_size=Sample_size,G=G )
ggplot(data = df_rho, mapping = aes(x = epsilon, y = G, colour = Sample_size)) + geom_line()+ 
  scale_colour_discrete(name="Sample Size",breaks = c('n=10','n=50','n=100','n=200','n=500','n=1000'))+xlab('Epsilon')+ylab('IMSE of V')+ labs(title = expression(sigma~'=0.5'))+theme(plot.title = element_text(hjust = 0.5)) +scale_y_log10()+scale_x_log10()
##############

