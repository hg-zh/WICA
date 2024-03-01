####
#Run this file can get the Figures in the appendix
####
library(ggplot2)
library(grid)
library(gridExtra)
library("cowplot")
load("~/Library/Mobile Documents/com~apple~CloudDocs/Research/CCA/JASA_submit_20211024/Code-Data-submit/Simulation/Case1/RidgeCV/RidgeCV0.5.RData")
M<-matrix(0,nrow = 200,ncol = 4)
for (r in 1:200) {
  M[r,]<-Res[[r]]$Ep
}
df <- data.frame(
  SelectEP= M[,1]
)
HN50<-ggplot(df, aes(x=SelectEP)) + geom_histogram(aes(y=..count../sum(..count..)),bins = 20)+ labs(title = "n=50")+ theme(plot.title = element_text(hjust = 0.5))+ scale_x_log10() + xlab(expression(paste("Selected ",epsilon)))+ylab("percentage")+ theme_bw()
df <- data.frame(
  SelectEP= M[,2]
)
HN100<-ggplot(df, aes(x=SelectEP)) + geom_histogram(aes(y=..count../sum(..count..)),bins = 20)+ labs(title = "n=100")+ theme(plot.title = element_text(hjust = 0.5))+ scale_x_log10()  + xlab(expression(paste("Selected ",epsilon)))+ylab("percentage")+ theme_bw()
df <- data.frame(
  SelectEP= M[,3]
)
HN200<-ggplot(df, aes(x=SelectEP)) + geom_histogram(aes(y=..count../sum(..count..)),bins = 20)+ labs(title = "n=200")+ theme(plot.title = element_text(hjust = 0.5))+ scale_x_log10()  + xlab(expression(paste("Selected ",epsilon)))+ylab("percentage")+ theme_bw()
df <- data.frame(
  SelectEP= M[,4]
)
HN500<-ggplot(df, aes(x=SelectEP)) + geom_histogram(aes(y=..count../sum(..count..)),bins = 20)+ labs(title = "n=500")+ theme(plot.title = element_text(hjust = 0.5))+ scale_x_log10()  + xlab(expression(paste("Selected ",epsilon)))+ylab("percentage")+ theme_bw()
#grid.arrange(N50, N100, N200, N500, ncol=2)
####
load("~/Library/Mobile Documents/com~apple~CloudDocs/Research/CCA/JASA_submit_20211024/Code-Data-submit/Simulation/Case1/KCV/K0.5CV.RData")
SelectEP<-as.vector(M)
N =c(rep('N=50',200),rep('N=100',200),rep('N=200',200),rep('N=500',200))
Data.df<-data.frame(SelectEP,N)
Data.df$N <- factor(Data.df$N , levels=c('N=50','N=100','N=200','N=500'))
ggplot(Data.df, aes(x=N, y=SelectEP)) + geom_boxplot()+ scale_y_log10()
####K
M<-matrix(0,nrow = 200,ncol = 4)
for (r in 1:200) {
  M[r,]<-Res[[r]]$SelectM
}
df <- data.frame(
  SelectM= M[,1]
)
HN50E<-ggplot(df, aes(x=SelectM)) + geom_histogram(aes(y=..count../sum(..count..)),bins = 4)+ labs(title = "n=50")+ theme(plot.title = element_text(hjust = 0.5))+ xlab(expression("Selected k"))+ylab("percentage")+ theme_bw()
df <- data.frame(
  SelectM= M[,2]
)
HN100E<-ggplot(df, aes(x=SelectM)) + geom_histogram(aes(y=..count../sum(..count..)),bins = 4)+ labs(title = "n=100")+ theme(plot.title = element_text(hjust = 0.5))+ xlab(expression("Selected k"))+ylab("percentage")+ theme_bw()
df <- data.frame(
  SelectM= M[,3]
)
HN200E<-ggplot(df, aes(x=SelectM)) + geom_histogram(aes(y=..count../sum(..count..)),bins = 4)+ labs(title = "n=200")+ theme(plot.title = element_text(hjust = 0.5))+ xlab(expression("Selected k"))+ylab("percentage")+ theme_bw()
df <- data.frame(
  SelectM= M[,4]
)
HN500E<-ggplot(df, aes(x=SelectM)) + geom_histogram(aes(y=..count../sum(..count..)),bins = 4)+ labs(title = "n=500")+ theme(plot.title = element_text(hjust = 0.5))+ xlab(expression("Selected k"))+ylab("percentage")+ theme_bw()
grid.arrange(HN50E,HN50, HN100E,HN100, HN200E,HN200, HN500E, HN500, ncol=2)
####
load("~/Library/Mobile Documents/com~apple~CloudDocs/Research/CCA/JASA_submit_20211024/Code-Data-submit/Simulation/Case1/RidgeM/Ridge0.5.RData")
epsilon<-10^rep(seq(-8,-1,length=41),times=4)
Sample_size<-rep(c('n=50','n=100','n=200','n=500'),each=41)
rho<-as.vector(t(E[2:5,]))
df_rho<-data.frame(epsilon=epsilon,Sample_size=Sample_size,rho=rho )
df_rho_new<-df_rho
df_rho_new$Sample_size<- factor(df_rho_new$Sample_size, levels = c('n=50','n=100','n=200','n=500'))
R_rho<-ggplot(data = df_rho_new, mapping = aes(x = epsilon, y = rho, group = Sample_size)) + geom_line(aes(linetype=Sample_size, color=Sample_size))+xlab(expression(epsilon))+ylab(expression(paste(" Error of ",rho)))+ scale_y_log10()+scale_x_log10()+theme(legend.position="top")+ guides(linetype = guide_legend(nrow = 1))
R_rho1<-ggplot(data = df_rho_new, mapping = aes(x = epsilon, y = rho, group = Sample_size)) + geom_line(aes(linetype=Sample_size, color=Sample_size))+xlab(expression(epsilon))+ylab('Absolute Error')+ scale_y_log10()+scale_x_log10()+
  theme(legend.key.height= unit(1, 'cm'), legend.key.width= unit(2, 'cm'),legend.text = element_text( size = 16))+ guides(linetype = guide_legend(nrow = 1))
#####
F<-as.vector(t(M11[2:5,]))
df_rho<-data.frame(epsilon=epsilon,Sample_size=Sample_size,F=F )
df_rho_new<-df_rho
df_rho_new$Sample_size<- factor(df_rho_new$Sample_size, levels = c('n=50','n=100','n=200','n=500'))
R_F<-ggplot(data = df_rho_new, mapping = aes(x = epsilon, y = F, group = Sample_size)) + geom_line(aes(linetype=Sample_size, color=Sample_size))+xlab(expression(epsilon))+ylab('IMSE of U')+ scale_y_log10()+scale_x_log10()+theme(legend.position="top")+ guides(linetype = guide_legend(nrow = 1))
#####
G<-as.vector(t(M22[2:5,]))
df_rho<-data.frame(epsilon=epsilon,Sample_size=Sample_size,G=G )
df_rho_new<-df_rho
df_rho_new$Sample_size<- factor(df_rho_new$Sample_size, levels = c('n=50','n=100','n=200','n=500'))
R_G<-ggplot(data = df_rho_new, mapping = aes(x = epsilon, y = G, group = Sample_size)) + geom_line(aes(linetype=Sample_size, color=Sample_size))+xlab(expression(epsilon))+ylab('IMSE of V')+ scale_y_log10()+scale_x_log10()+theme(legend.position="top")+ guides(linetype = guide_legend(nrow = 1))

#####
load("~/Library/Mobile Documents/com~apple~CloudDocs/Research/CCA/JASA_submit_20211024/Code-Data-submit/Simulation/Case1/KM/K0.5.RData")
epsilon<-rep(1:5,times=4)
Sample_size<-rep(c('n=50','n=100','n=200','n=500'),each=5)
rho<-as.vector(t(E[1:4,]))
df_rho<-data.frame(epsilon=epsilon,Sample_size=Sample_size,rho=rho )
df_rho_new<-df_rho
df_rho_new$Sample_size<- factor(df_rho_new$Sample_size, levels = c('n=50','n=100','n=200','n=500'))
FPC_rho<-ggplot(data = df_rho_new, mapping = aes(x = epsilon, y = rho, group = Sample_size)) + geom_line(aes(linetype=Sample_size, color=Sample_size))+xlab('k')+ylab(expression(paste(" Error of ",rho)))+ scale_y_log10()+theme(legend.position="top")+ guides(linetype = guide_legend(nrow = 1))
#####
F<-as.vector(t(M1[1:4,]))
df_rho<-data.frame(epsilon=epsilon,Sample_size=Sample_size,F=F )
df_rho_new<-df_rho
df_rho_new$Sample_size<- factor(df_rho_new$Sample_size, levels = c('n=50','n=100','n=200','n=500'))
FPC_F<-ggplot(data = df_rho_new, mapping = aes(x = epsilon, y = F, group = Sample_size)) + geom_line(aes(linetype=Sample_size, color=Sample_size))+xlab('k')+ylab('IMSE of U')+ scale_y_log10()+theme(legend.position="top")+ guides(linetype = guide_legend(nrow = 1))
#####
G<-as.vector(t(M2[1:4,]))
df_rho<-data.frame(epsilon=epsilon,Sample_size=Sample_size,G=G )
df_rho_new<-df_rho
df_rho_new$Sample_size<- factor(df_rho_new$Sample_size, levels = c('n=50','n=100','n=200','n=500'))
FPC_G<-ggplot(data = df_rho_new, mapping = aes(x = epsilon, y = G, group = Sample_size)) + geom_line(aes(linetype=Sample_size, color=Sample_size))+xlab('k')+ylab('IMSE of V')+ scale_y_log10()+theme(legend.position="top")+ guides(linetype = guide_legend(nrow = 1))
########
R_rho<- R_rho + theme(legend.position="none")
R_F<- R_F + theme(legend.position="none")
R_G<- R_G + theme(legend.position="none")
FPC_rho<- FPC_rho + theme(legend.position="none")
FPC_F<- FPC_F + theme(legend.position="none")
FPC_G<- FPC_G + theme(legend.position="none")
blankPlot <- ggplot()+geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()
legend <- get_legend(R_rho1)
grid.arrange(legend, R_rho , R_F , R_G ,FPC_rho,FPC_F, FPC_G,HN50E,HN50, HN100E,HN100, HN200E,HN200, HN500E, HN500,
             ncol=2, nrow = 8, layout_matrix = cbind(c(5,6,7,8,10,12,14,1),c(2,3,4,9,11,13,15,1)),
             widths = c(3, 3), heights = c( 3,3,3,3,3,3,3,1))
