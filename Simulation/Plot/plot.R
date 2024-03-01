####
#Run this file can get Figure 2.
####
library(grid)
library(gridExtra)
library("cowplot")
library(ggplot2)
##########
load("~/Documents/research/CCA/Simulation/1_22/RidgeM/Ridge0.05.RData")
epsilon<-10^rep(seq(-8,-1,length=41),times=4)
Sample_size<-rep(c('n=50','n=100','n=200','n=500'),each=41)
rho<-as.vector(t(E[2:5,]))
df_rho<-data.frame(epsilon=epsilon,Sample_size=Sample_size,rho=rho )
df_rho_new<-df_rho
df_rho_new$Sample_size<- factor(df_rho_new$Sample_size, levels = c('n=50','n=100','n=200','n=500'))
R_rho<-ggplot(data = df_rho_new, mapping = aes(x = epsilon, y = rho, group = Sample_size)) + geom_line(aes(linetype=Sample_size, color=Sample_size))+xlab(expression(epsilon))+ylab(expression(paste("Absolute Error of ",rho)))+ scale_y_log10()+scale_x_log10()+theme(legend.position="top")+ guides(linetype = guide_legend(nrow = 1))
R_rho1<-ggplot(data = df_rho_new, mapping = aes(x = epsilon, y = rho, group = Sample_size)) + geom_line(aes(linetype=Sample_size, color=Sample_size))+xlab(expression(epsilon))+ylab('Absolute Error')+ scale_y_log10()+scale_x_log10()+
  theme(legend.key.height= unit(1, 'cm'), legend.key.width= unit(2, 'cm'),legend.title = element_text(size = 18),legend.text = element_text( size = 16))+ guides(linetype = guide_legend(nrow = 1))
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
load("~/Documents/research/CCA/Simulation/1_22/KM/K0.05.RData")
epsilon<-rep(1:5,times=4)
Sample_size<-rep(c('n=50','n=100','n=200','n=500'),each=5)
rho<-as.vector(t(E[1:4,]))
df_rho<-data.frame(epsilon=epsilon,Sample_size=Sample_size,rho=rho )
df_rho_new<-df_rho
df_rho_new$Sample_size<- factor(df_rho_new$Sample_size, levels = c('n=50','n=100','n=200','n=500'))
FPC_rho<-ggplot(data = df_rho_new, mapping = aes(x = epsilon, y = rho, group = Sample_size)) + geom_line(aes(linetype=Sample_size, color=Sample_size))+xlab('k')+ylab(expression(paste("Absolute Error of ",rho)))+ scale_y_log10()+theme(legend.position="top")+ guides(linetype = guide_legend(nrow = 1))
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
grid.arrange(legend, R_rho , R_F , R_G ,FPC_rho,FPC_F, FPC_G,
             ncol=2, nrow = 4, layout_matrix = cbind(c(5,6,7,1),c(2,3,4,1)),
             widths = c(3, 3), heights = c( 3,3,3,1))
##########
grid.arrange(legend, R_rho , R_F  ,FPC_rho,FPC_F, 
             ncol=2, nrow = 3, layout_matrix = cbind(c(4,5,1),c(2,3,1)),
             widths = c(3, 3), heights = c( 3,3,1))

