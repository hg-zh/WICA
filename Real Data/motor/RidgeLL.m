load '/Volumes/Hang/matlab/2_16Gen/Data_struct2_20LL6.mat'
%%
%%%%%Mean Estimation%%%%%
CAU_mean_d=zeros(284,100);PUT_mean_d=zeros(284,100);
CAU_mean_p=zeros(284,101);PUT_mean_p=zeros(284,101);
CAU_mean_q=zeros(284,101);PUT_mean_q=zeros(284,101);
for i=1:209
    CAU_mean_q=CAU_mean_q+CAU{i}.quantile/209;
    PUT_mean_q=PUT_mean_q+PUT{i}.quantile/209;
    clearvars i
end
for j=1:284
    for k=1:100
        fun = @(x) CAU_F(x,j,k);
        CAU_mean_p(j,k)=fzero(fun,0.5);
        fun = @(x) PUT_F(x,j,k);
        PUT_mean_p(j,k)=fzero(fun,0.5);
        clearvars k
    end
    CAU_mean_p(j,101)=1;PUT_mean_p(j,101)=1;
    clearvars j fun
end
for j=1:284
    for k=1:100
        CAU_mean_d(j,k)=(CAU_mean_p(j,k+1)-CAU_mean_p(j,k))*100;
        PUT_mean_d(j,k)=(PUT_mean_p(j,k+1)-PUT_mean_p(j,k))*100;
        clearvars k
    end
    clearvars j
end
CAU_hat=struct('pdf',CAU_mean_d,'cdf',CAU_mean_p,'quantile',CAU_mean_q);
PUT_hat=struct('pdf',PUT_mean_d,'cdf',PUT_mean_p,'quantile',PUT_mean_q);
clearvars  CAU_mean_d CAU_mean_p CAU_mean_q PUT_mean_d PUT_mean_p PUT_mean_q
%%
%%%%%Log mapping%%%%%
CAU_T=cell([209,1]);PUT_T=cell([209,1]);
for i=1:209
    CAU_tmp=zeros(284,101);PUT_tmp=zeros(284,101);
    for j=1:284
        CAU_tmp(j,:)=interp1([0:0.01:1],CAU{i}.quantile(j,:) ,CAU_hat.cdf(j,:))-[0:0.01:1]; %#ok<*NBRAK>
        PUT_tmp(j,:)=interp1([0:0.01:1],PUT{i}.quantile(j,:) ,PUT_hat.cdf(j,:))-[0:0.01:1]; %#ok<*NBRAK>
    end
    CAU_T{i}=CAU_tmp;PUT_T{i}=PUT_tmp;
    clearvars  CAU_tmp PUT_tmp i j
end
CAU_hat.tangent=CAU_T;PUT_hat.tangent=PUT_T;
clearvars CAU_T PUT_T
%%
%%%%%Projection%%%%%
Basis_CAU=cell([50,1]);Basis_PUT=cell([50,1]);
for k=1:50
    tmp_C=zeros(284,101);tmp_P=zeros(284,101);
    for j=1:284
        tmp_C(j,:)=sqrt(2)*cos(k*pi*CAU_hat.cdf(j,:));
        tmp_P(j,:)=sqrt(2)*cos(k*pi*PUT_hat.cdf(j,:));
        clearvars j
    end
    Basis_CAU{k}=tmp_C;Basis_PUT{k}=tmp_P;
    clearvars  tmp_C tmp_P k
end
CAU_hat.basis=Basis_CAU;PUT_hat.basis=Basis_PUT;
clearvars Basis_CAU Basis_PUT
Xi_C=zeros(209,50);Xi_P=zeros(209,50);
for i=1:209
    for k=1:50
        tmp1=0;tmp2=0;
        for j=1:284
            tmp1=tmp1+sum(CAU_hat.tangent{i}(j,:).*CAU_hat.basis{k}(j,:).*[CAU_hat.pdf(j,:),0]);
            tmp2=tmp2+sum(PUT_hat.tangent{i}(j,:).*PUT_hat.basis{k}(j,:).*[PUT_hat.pdf(j,:),0]);
        end
        Xi_C(i,k)=tmp1*0.01/284; Xi_P(i,k)=tmp2*0.01/284;
    end
    clearvars  tmp1 tmp2 i k j
end
CAU_hat.xi=Xi_C;PUT_hat.xi=Xi_P;
clearvars Xi_C Xi_P
%%
[Train,Test_CAU,Test_PUT]=Split(5);
options = optimset('fminbnd');
options.TolX=1e-8;
epsilon=fminbnd(@CV,10^-8,10^-5,options);
% e=linspace(10^-8,10^-5,100);E=zeros(100,1);
% for i=1:100
%     E(i)=CV(e(i));
% end
% find(E==min(E))
% e=linspace(e(34),e(36),100);E=zeros(100,1);
% for i=1:100
%     E(i)=CV(e(i));
% end
% find(E==min(E))
% epsilon=e(45);
Ax=zeros(50,50);Ay=zeros(50,50);
Axy=zeros(50,50);Ayx=zeros(50,50);
for i=1:209
    for j=1:50
        for k=1:50
            Ax(j,k)=Ax(j,k)+CAU_hat.xi(i,j)*CAU_hat.xi(i,k)/209;
            Ay(j,k)=Ay(j,k)+PUT_hat.xi(i,j)*PUT_hat.xi(i,k)/209;
            Ayx(j,k)=Ayx(j,k)+PUT_hat.xi(i,j)*CAU_hat.xi(i,k)/209;
            Axy(j,k)=Axy(j,k)+CAU_hat.xi(i,j)*PUT_hat.xi(i,k)/209;
        end
    end
end
C=inv(Ax+epsilon*eye(50))*Axy*inv(Ay+epsilon*eye(50))*Ayx;
[V,D]=eig(C);D=sqrt(D);
R=real(V(:,1));
G=inv(Ay+epsilon*eye(50))*Ayx*R/norm(inv(Ay+epsilon*eye(50))*Ayx*R);
F_hat=zeros(284,101);G_hat=zeros(284,101);
for k=1:50
    F_hat=F_hat+R(k)*CAU_hat.basis{k};
    G_hat=G_hat+G(k)*PUT_hat.basis{k};
end
surf(F_hat)
heatmap(G_hat)
run('picdemo.m')
%%
m=KCV(15);
Ax_m=Ax(1:m,1:m);Ay_m=Ay(1:m,1:m);
Axy_m=Axy(1:m,1:m);Ayx_m=Ayx(1:m,1:m);
C_m=inv(Ax_m)*Axy_m*inv(Ay_m)*Ayx_m;
[V_m,D_m]=eig(C_m);D_m=sqrt(D_m);
R=real(V_m(:,1));
G=inv(Ay_m)*Ayx_m*R/norm(inv(Ay_m)*Ayx_m*R);
F_hat=zeros(284,101);G_hat=zeros(284,101);
for k=1:m
    F_hat=F_hat+R(k)*CAU_hat.basis{k};
    G_hat=G_hat+G(k)*PUT_hat.basis{k};
end
heatmap(G_hat)
run('picdemo.m')
ylabel('time')


