function [out1,out2,out3]=Split(Fold)
    global CAU_mean_q
    global PUT_mean_q
    CAU= evalin('base', 'CAU');PUT= evalin('base', 'PUT');
    Train=cell([Fold,1]);TestT_CAU=cell([Fold,1]);TestT_PUT=cell([Fold,1]);
    for fold=1:Fold
        bin=200/Fold;
        Test_id=((fold-1)*bin+1):(fold*bin);
        Train_id=setdiff(1:209,Test_id);
        Train_CAU=cell([length(Train_id ),1]);
        Train_PUT=cell([length(Train_id ),1]);
        Test_CAU=cell([length(Test_id ),1]);Test_PUT=cell([length(Test_id ),1]);
        for i=1:length(Train_id)
            Train_CAU{i}=CAU{Train_id(i)};
            Train_PUT{i}=PUT{Train_id(i)};
        end
        for i=1:length(Test_id)
            Test_CAU{i}=CAU{Test_id(i)};
            Test_PUT{i}=PUT{Test_id(i)};
        end
       %%
        %Mean estimation
        CAU_mean_d=zeros(284,100);PUT_mean_d=zeros(284,100);
        CAU_mean_p=zeros(284,101);PUT_mean_p=zeros(284,101);
        CAU_mean_q=zeros(284,101);PUT_mean_q=zeros(284,101);
        for i=1:length(Train_id)
            CAU_mean_q=CAU_mean_q+Train_CAU{i}.quantile/length(Train_id);
            PUT_mean_q=PUT_mean_q+Train_PUT{i}.quantile/length(Train_id);
        end
        for j=1:284
            for k=1:100
                fun = @(x) CAU_FS(x,j,k);
                CAU_mean_p(j,k)=fzero(fun,0.5);
                fun = @(x) PUT_FS(x,j,k);
                PUT_mean_p(j,k)=fzero(fun,0.5);
            end
            CAU_mean_p(j,101)=1;PUT_mean_p(j,101)=1;
            clearvars fun
        end
        for j=1:284
            for k=1:100
                CAU_mean_d(j,k)=(CAU_mean_p(j,k+1)-CAU_mean_p(j,k))*100;
                PUT_mean_d(j,k)=(PUT_mean_p(j,k+1)-PUT_mean_p(j,k))*100;
            end
        end
        CAU_hat=struct('pdf',CAU_mean_d,'cdf',CAU_mean_p,'quantile',CAU_mean_q);
        PUT_hat=struct('pdf',PUT_mean_d,'cdf',PUT_mean_p,'quantile',PUT_mean_q);
        clearvars  CAU_mean_d CAU_mean_p CAU_mean_q PUT_mean_d PUT_mean_p PUT_mean_q
       %%
        %%%%%Log mapping%%%%%
        CAU_T=cell([length(Train_id),1]);PUT_T=cell([length(Train_id),1]);
        for i=1:length(Train_id)
            CAU_tmp=zeros(284,101);PUT_tmp=zeros(284,101);
            for j=1:284
                CAU_tmp(j,:)=interp1([0:0.01:1],Train_CAU{i}.quantile(j,:) ,CAU_hat.cdf(j,:))-[0:0.01:1]; %#ok<*NBRAK>
                PUT_tmp(j,:)=interp1([0:0.01:1],Train_PUT{i}.quantile(j,:) ,PUT_hat.cdf(j,:))-[0:0.01:1]; %#ok<*NBRAK>
            end
            CAU_T{i}=CAU_tmp;PUT_T{i}=PUT_tmp;
            clearvars  CAU_tmp PUT_tmp i j
        end
        CAU_hat.tangent=CAU_T;PUT_hat.tangent=PUT_T;
        clearvars CAU_T PUT_T
        CAU_T=cell([length(Test_id),1]);PUT_T=cell([length(Test_id),1]);
        for i=1:length(Test_id)
            CAU_tmp=zeros(284,101);PUT_tmp=zeros(284,101);
            for j=1:284
                CAU_tmp(j,:)=interp1([0:0.01:1],Test_CAU{i}.quantile(j,:) ,CAU_hat.cdf(j,:))-[0:0.01:1]; %#ok<*NBRAK>
                PUT_tmp(j,:)=interp1([0:0.01:1],Test_PUT{i}.quantile(j,:) ,PUT_hat.cdf(j,:))-[0:0.01:1]; %#ok<*NBRAK>
            end
            CAU_T{i}=CAU_tmp;PUT_T{i}=PUT_tmp;
            clearvars  CAU_tmp PUT_tmp i j
        end
        TestT_CAU{fold}=CAU_T;TestT_PUT{fold}=PUT_T;
        clearvars CAU_T PUT_T
       %%
        %%%%%Projection%%%%%
        Basis_CAU=cell([50,1]);Basis_PUT=cell([50,1]);
        for k=1:50
               %Basis_CAU{k}=ones(284,101);Basis_PUT{k}=ones(284,101);
                tmp_C=zeros(284,101);tmp_P=zeros(284,101);
                for j=1:284
                    tmp_C(j,:)=sqrt(2)*cos((k)*pi*CAU_hat.cdf(j,:));
                    tmp_P(j,:)=sqrt(2)*cos((k)*pi*PUT_hat.cdf(j,:));
                    clearvars j
                end
                Basis_CAU{k}=tmp_C;Basis_PUT{k}=tmp_P;
                clearvars  tmp_C tmp_P k
        end
        CAU_hat.basis=Basis_CAU;PUT_hat.basis=Basis_PUT;
        clearvars Basis_CAU Basis_PUT
        Xi_C=zeros(length(Train_id),50);Xi_P=zeros(length(Train_id),50);
        for i=1:length(Train_id)
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
        Ax=zeros(50,50);Ay=zeros(50,50);
        Axy=zeros(50,50);Ayx=zeros(50,50);
        for i=1:length(Train_id)
            for j=1:50
                for k=1:50
                    Ax(j,k)=Ax(j,k)+CAU_hat.xi(i,j)*CAU_hat.xi(i,k)/length(Train_id);
                    Ay(j,k)=Ay(j,k)+PUT_hat.xi(i,j)*PUT_hat.xi(i,k)/length(Train_id);
                    Ayx(j,k)=Ayx(j,k)+PUT_hat.xi(i,j)*CAU_hat.xi(i,k)/length(Train_id);
                    Axy(j,k)=Axy(j,k)+CAU_hat.xi(i,j)*PUT_hat.xi(i,k)/length(Train_id);
                end
            end
        end
        Train{fold}=struct('Ax',Ax,'Ay',Ay,'Axy',Axy,'Ayx',Ayx,'CAU',CAU_hat,'PUT',PUT_hat);
    end
    out1=Train;out2=TestT_CAU;out3=TestT_PUT;
end