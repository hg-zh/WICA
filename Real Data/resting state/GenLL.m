load 'E:\matlab\Data_list\ID_LIST.mat'
Data=cell([length(ID_LIST),1]);
for i=1:length(ID_LIST)
    id=ID_LIST(i);
    DIR=['E:\matlab\rfMRI\Data_extrac\',num2str(id),'.mat'];
	load(DIR);
    Data{i}=struct('ID',id,'Caudate_L',Caudate_L,'Caudate_R',Caudate_R,'Putamen_L',Putamen_L,'Putamen_R',Putamen_R);
    clearvars -except ID_LIST Data
end
%%
%ALL
%%%%density function
CAU_d=cell([length(ID_LIST),1]);PUT_d=cell([length(ID_LIST),1]);
for i=1:length(ID_LIST)
    CAU_tmp=[Data{i}.Caudate_L];
    PUT_tmp=[Data{i}.Putamen_L];
    CAUd_tmp=zeros(284,100);PUTd_tmp=zeros(284,100);
    for j=1:284
        CAUd_tmp(j,1:100)=histcounts(CAU_tmp(CAU_tmp(:,j)>9000&(CAU_tmp(:,j)<=16000),j) ,100);
        PUTd_tmp(j,1:100)=histcounts(PUT_tmp(PUT_tmp(:,j)>8000&(PUT_tmp(:,j)<=15000),j) ,100);
        CAUd_tmp(j,CAUd_tmp(j,:)==0)=sqrt(eps);PUTd_tmp(j,PUTd_tmp(j,:)==0)=sqrt(eps);
        CAUd_tmp(j,:)=100*CAUd_tmp(j,:)/sum(CAUd_tmp(j,:));PUTd_tmp(j,:)=100*PUTd_tmp(j,:)/sum(PUTd_tmp(j,:));        
    end
    CAU_d{i}=CAUd_tmp;PUT_d{i}=PUTd_tmp;
    clearvars CAUd_tmp PUTd_tmp CAU_tmp PUT_tmp
end
figure(1)
surf(CAU_d{1})
figure(2)
surf(PUT_d{1})
clearvars Data
%%%%probability distrition
CAU_p=cell([length(ID_LIST),1]);PUT_p=cell([length(ID_LIST),1]);
for i=1:length(ID_LIST)
    CAUp_tmp=zeros(284,101);PUTp_tmp=zeros(284,101);
    for j=1:284       
        CAUp_tmp(j,2:101)=cumsum(CAU_d{i}(j,:))/100;
        PUTp_tmp(j,2:101)=cumsum(PUT_d{i}(j,:))/100;
    end
    CAU_p{i}=CAUp_tmp;PUT_p{i}=PUTp_tmp;
    clearvars CAUp_tmp PUTp_tmp
end
%%%%quantile function
CAU_q=cell([length(ID_LIST),1]);PUT_q=cell([length(ID_LIST),1]);
for i=1:length(ID_LIST)
    CAUq_tmp=zeros(284,101);PUTq_tmp=zeros(284,101);
    for j=1:284       
        for k=2:100
            fun = @(x) CAU(x,i,j,k);
            CAUq_tmp(j,k)=fzero(fun,0.5);
            fun = @(x) PUT(x,i,j,k);
            PUTq_tmp(j,k)=fzero(fun,0.5);
        end
        CAUq_tmp(j,101)=1;PUTq_tmp(j,101)=1;
    end
    CAU_q{i}=CAUq_tmp;PUT_q{i}=PUTq_tmp;
end
clearvars CAUq_tmp PUTq_tmp
CAU=cell([209,1]);PUT=cell([209,1]);
for i=1:209
    CAU{i}=struct('pdf',CAU_d{i},'cdf',CAU_p{i},'quantile',CAU_q{i});
    PUT{i}=struct('pdf',PUT_d{i},'cdf',PUT_p{i},'quantile',PUT_q{i});
end
clearvars -except ID_LIST CAU PUT
save('LL916815')