load 'E:\matlab\Data_list\ID_LIST.mat'
for i=1:length(ID_LIST)
    id=ID_LIST(i);
	Data_path=['D:\motor\',num2str(id),'_3T_tfMRI_MOTOR_preproc.zip'];
	unzip(Data_path,'temp');
	data=ft_read_cifti(['.\temp\',num2str(id),'\MNINonLinear\Results\tfMRI_MOTOR_LR\tfMRI_MOTOR_LR_Atlas_MSMAll.dtseries.nii']);
    Caudate_L_Lable=find(data.brainstructure==8);
    Caudate_R_Lable=find(data.brainstructure==9);
    Putamen_L_Lable=find(data.brainstructure==18);
    Putamen_R_Lable=find(data.brainstructure==19);
    Caudate_L=data.dtseries(Caudate_L_Lable,:);
    Caudate_R=data.dtseries(Caudate_R_Lable,:);
    Putamen_L=data.dtseries(Putamen_L_Lable,:);
    Putamen_R=data.dtseries(Putamen_R_Lable,:);
    Caudate_L=Caudate_L(all(~isnan(Caudate_L),2),:);
    Caudate_R=Caudate_R(all(~isnan(Caudate_R),2),:);
    Putamen_L=Putamen_L(all(~isnan(Putamen_L),2),:);
    Putamen_R=Putamen_R(all(~isnan(Putamen_R),2),:);
    save(num2str(id),'Caudate_L','Caudate_R','Putamen_L','Putamen_R','Caudate_L_Lable','Caudate_R_Lable','Putamen_L_Lable','Putamen_R_Lable');
    rmdir('temp','s')
    clearvars -except ID_LIST id
end
