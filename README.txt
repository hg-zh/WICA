For simulation:
1. The folder 'Case1' and 'Case2' contain all the codes and results of Case 1 and Case 2 in the simulation section.
2. The folder 'Case1' contains 4 folders, 'KM' and 'RidgeM' are used to derive Figure 2, while 'KCV' and 'RidgeCV' are used for Tables 1-3 and Figure 3. The folder 'Case2' contains 2 folders, 'KCV' and 'RidgeCV' are used for Tables 1-3. All these aforementioned folders contains a .R file and five .RData files that are obtained by setting different tuning parameters in the corresponding code. 
3. The folder 'Plot' contains the code to obtain Figure 2 and 3. In addition, before ploting the result you need to save the data obtained by the code in 'Case1' and 'Case2' and change the path in the .R file.

For real data analysis:
1. The data are downloaded from https://db.humanconnectome.org/app/template/Index.vm.
2. We choose the subjects with complete 3T motor task fMRI in the S1200 quarter released, and 'ID_LIST.mat' contains the id numbers of subjects we used.
3. Run ‘extrac.m’ and ‘GenLL’ to pre-process the data.
4. ’RidgeLL’ is the main file, run this file to get the results in Section 5.