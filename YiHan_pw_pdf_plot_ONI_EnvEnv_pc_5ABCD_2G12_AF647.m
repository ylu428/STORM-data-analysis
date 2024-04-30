clc;clear;
addpath(genpath('M:\Yen-Cheng\YCC_Matlab'));

[Myfiles1,Mypath1] = uigetfile({'*.mat';'*.*'}, 'Select "_pdist_af647_af647.mat" files', 'MultiSelect','on');
load(strcat(Mypath1, Myfiles1));
dataA = cumdist_AF647;
clear cumdist_AF647

[Myfiles2,Mypath2] = uigetfile({'*.mat';'*.*'}, 'Select "_pdist_af647_af647.mat" files', 'MultiSelect','on');
load(strcat(Mypath2, Myfiles2));
dataB = cumdist_AF647;
clear cumdist_AF647

% load('5C_HXB2_TM1_SMLs-DB15nm60SMLs_pdist_af647_af647.mat');
% dataC = cumdist_AF647;
% clear cumdist_AF647
% 
% load('5D_HXB2_TM3_SMLs-DB15nm60SMLs_pdist_af647_af647.mat');
% dataD = cumdist_AF647;
% clear cumdist_AF647

%% normalized
low_end=0;
interval=1;
high_end=100;
edges=low_end:interval:high_end;

binCounts1  =  histcounts (dataA , edges);
binCounts2  =  histcounts (dataB , edges);
% binCounts3  =  histcounts (dataC , edges);
% binCounts4  =  histcounts (dataD , edges);

sampleHistc1 = binCounts1./sum(binCounts1);
sampleHistc2 = binCounts2./sum(binCounts2);
% sampleHistc3 = binCounts3./sum(binCounts3);
% sampleHistc4 = binCounts4./sum(binCounts4);

Counts1  =  binCounts1./sum(binCounts1);
Counts2  =  binCounts2./sum(binCounts2);
% Counts3  =  binCounts3./sum(binCounts3);
% Counts4  =  binCounts4./sum(binCounts4);

samplePDF1  =  Counts1(1:end);
samplePDF2  =  Counts2(1:end);
% samplePDF3  =  Counts3(1:end);
% samplePDF4  =  Counts4(1:end);

ref = rand([1 100000])*(high_end);
binCounts_ref  =  histcounts (ref , edges);
Counts_ref  =  binCounts_ref./sum(binCounts_ref);
samplePDF_ref  =  Counts_ref(1:end);

midpoint=(low_end+interval/2):interval:(high_end-interval/2);
Amax_coord=midpoint(find(samplePDF1==max(samplePDF1)));
Bmax_coord=midpoint(find(samplePDF2==max(samplePDF2)));
% Cmax_coord=midpoint(find(samplePDF3==max(samplePDF3)));
% Emax_coord=midpoint(find(samplePDF4==max(samplePDF4)));
max_coord=[Amax_coord,Bmax_coord];

figure('Name','Histogram, photon thr>1000')
plot(midpoint,samplePDF1,"-k",midpoint,samplePDF2,"-g",...
     'LineWidth',3);

% plot(midpoint,samplePDF1,"-k",midpoint,samplePDF2,"-g",...
%      midpoint,samplePDF3,"-b",midpoint,samplePDF4,"-r",...
%      'LineWidth',3);
% refline
ax = gca;
ax.FontSize = 18;
title('Env-Env (2G12-AF647)','FontSize',20)
xlabel('Pairwise distance (nm)','FontSize',20)
ylabel('Empirical PDF','FontSize',20)
legend({'Control','+SQV'},'FontSize',14,'Location','best','Orientation','Vertical')
% legend({'Control','+SQV','+IFITM1','+IFITM3'},'FontSize',14,'Location','best','Orientation','Vertical')

%% comparison


raw_dataA=quantile(dataA,[0.25 0.5 0.75]);
raw_dataB=quantile(dataB,[0.25 0.5 0.75]);
% raw_dataC=quantile(dataC,[0.25 0.5 0.75]);
% raw_dataD=quantile(dataD,[0.25 0.5 0.75]);
raw_dataD=[raw_dataA;raw_dataB];

%threshold
thr=[0,100];
dataAthr=dataA(dataA>thr(1) & dataA<thr(2));
dataBthr=dataB(dataB>thr(1) & dataB<thr(2));
% dataCthr=dataC(dataC>thr(1) & dataC<thr(2));
% dataDthr=dataD(dataD>thr(1) & dataD<thr(2));

Q_A=quantile(dataAthr,[0.25 0.75]);
Q_B=quantile(dataBthr,[0.25 0.75]);
% Q_C=quantile(dataCthr,[0.25 0.75]);
% Q_D=quantile(dataDthr,[0.25 0.75]);

IQR_A=Q_A(2)-Q_A(1);
IQR_B=Q_B(2)-Q_B(1);
% IQR_C=Q_C(2)-Q_C(1);
% IQR_D=Q_D(2)-Q_D(1);

W_Iz_A=2*IQR_A*size(dataAthr,1)^(-1/3);
W_Iz_B=2*IQR_B*size(dataBthr,1)^(-1/3);
% W_Iz_C=2*IQR_C*size(dataCthr,1)^(-1/3);
% W_Iz_D=2*IQR_D*size(dataDthr,1)^(-1/3);

[h_opt,p_AB,df_AB] = KStest_optbin_N100(dataAthr, dataBthr, W_Iz_A, W_Iz_B, 0.05);
% [h_opt,p_CD,df_CD] = KStest_optbin_N100(dataCthr, dataDthr, W_Iz_C, W_Iz_D, 0.05);
% [h_opt,p_AC,df_AC] = KStest_optbin_N100(dataAthr, dataCthr, W_Iz_A, W_Iz_C, 0.05);
% [h_opt,p_AD,df_AD] = KStest_optbin_N100(dataAthr, dataDthr, W_Iz_A, W_Iz_D, 0.05);
% [h_opt,p_BC,df_BC] = KStest_optbin_N100(dataBthr, dataCthr, W_Iz_B, W_Iz_C, 0.05);
% [h_opt,p_BD,df_BD] = KStest_optbin_N100(dataBthr, dataDthr, W_Iz_B, W_Iz_D, 0.05);
p_sum=[p_AB];
% p_sum=[p_AB;p_CD;p_AC;p_AD;p_BC;p_BD];
disp(strcat('p value = ', num2str(p_sum)))