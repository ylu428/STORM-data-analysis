addpath(genpath('M:\Yen-Cheng\YCC_Matlab'));
dataA=A; %
dataB=B;


Q_A=quantile(dataA,[0.25 0.75]);
Q_B=quantile(dataB,[0.25 0.75]);


IQR_A=Q_A(2)-Q_A(1);
IQR_B=Q_B(2)-Q_B(1);

W_Iz_A=2*IQR_A*size(dataA,1)^(-1/3);
W_Iz_B=2*IQR_B*size(dataB,1)^(-1/3);
%%%% W_Iz is histogram bin size . Ref:
%%%% https://www.fmrib.ox.ac.uk/datasets/techrep/tr00mj2/tr00mj2/node24.html 



[h_AB, p_AB, df_AB] = KStest_optbin_N100(dataA, dataB, W_Iz_A, W_Iz_B, 0.05);


[p_AB]

% x = {'A'};
% y = {dataA};
% swarmchart(dataA, '.'); % Swarmchart function is supported in Matlab 2020 or later edition
ax1=subplot(2,1,1);
histogram(dataA, 10)
title('data A')
ax2=subplot(2,1,2);
histogram(dataB, 20)
linkaxes([ax1 ax2], 'x')
title('data B')
