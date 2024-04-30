addpath(genpath('M:\Yen-Cheng\YCC_Matlab'));
dataA=A; %
dataB=B;


Q_A=quantile(dataA,[0.25 0.75]);
Q_B=quantile(dataB,[0.25 0.75]);


IQR_A=Q_A(2)-Q_A(1);
IQR_B=Q_B(2)-Q_B(1);

W_Iz_A=2*IQR_A*size(dataA,1)^(-1/3);
W_Iz_B=2*IQR_B*size(dataB,1)^(-1/3);
%%%% W_Iz = 2(IQR)*N^(-1/3) is the histogram bin size summerized by Izenman, A. J. (1991). 
%%%% Ref: https://www.fmrib.ox.ac.uk/datasets/techrep/tr00mj2/tr00mj2/node24.html 



[h_AB, p_AB, df_AB] = KStest_optbin_N100(dataA, dataB, W_Iz_A, W_Iz_B, 0.05);
% h = kstest(x) returns a test decision for the null hypothesis that the  
% data in vector x comes from a standard normal distribution, against the 
% alternative that it does not come from such a distribution, using the 
% one-sample Kolmogorov-Smirnov test. The result h is 1 if the test 
% rejects the null hypothesis at the 5% significance level, or 0 otherwise.

[p_AB]


%%%% Plotting Swarmchart %%%%
x1 = ones(length(dataA), 1);
x2 = ones(length(dataB), 1)*2;


swarmchart(x1, dataA, 'filled'); % Swarmchart function is supported in Matlab 2020 or later edition
hold on
swarmchart(x2, dataB, 'filled'); % Swarmchart function is supported in Matlab 2020 or later edition
hold on
boxplot([dataA; dataB], [x1; x2]) % With ";" is vertcat. Without ";" is horicat.
hold on

plot([x1(1) x2(1)], [1 1]*max(max(A), max(B))*1.2, '-k')
plot([x1(1) x1(1)], [max(A)*1.15  max(max(A), max(B))*1.2], '-k')
plot([x2(1) x2(1)], [max(B)*1.15  max(max(A), max(B))*1.2], '-k')
ylim([0  max(max(A), max(B))*1.35])
if p_AB > 0.05
    signif = 'ns';
elseif p_AB <= 0.05 && p_AB > 0.01
    signif = '*';
elseif p_AB <= 0.01 && p_AB > 0.001
    signif = '**';
elseif p_AB < 0.001 && p_AB > 0.0001
    signif = '***';
else
    signif = '****';
end

% Texts and labels
str = strcat(signif, ' (p = ', num2str(p_AB), ')');
dim = [0.45 0.8 0.15 0.1]; %[x_begin y_begin length height]
annotation("textbox", dim, 'String', str, FitBoxToText='on', LineStyle='none')

Labels = {'HXB2', 'HXB2-SQV'};
set(gca, 'XTick', 1:6, 'XTickLabel', Labels);
ylabel('SMLs per virus')
hold off
