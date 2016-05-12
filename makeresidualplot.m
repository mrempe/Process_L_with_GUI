% maketauplot.m
%
% This function makes a bar plot of the residuals computed from running the ProcessL on the strain data using lactate 
% as the signal (run using strain_data_master_script.m)
% 
% This script requires the Matlab workspace strain_study_lactate_tau_residuals_starting_at17_00.mat to be loaded


residuals3_means = [mean(residual_AKR3) mean(residual_BA3) mean(residual_BL3) mean(residual_CDJ3) mean(residual_DBA3) mean(residual_OLA3)];
residuals5_means = [mean(residual_AKR5) mean(residual_BA5) mean(residual_BL5) mean(residual_CDJ5) mean(residual_DBA5) mean(residual_OLA5)];

residuals3_errors = [std(residual_AKR3)/sqrt(length(residual_AKR3)) std(residual_BA3)/sqrt(length(residual_BA3)) ...
 					 std(residual_BL3)/sqrt(length(residual_BL3)) 	std(residual_CDJ3)/sqrt(length(residual_CDJ3)) ...
 					 std(residual_DBA3)/sqrt(length(residual_DBA3)) std(residual_OLA3)/sqrt(length(residual_OLA3))];
residuals5_errors = [std(residual_AKR5)/sqrt(length(residual_AKR5)) std(residual_BA5)/sqrt(length(residual_BA5)) ...
 					 std(residual_BL5)/sqrt(length(residual_BL5)) 	std(residual_CDJ5)/sqrt(length(residual_CDJ5)) ...
 					 std(residual_DBA5)/sqrt(length(residual_DBA5)) std(residual_OLA5)/sqrt(length(residual_OLA5))];					 


figure
b=bar(1:6,[residuals3_means' residuals5_means'],1);
b(1).EdgeColor = 'none';
b(2).EdgeColor = 'none';
b(1).FaceColor = [0.5 0.5 0.5];
b(2).FaceColor = [0 0 0];

axis([0 7 0 max([residuals3_means residuals5_means])+0.1*max([residuals3_means residuals5_means])])
set(gca,'box','off')
ax=gca;
ax.XTickLabel ={'AKR'; 'BA'; 'BL'; 'CDJ'; 'DBA'; 'OLA'};
ax.XTickLabelRotation = 45;


xlabel('Strain')
ylabel('NRMSE')
title('Normalized root-mean-square deviation by strain')
legend('3-state model','5-state model')
legend boxoff


hold on
% then add errorbars and asterisks
numgroups = 6;  
numbars = 2;    % number of bars per group

	x1=(1:numgroups)-groupwidth/2 +(2*1-1)*groupwidth/(2*numbars);
	x2=(1:numgroups)-groupwidth/2 +(2*2-1)*groupwidth/(2*numbars);
	errorbar(x1,residuals3_means',residuals3_errors','k','linestyle','none')
	errorbar(x2,residuals5_means',residuals5_errors','k','linestyle','none')

hold off


