% maketauplot.m
%
% This function makes a bar plot of the residuals computed from running the ProcessL on the strain data using lactate 
% as the signal (run using strain_data_master_script.m)
% 
% This script requires the Matlab workspace strain_study_lactate_tau_residuals_starting_at17_00.mat to be loaded

% ***
%   Taui
% ***


Taui3_means = [mean(Taui_AKR3) mean(Taui_BA3) mean(Taui_BL3) mean(Taui_CDJ3) mean(Taui_DBA3) mean(Taui_OLA3)];
Taui5_means = [mean(Taui_AKR5) mean(Taui_BA5) mean(Taui_BL5) mean(Taui_CDJ5) mean(Taui_DBA5) mean(Taui_OLA5)];

Taui3_errors = [std(Taui_AKR3)/sqrt(length(Taui_AKR3)) std(Taui_BA3)/sqrt(length(Taui_BA3)) ...
 					 std(Taui_BL3)/sqrt(length(Taui_BL3)) 	std(Taui_CDJ3)/sqrt(length(Taui_CDJ3)) ...
 					 std(Taui_DBA3)/sqrt(length(Taui_DBA3)) std(Taui_OLA3)/sqrt(length(Taui_OLA3))];
Taui5_errors = [std(Taui_AKR5)/sqrt(length(Taui_AKR5)) std(Taui_BA5)/sqrt(length(Taui_BA5)) ...
 					 std(Taui_BL5)/sqrt(length(Taui_BL5)) 	std(Taui_CDJ5)/sqrt(length(Taui_CDJ5)) ...
 					 std(Taui_DBA5)/sqrt(length(Taui_DBA5)) std(Taui_OLA5)/sqrt(length(Taui_OLA5))];					 


figure
subplot(1,2,1)
b=bar(1:6,[Taui3_means' Taui5_means'],1);
b(1).EdgeColor = 'none';
b(2).EdgeColor = 'none';
b(1).FaceColor = [0.5 0.5 0.5];
b(2).FaceColor = [0 0 0];

axis([0 7 0 max([Taui3_means Taui5_means])+0.1*max([Taui3_means Taui5_means])])
set(gca,'box','off')
ax=gca;
ax.XTickLabel ={'AKR'; 'BA'; 'BL'; 'CDJ'; 'DBA'; 'OLA'};
ax.XTickLabelRotation = 45;


xlabel('Strain')
ylabel('Tau_i')
title('\tau_i by strain')
legend('3-state model','5-state model')
legend boxoff


hold on
% then add errorbars and asterisks
numgroups = 6;  
numbars = 2;    % number of bars per group

	x1=(1:numgroups)-groupwidth/2 +(2*1-1)*groupwidth/(2*numbars);
	x2=(1:numgroups)-groupwidth/2 +(2*2-1)*groupwidth/(2*numbars);
	errorbar(x1,Taui3_means',Taui3_errors','k','linestyle','none')
	errorbar(x2,Taui5_means',Taui5_errors','k','linestyle','none')

hold off


% ***
% Taud
% ***



TauD3_means = [mean(TauD_AKR3) mean(TauD_BA3) mean(TauD_BL3) mean(TauD_CDJ3) mean(TauD_DBA3) mean(TauD_OLA3)];
TauD5_means = [mean(TauD_AKR5) mean(TauD_BA5) mean(TauD_BL5) mean(TauD_CDJ5) mean(TauD_DBA5) mean(TauD_OLA5)];

TauD3_errors = [std(TauD_AKR3)/sqrt(length(TauD_AKR3)) std(TauD_BA3)/sqrt(length(TauD_BA3)) ...
 					 std(TauD_BL3)/sqrt(length(TauD_BL3)) 	std(TauD_CDJ3)/sqrt(length(TauD_CDJ3)) ...
 					 std(TauD_DBA3)/sqrt(length(TauD_DBA3)) std(TauD_OLA3)/sqrt(length(TauD_OLA3))];
TauD5_errors = [std(TauD_AKR5)/sqrt(length(TauD_AKR5)) std(TauD_BA5)/sqrt(length(TauD_BA5)) ...
 					 std(TauD_BL5)/sqrt(length(TauD_BL5)) 	std(TauD_CDJ5)/sqrt(length(TauD_CDJ5)) ...
 					 std(TauD_DBA5)/sqrt(length(TauD_DBA5)) std(TauD_OLA5)/sqrt(length(TauD_OLA5))];					 


subplot(1,2,2)
b=bar(1:6,[TauD3_means' TauD5_means'],1);
b(1).EdgeColor = 'none';
b(2).EdgeColor = 'none';
b(1).FaceColor = [0.5 0.5 0.5];
b(2).FaceColor = [0 0 0];

axis([0 7 0 max([TauD3_means TauD5_means])+0.1*max([TauD3_means TauD5_means])])
set(gca,'box','off')
ax=gca;
ax.XTickLabel ={'AKR'; 'BA'; 'BL'; 'CDJ'; 'DBA'; 'OLA'};
ax.XTickLabelRotation = 45;


xlabel('Strain')
ylabel('Tau_D')
title('\tau_D by strain')
legend('3-state model','5-state model')
legend boxoff


hold on
% then add errorbars and asterisks
numgroups = 6;  
numbars = 2;    % number of bars per group

	x1=(1:numgroups)-groupwidth/2 +(2*1-1)*groupwidth/(2*numbars);
	x2=(1:numgroups)-groupwidth/2 +(2*2-1)*groupwidth/(2*numbars);
	errorbar(x1,TauD3_means',TauD3_errors','k','linestyle','none')
	errorbar(x2,TauD5_means',TauD5_errors','k','linestyle','none')

hold off


