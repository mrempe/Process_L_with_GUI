% script to plot the data comparing 5state model to 3state model using shift work
%
% All of the outputs from the script are in matrices where each row corresponds to a separate freq. range for delta. 
% Each column corresponds to a different animal

load ('FiveStateVsThreeStateVariableSWA.mat')



% Residuals of AW vs frequency band
figure
h1=plot(1:size(residual5stateAW,1),mean(residual5stateAW,2));
hold on 
h2=plot(1:size(residual5stateAW,1),mean(residual5stateAW,2),'.','MarkerSize',12);
h3=plot(1:size(residual3stateAW,1),mean(residual3stateAW,2));
h4=plot(1:size(residual3stateAW,1),mean(residual3stateAW,2),'o','MarkerSize',5);
errorbar(1:size(residual5stateAW,1),mean(residual5stateAW,2),std(residual5stateAW,'omitnan',2)/sqrt(size(residual5stateAW,1)));
errorbar(1:size(residual3stateAW,1),mean(residual3stateAW,2),std(residual3stateAW,'omitnan',2)/sqrt(size(residual3stateAW,1)));
hold off
set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend('5 state', '3 state')
xlabel('Frequency Range used for SWA')
ylabel('Residuals')
title('Active Phase Workers')
ax = gca;
ax.XTick = 1:length(residual5stateAW);
ax.XTickLabel = freq_ranges_vector;
axis([0 length(residual5stateAW)+1 0 max([mean(residual5stateAW,2);mean(residual3stateAW,2)])+0.1*max([mean(residual5stateAW,2);mean(residual3stateAW,2)])])


% Residuals of RW vs frequency band
figure
h1=plot(1:size(residual5stateRW,1),mean(residual5stateRW,2));
hold on 
h2=plot(1:size(residual5stateRW,1),mean(residual5stateRW,2),'.','MarkerSize',12);
h3=plot(1:size(residual3stateRW,1),mean(residual3stateRW,2));
h4=plot(1:size(residual3stateRW,1),mean(residual3stateRW,2),'o','MarkerSize',5);
errorbar(1:size(residual5stateRW,1),mean(residual5stateRW,2),std(residual5stateRW,'omitnan',2)/sqrt(size(residual5stateRW,1)));
errorbar(1:size(residual3stateRW,1),mean(residual3stateRW,2),std(residual3stateRW,'omitnan',2)/sqrt(size(residual3stateRW,1)));
hold off
set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend('5 state', '3 state')
xlabel('Frequency Range used for SWA')
ylabel('Residuals')
title('Resting Phase Workers')
ax = gca;
ax.XTick = 1:length(residual5stateRW);
ax.XTickLabel = freq_ranges_vector;
axis([0 length(residual5stateRW)+1 0 max([mean(residual5stateRW,2);mean(residual3stateRW,2)])+0.1*max([mean(residual5stateRW,2);mean(residual3stateRW,2)])])





