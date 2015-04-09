function [t_mdpt_SWS,data_at_SWS_midpoints,t_mdpt_indices]=find_all_SWS_episodes(datafile,signal,epoch_length)
%USAGE:
% [t_mdpt_SWS,t_midpt_indices,data_at_SWS_midpoints]=find_all_SWS_episodes(data,signal,epoch_length)
% 
% This function finds all the SWS episodes of length > 5 minutes
% and computes the median delta power (or lactate) in each of these SWS episodes
% and the time of episode midpoint. Following Franken et al. 2001
% Figure 1.  
%
% INPUTS: 
% datafile:    This is a data file containing 4 columns, where
%          sleep state is in the first column (0=wake,1=SWS,2=REM), lactate
%          is in the second column, average delta power in EEG1 is in the 
%          third column, and averge delta power in EEG2 is in the 4th column.
%
% signal: 'lactate' or 'delta1' or 'delta2'
%
% OUTPUTS:
% t_mdpt_SWS:  the times of the midpoints of the SWS episodes
% longer than 5 minutes. 
%
% t_midpt_SWS:  the times corresponding to the midpoints of the SWS
%               episodes
% 
% data_at_SWS_midpoints:  median delta power or lactate 
%                         signal at the midpoints of SWS episodes
%                         longer than 5 minutes
%
% t_mdpt_indices:  the indices of t (orignal time vector)
%                  corresponding to the midpoints of the SWS episodes.


% First set up a vector of time in hours
t_hours=0:1/(60*60/epoch_length):(1/(60*60/epoch_length))*(size(datafile,1)-1);  %converting seconds to hours
                                        


% Pick off the correct data from datafile
if strcmp(signal,'delta1')
  data=datafile(:,3);
elseif strcmp(signal,'delta2')
  data=datafile(:,4);
elseif strcmp(signal,'lactate')
  data=datafile(:,2);
end

% Use contiguous.m to find all SWS episodes longer than 5 minutes
runs=contiguous(datafile(:,1),1); %all contiguous "runs", even runs of only 1
allruns=runs{1,2};  % go from cell array to matrix



% initialize
t_mdpt_SWS=0;
data_at_SWS_midpoints=0;
t_mdpt_indices=0;
counter=0;  % counter for number of SWS episodes longer than 5 min.

% find all SWS episodes > 5 min
SWS_episode_length = 5;  % in minutes
rows_in_SWS_episode = SWS_episode_length*60/epoch_length;
for i=1:size(allruns,1)
  if (allruns(i,2)-allruns(i,1)) > rows_in_SWS_episode-1  % 30 rows means 5 minutes (29 inclusive)
    counter = counter+1;

    t_mdpt_SWS(counter) = t_hours(floor((allruns(i,1)+ ...
                                              allruns(i,2))/2));
    data_at_SWS_midpoints(counter) = median(data(allruns(i,1): ...
                                                 allruns(i,2)));
    t_mdpt_indices(counter)=floor((allruns(i,1)+allruns(i,2))/2);
  

  end
end





figure
only_sleep_indices=find(datafile(:,1)==1);
sleep_delta=datafile(only_sleep_indices,4);
scatter(t_hours(only_sleep_indices),sleep_delta,25,'r')
hold on
%plot(t_hours,data,t_mdpt_SWS,data_at_SWS_midpoints,'o')
plot(t_mdpt_SWS,data_at_SWS_midpoints,'go')
title('this is just midpoints of SWS episodes, not a model fit')
xlabel('Time (hours)')
hold off

% improve this plot to include the sleep state data 
