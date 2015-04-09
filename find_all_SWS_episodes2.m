function [t_mdpt_SWS,data_at_SWS_midpoints,t_mdpt_indices]=find_all_SWS_episodes2(datafile,epoch_length)
%USAGE:
% [t_mdpt_SWS,t_midpt_indices,data_at_SWS_midpoints]=find_all_SWS_episodes2(datafile)
% 
% This function finds all the episodes of length 5 minutes
% in which SWA makes up at least 90% of the activity in those 5 minutes
% and computes the median delta power (or lactate) in each of these SWS episodes
% and the time of episode midpoint. Similar but different from Franken et al. 2001
% Figure 1.  
%
% INPUTS: 
% datafile:  This is a data file containing 2 columns, where
%            sleep state is in the first column (0=wake,1=SWS,2=REM), lactate
%            or delta power is in the second colum
% epoch_length: the length of the scoring epoch (in seconds)
%

%
% OUTPUTS:
% t_mdpt_SWS:  the times (in hours) of the midpoints of the 5 minute episodes
% in which SWA makes up at least 90% of data. 
%
% data_at_SWS_midpoints:  median delta power or lactate 
%                         signal at the midpoints of SWS episodes
%                         longer than 5 minutes
%
% t_mdpt_indices:  the indices of t (original time vector)
%                  corresponding to the midpoints of the SWS episodes.


% First set up a vector of time in hours
t_hours=0:1/(60*60/epoch_length):(1/(60*60/epoch_length))*(size(datafile,1)-1);  % convert seconds to hours
                                        


data=datafile(:,2);
 

% initialize
t_mdpt_SWS=0;
data_at_SWS_midpoints=0;
t_mdpt_indices=0;
counter=0;  % counter for number of SWS episodes longer than 5 min.
starting_indices=0;

% for each 5 min sliding window check to see if 90% or more is SWA
SWSwindow_length = 5;  % minutes
percentage = .9;    % .9 means 90% of the window of length SWSwindow_length needs to be SWS to be counted in the analysis
rows_in_SWS_episode = SWSwindow_length*60/epoch_length;

for i=rows_in_SWS_episode:size(datafile,1)
  if length(find(datafile(i-(rows_in_SWS_episode-1):i,1)==1))>=percentage*rows_in_SWS_episode
    counter = counter+1;
    starting_indices(counter)=i-(rows_in_SWS_episode-1);
    % t_mdpt_SWS(counter) = mean([t_hours(i-15),t_hours(i-14)]);
    % data_at_SWS_midpoints(counter) = median(data(i-29:i));
    % t_mdpt_indices(counter)=i-15;
 
  end
end

%disp(['5 min windows with 90 percent SWA: ' num2str(counter)])

% now combine overlapping windows
first_index_of_streak=1;
streak_counter=1;
while first_index_of_streak < length(starting_indices)
  c = max(find(starting_indices<(rows_in_SWS_episode+starting_indices(first_index_of_streak))));

%end_of_streak=starting_indices(c)+30
 
  t_mdpt_SWS(streak_counter) = mean([t_hours(starting_indices(first_index_of_streak)),t_hours(starting_indices(c)+(rows_in_SWS_episode-1))]);
  data_at_SWS_midpoints(streak_counter) = median(data(starting_indices(first_index_of_streak):starting_indices(c)+(rows_in_SWS_episode-1)));
  t_mdpt_indices(streak_counter)=starting_indices(first_index_of_streak)+(rows_in_SWS_episode/2-1);
  first_index_of_streak=c+1;  
  streak_counter=streak_counter+1;  

end

if numel(t_mdpt_indices)==1 & t_mdpt_indices==0
  error('I couldn''t find any 5 minute episodes containing at least 90% SWS')
end


%figure
%plot(t_hours,data,t_mdpt_SWS,data_at_SWS_midpoints,'o')
%plot(t_mdpt_SWS,data_at_SWS_midpoints,'o')
%title('this is just midpoints of SWS episodes, not a model fit')
%xlabel('Time (hours)')
%hold off

% improve this plot to include the sleep state data 
