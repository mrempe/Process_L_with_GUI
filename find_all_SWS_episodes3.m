function [t_mdpt_SWS,data_at_SWS_midpoints,t_mdpt_indices]=find_all_SWS_episodes3(datafile,timestampvec,epoch_length)
%USAGE:
% [t_mdpt_SWS,t_midpt_indices,data_at_SWS_midpoints]=find_all_SWS_episodes2(datafile)
% 
% This function finds all the episodes of length 5 minutes
% in which SWA makes up at least 90% of the activity in those 5 minutes
% and computes the median delta power (or lactate) in each of these SWS episodes
% This version uses the matlab data structure datetime. 
% NOTE: the 5-minute episodes are not overlapping, so if half of one 
% 5-minute episode is SWS and half of the next one is, that 5-minute episode 
% of SWS will not be detected.  find_all_SWS_episodes4.m uses overlapping sliding
% 5-minute windows, but the resulting data points are much noisier.  
% and the time of episode midpoint. Similar but different from Franken et al. 2001
% Figure 1.  
%
% INPUTS: 
% datafile:  This is a data file containing 2 columns, where
%            sleep state is in the first column (0=wake,1=SWS,2=REM), lactate
%            or delta power is in the second column
% timestampvec: a datetime vector of the timestamps of each epoch that is included after throwing out artifacts
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
%t_hours=0:1/(60*60/epoch_length):(1/(60*60/epoch_length))*(size(datafile,1)-1);  % convert seconds to hours
                                        

sleepstate=datafile(:,1);
data=datafile(:,2);
 

% initialize
%t_mdpt_SWS=0;  % don't initialize this.  MATLAB will think it is a double, not a datetime object
data_at_SWS_midpoints=0;
%t_mdpt_indices=0;  % don't initialize this.  MATLAB will think it is a double, not a datetime object.
starting_indices=0;

% for each 5 min sliding window check to see if 90% or more is SWA
SWSwindow_length = 5;  % minutes
percentage = .9;    % .9 means 90% of the window of length SWSwindow_length needs to be SWS to be counted in the analysis
%rows_in_SWS_episode = SWSwindow_length*60/epoch_length;

num_epochs_in_window     = (SWSwindow_length*60)/epoch_length;
num_windows_in_recording = floor(length(data)/num_epochs_in_window);

d=diff(timestampvec);
 % loop through num_windows_in_recording and determine if there is a jump in the timestamps or not, indicating missing data
 for i=1:num_windows_in_recording
  if max(d( (i-1)*num_epochs_in_window+1:i*num_epochs_in_window-1)) == seconds(epoch_length)
      does_not_contain_jump(i)=1;
    else
      does_not_contain_jump(i)=0;
  end 
end

non_missing_data_windows = find(does_not_contain_jump);

SWS_episode_index=1;
for i=non_missing_data_windows
  if length(find(sleepstate((i-1)*num_epochs_in_window+1:i*num_epochs_in_window)==1))>=percentage*num_epochs_in_window
    t_mdpt_SWS(SWS_episode_index) = timestampvec((i-1)*num_epochs_in_window+1)+seconds(epoch_length*num_epochs_in_window/2);
    data_at_SWS_midpoints(SWS_episode_index) = median(data((i-1)*num_epochs_in_window+1:i*num_epochs_in_window));
    t_mdpt_indices(SWS_episode_index)=round(mean([(i-1)*num_epochs_in_window+1,i*num_epochs_in_window]));  % these are indices after removing artifacts
    SWS_episode_index = SWS_episode_index + 1;
    end
end

if (numel(t_mdpt_indices)==1 & t_mdpt_indices==0) 
  error('I couldn''t find any 5 minute episodes containing at least 90% SWS')
end


%figure
%plot(t_hours,data,t_mdpt_SWS,data_at_SWS_midpoints,'o')
%plot(t_mdpt_SWS,data_at_SWS_midpoints,'o')
%title('this is just midpoints of SWS episodes, not a model fit')
%xlabel('Time (hours)')
%hold off

% improve this plot to include the sleep state data 
