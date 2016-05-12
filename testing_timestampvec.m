 % testing some of the new find_all_SWS_episodes stuff using timestampvec

 num_epochs_in_window     = (SWSwindow_length*60)/epoch_length_in_seconds;
 num_windows_in_recording = floor(length(data)/num_epochs_in_window);

d=diff(timestampvec);
 % loop through num_windows_in_recording and determine if there is a jump in the timestamps or not, indicating missing data
 for i=1:num_windows_in_recording
 	if max(d( (i-1)*num_epochs_in_window+1:i*num_epochs_in_window-1)) == seconds(epoch_length_in_seconds)
			does_not_contain_jump(i)=1;
		else
			does_not_contain_jump(i)=0;
	end 
end

non_missing_data_windows = find(does_not_contain_jump);

SWS_episode_index=1
for i=non_missing_data_windows
	if length(find(sleepstate((i-1)*num_epochs_in_window+1:i*num_epochs_in_window)==1))>=percentage*num_epochs_in_window
		t_mdpt_SWS(SWS_episode_index) = timestampvec((i-1)*num_epochs_in_window+1)+seconds(epoch_length_in_seconds*num_epochs_in_window/2);
  		data_at_SWS_midpoints(SWS_episode_index) = median(data((i-1)*num_epochs_in_window+1:i*num_epochs_in_window));
  		t_mdpt_indices(SWS_episode_index)=round(mean([(i-1)*num_epochs_in_window+1,i*num_epochs_in_window]));  % these are indices after removing artifacts
  		SWS_episode_index = SWS_episode_index + 1;
  	end
end