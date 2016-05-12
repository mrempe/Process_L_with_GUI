% script to run ProcessLBatchmode_NOGUI on the shift work data during recovery
% using a moving 24-hour window


% directory where the data files are

AWLD_directory = 'C:\Users\wisorlab\Desktop\testAWLD\';




% % make cell arrays of the different filenams, based on protocol
% % AW LD  Active phase workers in Light/Dark conditions
% AWLD_files = dir(strcat(directory,'*AW_LD*.txt'));
% RWLD_files = dir(strcat(directory,'*RW_LD*.txt'));
% AWDD_files = dir(strcat(directory,'*AW_DD*.txt'));
% RWDD_files = dir(strcat(directory,'*RW_DD*.txt'));

% AWLD_files = struct2cell(AWLD_files);
% AWLD_files(2:end,:)=[];   %remove all other fields like date, bytes, isdir, datenum, etc.  
% if ~iscell(AWLD_files), AWLD_files = {AWLD_files}; end

% RWLD_files = struct2cell(RWLD_files);
% RWLD_files(2:end,:)=[];   %remove all other fields like date, bytes, isdir, datenum, etc.  
% if ~iscell(RWLD_files), RWLD_files = {RWLD_files}; end

% AWDD_files = struct2cell(AWDD_files);
% AWDD_files(2:end,:)=[];   %remove all other fields like date, bytes, isdir, datenum, etc.  
% if ~iscell(AWDD_files), AWDD_files = {AWDD_files}; end

% RWDD_files = struct2cell(RWDD_files);
% RWDD_files(2:end,:)=[];   %remove all other fields like date, bytes, isdir, datenum, etc.  
% if ~iscell(RWDD_files), RWDD_files = {RWDD_files}; end

signal = 'EEG1';
model = '5state';

% calculate total number of windows I will need.
% total windows needed = total number of hours - window width + 1
window_width = 24;  % hours  
total_windows = 144 - window_width + 1;



% AWLD
for i=1:total_windows
	start_index = 8650 + 360*(i-1);
	end_index   = start_index + window_width*60*6;  % 8640 for a 24 hour window
	[signal_data,state_data,timestampvec,residual,best_S,UppA,LowA,dynamic_range,Timer,TauiAWLD{i},TauDAWLD{i}]=PROCESSLBATCHMODE_NOGUI(AWLD_directory,signal,model,start_index,end_index); 
end 



% Get averages (across recordings) of tau values for each moving window
mean_across_time_TauiAWLD = mean([TauiAWLD{:}],1);


figure
plot(mean_across_time_TauiAWLD,'.','MarkerSize',10)





% Save all the data into a matlab workspace
% save moving_window_in_recovery_data.mat