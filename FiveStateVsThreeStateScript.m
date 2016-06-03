% Script to run PROCESSLBATCHMODE_NOGUI using various frequency ranges 
% to define SWA (1-4, 2-4, 3-5, etc. )
clear 

% turn off new figure generation
set(0,'DefaultFigureVisible','off')


% Directories
AWdirectory = '\\FS1\WisorData\Gronli\Night work, animal model\txt files\Work period\Baseline and workdays\txt files\ActivePhaseWorkers\';
RWdirectory = '\\FS1\WisorData\Gronli\Night work, animal model\txt files\Work period\Baseline and workdays\txt files\RestingPhaseWorkers\';

% testing:
% AWdirectory = '\\FS1\WisorData\Gronli\Night work, animal model\txt files\Work period\Baseline and workdays\txt files\testing\';
% RWdirectory = '\\FS1\WisorData\Gronli\Night work, animal model\txt files\Work period\Baseline and workdays\txt files\testing\';

% vector of frequency ranges
freq_ranges_vector = ['1-4'; '1-5'; '1-6'; '1-7'; '1-8'; '2-4'; '2-5'; '2-6'; '2-7'; '2-8'; '3-4'; '3-5'; '3-6'; '3-7'; '3-8'];
%freq_ranges_vector = ['1-4';  '2-4';  '3-4' ]  % for testing
freq_ranges_vector = cellstr(freq_ranges_vector);  % convert the character array to a cell array where each element is a freq. range

% 5state model with AW
% --------------------
for i=1:length(freq_ranges_vector)
	[~,~,~,residual5stateAW(i,:),~,~,~,~,~,Taui5stateAW(i,:),TauD5stateAW(i,:)]=PROCESSLBATCHMODE_NOGUI(AWdirectory,'EEG1',freq_ranges_vector{i},'5state');
end
disp('-------------------------------------- 25% finished -------------------------------------------------------------')
 
% 5state model with RW
% --------------------
for i=1:length(freq_ranges_vector)
	[~,~,~,residual5stateRW(i,:),~,~,~,~,~,Taui5stateRW(i,:),TauD5stateRW(i,:)]=PROCESSLBATCHMODE_NOGUI(RWdirectory,'EEG1',freq_ranges_vector{i},'5state');
end
disp('------------------------------------- Halfway done -------------------------------------------------------------')
 
% 3state model with AW
% --------------------
for i=1:length(freq_ranges_vector)
	[~,~,~,residual3stateAW(i,:),~,~,~,~,~,Taui3stateAW(i,:),TauD3stateAW(i,:)]=PROCESSLBATCHMODE_NOGUI(AWdirectory,'EEG1',freq_ranges_vector{i},'3state');
end
disp('------------------------------------- 75% finished -------------------------------------------------------------')


% 3state model with RW
% --------------------
for i=1:length(freq_ranges_vector)
	[~,~,~,residual3stateRW(i,:),~,~,~,~,~,Taui3stateRW(i,:),TauD3stateRW(i,:)]=PROCESSLBATCHMODE_NOGUI(RWdirectory,'EEG1',freq_ranges_vector{i},'3state');
end


% turn figures back on
set(0,'DefaultFigureVisible','on')

% Save everything to a workspace
disp('----------------- Saving data to workspace -------------------------------------------------------------------------')
save('FiveStateVsThreeStateVariableSWA.mat','freq_ranges_vector','residual*','Tau*')

