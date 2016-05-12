%plotEMGforWakeandREM.m
%
% this script reads in a .txt file and makes a plot of the EMG values for every epoch marked
% as wake or REMS, ordered by EMG value from largest to smallest


% Select the .txt file
% Pop up a window 
[files,directory] = uigetfile('Multiselect','on','\\FS1\WisorData\Rempe\Data\strain_study_data\*.txt','PROCESS_L   Please Select .txt file(s)');  %last parameter sent to uigetfile ('*.edf*) specifies that only edf files will be displayed in the user interface.
if ~iscell(files), files = {files}; end