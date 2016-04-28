function write_tau_values_to_file(txtfiles,directory,model,signal, algorithm,do_rescore, do_restrict_hours_from_start, ...
								  do_restrict_start_time_end_time, do_restrict_using_epochs,restrict_hours_from_start, restrict_start_clock_time, ...
								  restrict_end_clock_time, restrict_epoch_start, restrict_epoch_end,Taui, TauD, LA, UA)
% USAGE: write_tau_values_to_file(txtfiles,directory,model,signal,algorithm,do_rescore,do_restrict,restrict,Taui,TauD,LA,UA)
%
%
%
addpath 'C:\Users\wisorlab\Documents\MATLAB\Brennecke\matlab-pipeline\Matlab\etc\matlab-utils\'

if length(LA{1})==1
	LAandUAareScalar = 1;
else
	LAandUAareScalar = 0;
end






% If Taui and Taud and not column vectors, make them column vectors
if size(Taui,2) ~= 1
	Taui=Taui';
end
if size(TauD,2) ~= 1
	TauD=TauD';
end


% first find the most recent directory located in "directory"
files = dir(directory);
directories = files([files.isdir]);
isgood = ones(length(directories),1);
for i=1:length(directories)
	if strncmpi(directories(i).name,'.',1)  % if directory name starts with a period, disregard it
		isgood(i) = 0;
	end
end

directories = directories(find(isgood));
dates = [directories.datenum];
[~,newestIndex] = max(dates);
most_recent_directory = directories(newestIndex).name;

% now write an excel spreadsheet to most_recent_directory (this is just the name of the directory, not the entire path)
%xl=XL(strcat(directory,most_recent_directory,'\','agreement_stats.xls')); % this didn't work because the .xls file didn't exist yet. 
xl=XL;
%sheet = xl.Sheets.Item(1);
%[numcols,numrows] = xl.sheetSize(sheet);
NA = 'N/A';

sheet_agree = xl.addSheets({'Tau Values and LA + UA'});
sheet_params = xl.addSheets({'Calling Parameters'});

xl.setCells(sheet_agree{1},[1,2],txtfiles','false','true');
xl.setCells(sheet_agree{1},[2,2],Taui);
xl.setCells(sheet_agree{1},[3,2],TauD);
if LAandUAareScalar
	xl.setCells(sheet_agree{1},[4,2],LA);
	xl.setCells(sheet_agree{1},[5,2],UA);
else 
	xl.setCells(sheet_agree{1},[4,2],repmat(NA,length(Taui),1));
	xl.setCells(sheet_agree{1},[5,2],repmat(NA,length(Taui),1));	
end 
xl.setCells(sheet_agree{1},[1,1],{'Files'},'669999');
xl.setCells(sheet_agree{1},[2,1],{'Ti'},'669999','true');
xl.setCells(sheet_agree{1},[3,1],{'Td'},'669999','true');
xl.setCells(sheet_agree{1},[4,1],{'LA'},'669999','true');
xl.setCells(sheet_agree{1},[5,1],{'UA'},'669999','true');


xl.setCells(sheet_params{1},[1,2],{'Model'});
xl.setCells(sheet_params{1},[1,3],{'Signal'});
xl.setCells(sheet_params{1},[1,4],{'Algorithm'});
xl.setCells(sheet_params{1},[1,5],{'Rescored into Quiet Wake and Active Wake?'});
xl.setCells(sheet_params{1},[1,6],{'Was the dataset restricted in hours from start of recording?'});
xl.setCells(sheet_params{1},[1,7],{'Did the user input a start time and end time to restrict the recording?'})
xl.setCells(sheet_params{1},[1,8],{'Did the user input a starting epoch and ending epoch to restrict the recording?'})
xl.setCells(sheet_params{1},[1,9],{'Start Time in hours from beginning of recording:'});
xl.setCells(sheet_params{1},[1,10],{'End Time in hours from the beginning of the recording:'});
xl.setCells(sheet_params{1},[1,11],{'Start Time (in clock time):'});
xl.setCells(sheet_params{1},[1,12],{'End Time (in clock time):'});
xl.setCells(sheet_params{1},[1,13],{'Starting Epoch (in epochs from beginning of recording):'});
xl.setCells(sheet_params{1},[1,14],{'Ending Epoch (in epochs from beginning of recording):'});
xl.setCells(sheet_params{1},[2,2],{model});
xl.setCells(sheet_params{1},[2,3],{signal});
xl.setCells(sheet_params{1},[2,4],{algorithm});
xl.setCells(sheet_params{1},[2,5],do_rescore);
xl.setCells(sheet_params{1},[2,6],do_restrict_hours_from_start);
xl.setCells(sheet_params{1},[2,7],do_restrict_start_time_end_time);
xl.setCells(sheet_params{1},[2,8],do_restrict_using_epochs);


if do_restrict_hours_from_start
	xl.setCells(sheet_params{1},[2,9], restrict_hours_from_start(1));
	xl.setCells(sheet_params{1},[2,10],restrict_hours_from_start(2));
	xl.setCells(sheet_params{1},[2,11],{'N/A'});
	xl.setCells(sheet_params{1},[2,12],{'N/A'});
	xl.setCells(sheet_params{1},[2,13],{'N/A'}); 
	xl.setCells(sheet_params{1},[2,14],{'N/A'}); 

elseif do_restrict_start_time_end_time
	xl.setCells(sheet_params{1},[2,9],{'N/A'});
	xl.setCells(sheet_params{1},[2,10],{'N/A'});
	xl.setCells(sheet_params{1},[2,11],restrict_start_clock_time);
	xl.setCells(sheet_params{1},[2,12],restrict_end_clock_time);
	xl.setCells(sheet_params{1},[2,13],{'N/A'});
	xl.setCells(sheet_params{1},[2,14],{'N/A'}); 

elseif do_restrict_using_epochs
	xl.setCells(sheet_params{1},[2,9],{'N/A'});
	xl.setCells(sheet_params{1},[2,10],{'N/A'});
	xl.setCells(sheet_params{1},[2,11],{'N/A'});
	xl.setCells(sheet_params{1},[2,12],{'N/A'});
	xl.setCells(sheet_params{1},[2,13],restrict_epoch_start);
	xl.setCells(sheet_params{1},[2,14],restrict_epoch_end);
else
	xl.setCells(sheet_params{1},[2,9],{'N/A'});
	xl.setCells(sheet_params{1},[2,10],{'N/A'});
	xl.setCells(sheet_params{1},[2,11],{'N/A'});
	xl.setCells(sheet_params{1},[2,12],{'N/A'});
	xl.setCells(sheet_params{1},[2,13],{'N/A'});
	xl.setCells(sheet_params{1},[2,14],{'N/A'});
end



xl.sourceInfo(mfilename('fullpath'));
xl.rmDefaultSheets();

xl.saveAs('tau_values.xls',strcat(directory,most_recent_directory,'\'));
fclose('all');  %so Excel doesn't think MATLAB st