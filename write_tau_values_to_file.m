function write_tau_values_to_file(txtfiles,directory,model,signal,algorithm,do_rescore,do_restrict,restrict,Taui,TauD)
% USAGE: write_tau_values_to_file(txtfiles,directory,model,signal,algorithm,do_rescore,do_restrict,restrict,Taui,TauD)
%
%
%
addpath 'C:\Users\wisorlab\Documents\MATLAB\Brennecke\matlab-pipeline\Matlab\etc\matlab-utils\'



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


sheet_agree = xl.addSheets({'Tau Values'});
sheet_params = xl.addSheets({'Calling Parameters'});

xl.setCells(sheet_agree{1},[1,2],txtfiles','false','true');
xl.setCells(sheet_agree{1},[2,2],Taui');
xl.setCells(sheet_agree{1},[3,2],TauD');
xl.setCells(sheet_agree{1},[1,1],{'Files'},'669999');
xl.setCells(sheet_agree{1},[2,1],{'Ti'},'669999','true');
xl.setCells(sheet_agree{1},[3,1],{'Td'},'669999','true');

xl.setCells(sheet_params{1},[1,2],{'Model'});
xl.setCells(sheet_params{1},[1,3],{'Signal'});
xl.setCells(sheet_params{1},[1,4],{'Algorithm'});
xl.setCells(sheet_params{1},[1,5],{'Rescored into Quiet Wake and Active Wake?'});
xl.setCells(sheet_params{1},[1,6],{'Was the dataset restricted?'});
xl.setCells(sheet_params{1},[1,7],{'Start Time in hours from beginning of recording:'});
xl.setCells(sheet_params{1},[1,8],{'End Time in hours from the beginning of the recording:'});
xl.setCells(sheet_params{1},[2,2],{model});
xl.setCells(sheet_params{1},[2,3],{signal});
xl.setCells(sheet_params{1},[2,4],{algorithm});
xl.setCells(sheet_params{1},[2,5],do_rescore);
xl.setCells(sheet_params{1},[2,6],do_restrict);

if do_restrict ==1
	xl.setCells(sheet_params{1},[2,7],restrict(1));
	xl.setCells(sheet_params{1},[2,8],restrict(2));
else
	xl.setCells(sheet_params{1},[2,7],{'N/A'});
	xl.setCells(sheet_params{1},[2,8],{'N/A'});
end



xl.sourceInfo(mfilename('fullpath'));
xl.rmDefaultSheets();

xl.saveAs('tau_values.xls',strcat(directory,most_recent_directory,'\'));
fclose('all');  %so Excel doesn't think MATLAB st