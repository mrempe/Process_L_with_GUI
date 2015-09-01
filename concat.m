% ======================================================
% Concatenation script for TXT files
% 
% Concatenates two TXT files into an Excel workbook
% and interpolates data in the time lapsed between the files
% 
% created by @jonbrennecke
% ======================================================
clear

% instantiate an ActiveX connection with Excel
xl = XL;

% ***type this in before you start*** addpath Brennecke\matlab-pipeline\Matlab\etc\matlab-utils\

% open file dialog
[files,path] = uigetfile({'*txt','Text Files (*.txt)';'*.*','All Files';},'Select Two Files to Concatenate','MultiSelect','On');
if ischar(files), files = {files}; end

% file name
t1 = regexp(files{1},'(\w*)','match');
t2 = regexp(files{2},'(\w*)','match');
% t1 = regexp(files{1},'(\w+\d+)\s+(\d+)_(\d+)_(\d+)','tokens');
% t2 = regexp(files{2},'(\w+\d+)\s+(\d+)_(\d+)_(\d+)','tokens');
sheets = xl.addSheets({[ t1{1} '-' t2{1} '_' t1{2} ]});

% loop through files
for i=1:length(files)
    
    file = fread(fopen([path files{i}],'r'),'*char');
    lines{i} = textscan(file,'%s','delimiter','\n');

    % after the first file, compute the difference between the files
    if i==2 
    	% calculate the time difference between the two files (in seconds)
    	% and interpolate times between the two

        tlast = regexp(lines{1,1}{1,1}{j},'(\d+/\d+/\d+)','match');
        tstamp1 = regexp(lines{1,1}{1,1}{j},'(\d+):(\d+):(\d+)','tokens');
        tstamp2 = regexp(lines{1,2}{1,1}{3},'(\d+):(\d+):(\d+)','tokens');
        tstampprevious = regexp(lines{1,1}{1,1}{j-1},'(\d+):(\d+):(\d+)','tokens');
        tsec = [str2num([tstamp1{1,1}{1};tstamp2{1,1}{1}])*60^2,str2num([tstamp1{1,1}{2};tstamp2{1,1}{2}])*60,str2num([tstamp1{1,1}{3};tstamp2{1,1}{3}])];
        tdelta = sum(tsec(2,:)-tsec(1,:)); % time gap that needs to be filled in (in seconds)

        % compute the epoch length 
        dt = str2num([tstampprevious{1,1}{3};tstamp1{1,1}{3}])
        if dt(2)==0 dt(2)=60; end % handle the case where the second time stamp is 0, indicating the turn of the minute
        epoch_length = dt(2)-dt(1) % epoch length (in seconds)
        

        if tdelta > 300
            disp(['WARNING: There is a gap of ' num2str(tdelta/60) ' minutes between the two files.'])        
        end

        % create duplicate rows bridging the gap between the files

        for k=1:floor(tdelta/epoch_length)
        	% convert tdelta (seconds) back into sec, min, hour
        	t = sum(tsec(1,:),2) + epoch_length*k;
        	s = mod(t,60); m = floor(mod(t-s,60^2)/60); h = floor((t-s-m)/60^2);


        	% write to Excel
        	xl.setCells(sheets{1},[1,size(lines{1,1}{1,1},1)+k],{[cell2mat(tlast) ',' num2str(h,'%02d') ':' num2str(m,'%02d') ':' num2str(s,'%02d')],line{1,1}{2:end}}, 'FFEE00');
        end
    end
    
%     outlines = cell(size(lines{1,i}{1,1},1),1);
%     
%     for j=1:length(lines{1,i}{1,1})
%         line = textscan(lines{1,i}{1,1}{j},'%s','delimiter','\t');
%         outlines{j,:} = line{1,1}';
%         
%        
% %         lines{1,i}{j,:} = line{1,1}';
%         lines{1,i}{1,1}{j,1} = line{1,1}';
%     end

    % copy the first and second files into Excel
    for j=1:length(lines{1,i}{1,1})
        line = textscan(lines{1,i}{1,1}{j},'%s','delimiter','\t');
        if i>1 & j>2
            xl.setCells(sheets{1},[1,size(lines{1,1}{1,1},1)+j+k-3],line{1,1}');
        else
            xl.setCells(sheets{1},[1,j],line{1,1}');
        end
    end
end
fclose('all')

