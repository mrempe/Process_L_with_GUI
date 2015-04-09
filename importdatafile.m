function [data,textdata]=importdatafile(FileToRead,directory)
% usage: [data,textdata]=importdatafile(FileToRead,directory)

% this function uses textscan to read in a .txt data file
% and convert it to two matrices, one for numeric data and 
% one for text
% textscan is more robust than importdata.m as textscan 
% will not quit like importdata does if there is missing data
%
% Also, I'm adding some code to check for negative lactate values 
% at the end of the file, and to cut off the file before the 
% negative values if they exist.  


%OUTPUTS:
%      data     This is a matrix containing 2 fewer rows than
%               the original .txt file with only the numeric
%               data.  

%      textdata This is a cell array containing the strings in 
%               columns 1 and 2 of the data file (timestamp and sleep state)
%               NOTE1: use curly brackets to access elements in 
%               textdata. i.e. textdata{3,2} for the sleep state 
%               on line 3.  
%               NOTE2: textdata has already had the first two
%               rows removed.  It is just strings on lines 
%               containing EEG data. 

if nargin ==1 
	directory='';
end


DELIMITER = sprintf('\t','');     % tab delimited
HEADERLINES = 2-1;               % there are 2 header lines, but one gets read in the call to fgetl before textscan is called



filename=strcat(directory,FileToRead);

fid=fopen(filename);                          % compute number of columns present
tLines = fgetl(fid);
NUM_COLUMNS = numel(strfind(tLines,DELIMITER)) + 1;    


format=['%s%s' repmat('%f', [1 NUM_COLUMNS-2])];    % minus 2 for the 2 columns of text
%format=['%s',repmat('%f',1,NUM_COLUMNS)];
c=textscan(fid,format,'Headerlines',HEADERLINES,'Delimiter',DELIMITER,'CollectOutput',1,'ReturnOnError',0);
textdata=c{1};
data=c{2};
fclose(fid);



% check for negative values in the lactate signal (first
% column of data) at the beginning or the end of 
% the file.  If so, cut off that data

% first find all locations where lactate < 0 and
% remove them from data and textdata
locations=find(data(:,1)<0);
data(locations,:)=[];
textdata(locations,:)=[];

