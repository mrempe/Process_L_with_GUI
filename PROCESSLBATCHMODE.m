function [signal_data,state_data,residual,best_S,UppA,LowA,dynamic_range,Timer,Taui,TauD]=PROCESSLBATCHMODE
%
% PREVIOUS USAGE: [signal_data,state_data,residual,best_S,UppA,LowA,dynamic_range,Timer,Taui,Taud] =PROCESSLBATCHMODE(directory,signal,algorithm,keyword,restrict)
%
% INPUTS:
% directory is where the data files are (include the entire directory. i.e. 'D:\mrempe\strain_study_data\BL\long_files\'
%
% signal can be one of three options: 'lactate', 'delta1' or 'delta2',
% where delta1 means the delta power in EEG1 and delta2 is delta power in EEG2
%
% algorithm is either 'NelderMead' or 'BruteForce' and determines how the optimum 
% values of taui and taud are found. 
% 
% keyword: optional argument to specify a keyword in each .txt file you want to process. 
% For example, if instead of processing all files that end in .txt if you want to process 
% all files that have the work AUTO in the filename, let keyword be 'AUTO'
%
% restrict: optional argument to restrict the simulation to only use data during a baseline period,
%           a sleep deprivation period, or a recovery period.  Valid options are a two-element vector containing 
%           the start time and end time of the data you want to keep, measured in hours from the beginning of the 
%           recording. Other options are  'baseline','SD','recovery', or 'none'.
%           The 'none' option means to use the entire recording.
%  
% OUTPUTS:
% signal_data: a cell array containing vectors of either delta power or lactate, one for each file in the directory
% state_data: a cell array containing vectors of the sleep state (0,1,2,3,4) for each file in the directory
% residual:   a vector containing a vector of mean sum of squares error between best fit and data for each file in the directory
% best_S:  a cell array containing the best fit curve S, one for each file in the directory
% UppA:    upper asymptote
% LowA:    lower asymptote
% Taui:    a vector of the rise time time constant (during wake, REMS, or active wake), one value for each file in the directory
% Taud:    a vector of the fall time time constant, (during SWS and quiet wake) one value for each file in the directory
% 

%profile -memory on



addpath 'C:\Users\wisorlab\Documents\MATLAB\Brennecke\matlab-pipeline\Matlab\etc\matlab-utils\';  %where importdatafile.m XL.m and create_TimeStampMatrix_from_textdata.m live


% Pop up a window 
[files,directory] = uigetfile('Multiselect','on','D:\*.txt','Please Select .txt file(s)');  %last parameter sent to uigetfile ('*.edf*) specifies that only edf files will be displayed in the user interface.
if ~iscell(files), files = {files}; end



prompt = {'Do you want to use EEG1, EEG2, or lactate?', ...
'Do you want to use BruteForce or NelderMead?','Do you want to rescore the data into Quiet Wake and Active Wake and use a 5-state model rather than a 3-state model?  (1 for yes, 0 for no)', ...
'Do you want to restrict the data? (1 for yes, 0 for no)','Would you like to write the optimal tau values to a file with a Data Source Info Tab? (1 for yes, 0 for no)'};
defaults = {'EEG2','NelderMead','0','0','0'}; 
dlg_title = 'Input';
inputs = inputdlg(prompt,dlg_title,1,defaults,'on');


signal=inputs{1};
algorithm = inputs{2};
do_rescore = str2double(inputs{3});
do_restrict = str2double(inputs{4});
do_write_tau_values_to_file = str2double(inputs{5});

if do_rescore == 1
  model = '5state';
elseif do_rescore == 0
  model = '3state';
else
  error('You entered something other than 0 or 1 for rescoring into QW and AW')
end

% % Handle the case where a keyword is given
% if use_keyword == 1
%   prompt2 = {'Which Keyword Shall I look for in the filenames of the files you wish to process?'};
%   defaults2 = {'AUTO'};
%   dlg_title2 = 'Keyword';
%   keyword_input = inputdlg(prompt2,dlg_title2,1,defaults2,'on');
%   keyword = keyword_input{1};
% else 
%   keyword = 'None';
% end

% Handle the case where a restriction is given
if do_restrict ==1
  prompt3 = {'Please enter a start time in hours from the beginning of the recording', 'Ending time (in hours from beginning of recording)'};
  defaults3 = {'14','20'};
  dlg_title3 = 'How do you want to restrict the data? ';
  restrict_input = inputdlg(prompt3,dlg_title3,1,defaults3,'on');
  restrict = [str2double(restrict_input{1}) str2double(restrict_input{2})];
else
   restrict = 'none';
end 

% Set up directory_plus_extension, based on whether a keyword was entered
% if keyword ~= 'None'
%     directory_plus_extension = strcat(directory,'*',keyword,'.txt');
%   else 
directory_plus_extension = strcat(directory,'*.txt');
%end

% if nargin==3
% directory_plus_extension=strcat(directory,'*.txt');
% restrict = 'none';
% elseif nargin==4 
%   if  keyword ~= 'none' & keyword ~= 'None'
%     directory_plus_extension = strcat(directory,'*',keyword,'.txt');
%   else 
%     directory_plus_extension = strcat(directory,'*.txt');
%   end
%   restrict = 'none';

% elseif nargin==5 
%   if keyword ~= 'none' & keyword ~= 'None'
%     directory_plus_extension=strcat(directory,'*',keyword,'.txt');
%   else
%   directory_plus_extension=strcat(directory,'*.txt');
%   end 
% end


if ~strcmp(signal,'lactate') & ~strcmp(signal,'delta1') & ~strcmp(signal,'delta2') & ~strcmp(signal,'EEG1') & ~strcmp(signal,'EEG2')
  error('Input signal must be one of the following: ''lactate'', ''delta1'', ''delta2'', ''EEG1'', or ''EEG2''')
end




%files = dir(directory_plus_extension);     % the dir function returns a cell array containing the name, date, bytes and date-as-number as a single array for each txt file in this directory.

% %%%%%%%%%%%%%%%%%%%%%%%%%
% for i=length(files):-1:1                % don't include files that have already been autoscored
%   fname = files(i).name;
%   if strfind(fname,'AUTOSCORED')
%     files(i)=[];
%   end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%


HowManyFiles = length(files) % Need to know the number of files to process.  This number is encoded as the variable "HowManyFiles". 


% ---
% LOADING LOOP
% First loop through the files, load all the data and decide if 
% I will use each data set or not
% --- 
for FileCounter=1:length(files)  %this loop imports the data files one-by-one and processes the data in them into output files.   
clear PhysioVars EEG1 EEG2 EMG EMG_data
clear EEG1
%clear dynamic_range
%clear TimeStampMatrix   

  %[data,textdata]=importdatafile(files(FileCounter).name,directory);%importfile returns data (a matrix) and textdata (a cell array)
  [data,textdata]=importdatafile(files{FileCounter},directory);%importfile returns data (a matrix) and textdata (a cell array)

  display(files{FileCounter}) % One matrix (textdata) holds the date/time stamp and sleep state.  The other (data) holds the lactate and EEG data.
  

% Compute the length of one epoch here using TimeStampMatrix
% I used to just compute the difference in seconds between first time stamp and second,
% but that was not robust enough since sometimes there is a change in minutes between those two timestamps
% textdata has mm/dd/yyy,HH:MM:SS AM format, but may or may not include quotes or single quotes,
% so to be more robust, find the second colon and use that to get the seconds field of the time stamp
f=find(textdata{1,1}==':');   % Find all locations of the colon in the first time stamp
first_colon_loc = f(1);   
last_colon_loc = f(2);
hour_first_time_stamp    = str2num(textdata{1,1}(first_colon_loc-2))*10+str2num(textdata{1,1}(first_colon_loc-1));  
hour_second_time_stamp   = str2num(textdata{2,1}(first_colon_loc-2))*10+str2num(textdata{2,1}(first_colon_loc-1));  
minute_first_time_stamp  = str2num(textdata{1,1}(first_colon_loc+1))*10+str2num(textdata{1,1}(first_colon_loc+2));
minute_second_time_stamp = str2num(textdata{2,1}(first_colon_loc+1))*10+str2num(textdata{2,1}(first_colon_loc+2));
second_first_time_stamp  = str2num(textdata{1,1}(last_colon_loc+1))*10+str2num(textdata{1,1}(last_colon_loc+2));
second_second_time_stamp = str2num(textdata{2,1}(last_colon_loc+1))*10+str2num(textdata{2,1}(last_colon_loc+2));

epoch_length_in_seconds(FileCounter)=etime([2014 2 28 hour_second_time_stamp minute_second_time_stamp second_second_time_stamp],[2014 2 28 hour_first_time_stamp minute_first_time_stamp second_first_time_stamp]);

window_length = 4;      % length of moving window used to compute UA and LA if signal is lactate.  
  if strcmp(signal,'lactate')      % cut off data if using lactate sensor, smooth lactate data
    lactate_cutoff_time_hours=60;  % time in hours to cut off the lactate signal (lifetime of sensor)
    lactate_cutoff_time_rows=lactate_cutoff_time_hours*60*60/epoch_length_in_seconds(FileCounter);
    
    if epoch_length_in_seconds >=10
    LactateSmoothed=medianfiltervectorized(data(:,1),1);
    %[numchanged,LactateSmoothed]=SmootheLactate(data(:,1));
    %disp(['number of points smoothed:', num2str(numchanged)])
    % figure
    % plot(1:size(data(:,1),1),data(:,1),'r',1:size(data(:,1),1),LactateSmoothed,'b')
  elseif epoch_length_in_seconds < 10
    LactateSmoothed=medianfiltervectorized(data(:,1),2);
  end



   if size(data,1) > lactate_cutoff_time_rows
      data=data(1:lactate_cutoff_time_rows,:);
    end
  end

                                                       
 
  PhysioVars=zeros(size(data,1),4);
 

  missing_values=0;
  for i = 1: size(data,1)  
    
    if isempty(textdata{i,2})==1        % call unscored epochs wake
      missing_values=missing_values+1;
      PhysioVars(i,1)=0;  

    elseif textdata{i,2}=='W' 
      PhysioVars(i,1)=0;
    elseif textdata{i,2}=='S'
      PhysioVars(i,1)=1;
    elseif textdata{i,2}=='P'
      PhysioVars(i,1)=2;
    elseif textdata{i,2}=='R'
      PhysioVars(i,1)=2;
    elseif sum(textdata{i,2}=='Tr')==2
      PhysioVars(i,1)=0;                 % call transitions wake
    elseif textdata{i,2}=='X'            %artefact
      PhysioVars(i,1)=5; 
        else   
          error('I found a sleep state that wasn''t W,S,P,R,Tr, or X');
    end
  end
  disp(['There were ',num2str(missing_values), ' epochs that were not scored.'])

 

  if strcmp(signal,'lactate') 
    PhysioVars(:,2)=LactateSmoothed(1:size(data,1));
  disp(['Average lactate power: ' num2str(mean(PhysioVars(:,2)))])
  else PhysioVars(:,2)=data(:,1);
  end
  
  % Find the columns with EEG1 1-2 Hz and EEG2 1-2 Hz  
  fid = fopen(strcat(directory,files{FileCounter}));
  tLine1 = fgetl(fid);
  tLine2 = fgetl(fid);
  ColumnsHeads = textscan(tLine1,'%s','delimiter', sprintf('\t'));
  HeadChars=char(ColumnsHeads{1,1});
  for i=1:length(HeadChars)
    EEG1(i)=~isempty(strfind(HeadChars(i,:),'EEG 1'));
    EEG2(i)=~isempty(strfind(HeadChars(i,:),'EEG 2'));
    onetotwo(i)=~isempty(strfind(HeadChars(i,:),'1-2 '));
    EMG(i)=~isempty(strfind(HeadChars(i,:),'EMG'));
  end
  fclose('all');   %close all the files


  EEG1_1to2Hzcolumn = intersect(find(EEG1),find(onetotwo))-2; % subtract 2 to account for the fact that the first two columns are timestamp and lactate
  EEG2_1to2Hzcolumn = intersect(find(EEG2),find(onetotwo))-2;
  EMG_column = find(EMG)-2;

PhysioVars(:,3) = sum(data(:,EEG1_1to2Hzcolumn:EEG1_1to2Hzcolumn+2),2);  %the plus 2 means add the values in the columns for 1-2,2-3 and 3-4 Hz
PhysioVars(:,4) = sum(data(:,EEG2_1to2Hzcolumn:EEG2_1to2Hzcolumn+2),2);

% EMG data 
EMG_data = data(:,EMG_column);






  % PhysioVars(:,3) = mean(data(:,3:5),2);     % as many rows as there are rows in the input file, EEG1 delta power (1-4Hz) 
  % if size(data,2) == 82
  %   PhysioVars(:,4) = mean(data(:,43:45),2);   % EEG2 delta power (1-4Hz) if .txt file goes up to 40 Hz
  % elseif size(data,2) == 44
  %   PhysioVars(:,4) = mean(data(:,24:26),2); % EEG2 delta power (1-4Hz) if .txt file goes up to 20 Hz
  % end

  d1smoothed = medianfiltervectorized(PhysioVars(:,3),2); 
  d2smoothed = medianfiltervectorized(PhysioVars(:,4),2);
  
  PhysioVars(:,3) = d1smoothed;
  PhysioVars(:,4) = d2smoothed;
  


 % Handle artifacts 
  if length(find(PhysioVars(:,1)==5)) > 0
    PhysioVars = handle_artefacts(PhysioVars);
  end 

if strcmp(model,'5state')
  % re-score wake epochs into quiet wake vs. active wake, based on EMG. Wake=0,SWS=1,REM=2,quiet wake=3, active wake=4
  newstate = RescoreQuietVsActiveWake(PhysioVars(:,1),EMG_data,0.33,0.66,FileCounter,files);
  state_data{FileCounter} = newstate;
end 

if strcmp(model,'3state')
  state_data{FileCounter} = PhysioVars(:,1);
end


  if strcmp(signal,'lactate')
    signal_data{FileCounter} = PhysioVars(:,2);
  elseif strcmp(signal,'delta1') | strcmp(signal,'EEG1')
    signal_data{FileCounter} = PhysioVars(:,3);
  elseif strcmp(signal,'delta2') | strcmp(signal,'EEG2')
    signal_data{FileCounter} = PhysioVars(:,4);
  end

% Compute the dynamic range for each data file (90th percentile - 10th percentile)
  dynamic_range(FileCounter) = quantile(signal_data{FileCounter},.9)-quantile(signal_data{FileCounter},.1);

% Restrict the recording to only baseline (10AM-10AM), sleep deprivation (10AM-4PM day 2) or recovery (4PM-end)  
% first read in all the timestamp data into a matrix
% size(textdata)
% for i=1:length(textdata)
%   TimeStampMatrix{FileCounter}(:,i) = sscanf(textdata{i,1},'%f:%f:%f');
%   end

TimeStampMatrix{FileCounter} = create_TimeStampMatrix_from_textdata(textdata);
  
% if restrict is not a string, but is a two-element vector [tstart tend] instead (where tstart and tend are hours 
% from the beginning of the recording),
% restrict the recording to be only from tstart to tend. 
 

if ~isstr(restrict) & numel(restrict)==2
  dt = 1/(60*60/epoch_length_in_seconds(FileCounter));
  t  = 0:dt:dt*(size(signal_data{FileCounter},1)-1);
  if strcmp(signal,'lactate')             % include window_length/2 hours of data on either side 
    if restrict(1)-(window_length)/2 > 0      % check lower edge
      start_index = find(abs(t-(restrict(1)-(window_length/2)))<1e-12);
    else
      start_index = find(abs(t-restrict(1))<1e-12);  %where t=restrict(1) with some tolerance for round-off
    end
    if restrict(2) + (window_length/2 > t(end))   % check upper edge
      end_index = find(abs(t-(restrict(2)+(window_length/2)))<1e-12);
    else
      if restrict(2) > t(end)
        end_index = length(t);
      else 
        end_index = find(abs(t-restrict(2))<1e-12);
      end
    end
  end
  if ~strcmp(signal,'lactate')
    start_index = find(abs(t-restrict(1))<1e-12);  %where t=restrict(1) with some tolerance for round-off
    if restrict(2) > t(end)
      end_index = length(t);
    else 
      end_index = find(abs(t-restrict(2))<1e-12);
    end
  end
end

if strcmp(restrict,'baseline')
  locs_of_start_times = find([TimeStampMatrix{FileCounter}(:).Hour]==10 & [TimeStampMatrix{FileCounter}(:).Minute]==0 & [TimeStampMatrix{FileCounter}(:).Second]==0); %10AM
  start_index = locs_of_start_times(1);
  end_index = locs_of_start_times(2);  %second instance of 10AM is 24 hours into recording
end

if strcmp(restrict,'SD')
  locs_of_start_times = find([TimeStampMatrix{FileCounter}(:).Hour]==10 & [TimeStampMatrix{FileCounter}(:).Minute]==0 & [TimeStampMatrix{FileCounter}(:).Second]==0); %10AM
  start_index = locs_of_start_times(2);
  locs_of_end_times = find([TimeStampMatrix{FileCounter}(:).Hour]==16 & [TimeStampMatrix{FileCounter}(:).Minute]==0 & [TimeStampMatrix{FileCounter}(:).Second]==0); %4PM
  end_index = locs_of_end_times(2); %second instance of 4PM.  First is during baseline
end

if strcmp(restrict,'recovery')
  locs_of_start_times = find([TimeStampMatrix{FileCounter}(:).Hour]==16 & [TimeStampMatrix{FileCounter}(:).Minute]==0 & [TimeStampMatrix{FileCounter}(:).Second]==0); %4PM
  start_index = locs_of_start_times(2);  %second instance of 4PM.  First is during baseline
  end_index = size(signal_data,1);
end

if strcmp(restrict,'baseline') | strcmp(restrict,'SD') | strcmp(restrict,'recovery') | (~isstr(restrict) & numel(restrict)==2)
  state_data{FileCounter}  = state_data{FileCounter}(start_index:end_index,1);  %reset state_data and signal_data cell arrays
  signal_data{FileCounter} = signal_data{FileCounter}(start_index:end_index,1);
  TimeStampMatrix{FileCounter} = TimeStampMatrix{FileCounter}(start_index:end_index);
end 

% Now restrict the data is restrict is



% Cut off all data before 8:00PM 
 %  locs_of_start_times = find([TimeStampMatrix{FileCounter}(:).Hour]==20 & [TimeStampMatrix{FileCounter}(:).Minute]==0 & [TimeStampMatrix{FileCounter}(:).Second]==0); %the twenty is for 20:00, 8:00PM
 % start_index = locs_of_start_times(1);

 
  % state_data{FileCounter}  = state_data{FileCounter}(start_index:end,1);  %reset state_data and signal_data cell arrays to only include the data starting at 8:00PM
  % signal_data{FileCounter} = signal_data{FileCounter}(start_index:end,1);
  % TimeStampMatrix{FileCounter} = TimeStampMatrix{FileCounter}(start_index:end);


  % compute the length of the datafile in hours 
  % start_time = TimeStampMatrix{FileCounter}(:,1);
  % end_time   = TimeStampMatrix{FileCounter}(:,end);

  % start_time(1:3) = [start_time(3); start_time(1); start_time(2)];
  % end_time(1:3) = [end_time(3); end_time(1); end_time(2)];
  % length_of_recording = etime(end_time',start_time');
  % length_of_recording = length_of_recording/60/60; % convert from seconds to hours

 


% ----- 
% For lactate simulations, allow the lactate sensor to settle down:
% Find the first two NREM episodes of at least 1 minute, and start 
% the simulation at the end of the second 1-minute NREM episode
% ----
if strcmp(signal,'lactate') & (strcmp(restrict,'none') || restrict(1)==0)
   ind_of_second_NREM_episode_end = find_first_two_NREM_episodes(state_data{FileCounter});
   state_data{FileCounter}  = state_data{FileCounter}(ind_of_second_NREM_episode_end:end);
   signal_data{FileCounter} = signal_data{FileCounter}(ind_of_second_NREM_episode_end:end);
end






end % end of looping through files to load data, decide which files to exclude, and cut off lactate transient

% ------
% Exclusion criteria:
% Compute the dynamic range for each dataset and 
% include only the datasets with the 7 largest
% dynamic ranges.
%------
% [sorteddata,sortIndex]=sort(dynamic_range,'descend');
% if length(dynamic_range) >= 20
% Indices_of_largest = sortIndex(1:20);  % 20 largest dynamic ranges
% else
% Indices_of_largest = sortIndex;
% end

% % reset signal_data and state_data cell arrays to only include files that haven't been excluded 
% % by our exclusion rule
% state_data  = state_data(Indices_of_largest);
% signal_data = signal_data(Indices_of_largest);
% files       = files(Indices_of_largest);




%---
% COMPUTING LOOP
%---
% Now that I've loaded all the data and determined which datasets to keep (and simulate)
disp('Starting computing loop...')
for FileCounter=1:length(files)
  disp(['File number ', num2str(FileCounter), ' of ', num2str(length(files))])
  display(files{FileCounter})
  if strcmp(algorithm,'NelderMead')  
  [Ti,Td,LA,UA,best_error,error_instant,S,ElapsedTime] = Franken_like_model_with_nelder_mead([state_data{FileCounter} signal_data{FileCounter}],signal,files{FileCounter},model,epoch_length_in_seconds(FileCounter),window_length);
  end
  if strcmp(algorithm,'BruteForce')
  [Ti,Td,LA,UA,best_error,error_instant,S,ElapsedTime] = Franken_like_model([state_data{FileCounter} signal_data{FileCounter}],signal,files{FileCounter},model,epoch_length_in_seconds(FileCounter),window_length); %for brute-force 
  end

residual(FileCounter) = best_error;

%[LAnormalized,UAnormalized]=normalizeLAUA(LA,UA,[state_data{FileCounter} signal_data{FileCounter}],TimeStampMatrix{FileCounter});
LAnormalized = LA;
UAnormalized = UA;
  Taui(FileCounter) = Ti;
  TauD(FileCounter) = Td;
  

  LowA{FileCounter} = LAnormalized;  %normalized LA
  UppA{FileCounter} = UAnormalized;  %normalized UA
  Timer(FileCounter) = ElapsedTime;
  Error(FileCounter) = best_error;
  Error2(FileCounter) = error_instant;
  %delete(findall(0,'Type','figure')); %if you want to delete all figures before next run
  
  % populate cell array with S output
  best_S{FileCounter} = S;

 

end  %end of looping through files
    
%write a DataSourceInfo tab
% disp('Writing data to an Excel file...')
% addpath ../../../../../Brennecke/matlab-pipeline/Matlab/etc/matlab-utils/;
% xl=XL;
% sheet = xl.Sheets.Item(1);

% col1_name{1} = 'Taui Values';                                       % column headers
% first_header_cells=xl.setCells(sheet,[1,1],col1_name,'false','true');
% set(first_header_cells.Font, 'Bold', true)
% col2_name{1} = 'Taud Values';
% second_header_cells=xl.setCells(sheet,[2,1],col2_name,'false','true');
% set(second_header_cells.Font, 'Bold', true)
% col3_name{1} = 'CPU time';
% third_header_cells=xl.setCells(sheet,[3,1],col3_name,'false','true');
% set(third_header_cells.Font, 'Bold', true)


% xl.setCells(sheet,[1,2],Taui','false','false');

% xl.sourceInfo(mfilename('fullpath'));                                % DataSourceInfo tab
% xl.saveAs(strcat('PROCESSLBATCHMODE output', '.txt'));



% make a bar graph showing the error in the model and the error using 
% the instant model that follows the lactate upper and lower bounds and 
% switches in-between them if the sleep state changes
% if strcmp(signal,'lactate')
%   figure
%   bar([1,2],[mean(Error)/mean(Error) mean(Error2)/mean(Error)])
%   set(gca,'XTickLabel',{'Model fit','instant model'})
%   hold on
%   h=errorbar([1 2],[mean(Error)/mean(Error) mean(Error2)/mean(Error)],[std(Error)/sqrt(length(Error)) std(Error2)/sqrt(length(Error2))]); 
%   %set(h(2),'LineStyle','none','Marker','s','MarkerEdgeColor','k')
  
%   d=daspect;
%   daspect([d(1)*2 d(2) d(3)])
  
%   title('Error between model fit and instant model that follows UA or LA')
% end

if do_write_tau_values_to_file
  % set up directory for output files (if making them)
  date_time = datestr(now,'mm.dd.yyyy.hh.MM');
  output_directory = strcat(directory,'Tau_values_output_',date_time);
  mkdir(output_directory)

  write_tau_values_to_file(files,directory,model,signal,algorithm,do_rescore,do_restrict,restrict,Taui,TauD)
end




load chirp % this assigns chirp to the variable y
sound  (y)





%profile viewer