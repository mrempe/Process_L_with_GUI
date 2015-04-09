function [Ti,Td,LAnormalized,UAnormalized,best_error,error_instant,best_S,ElapsedTime]=Franken_like_model(datafile,signal,filename,epoch_length,window_length)
% USAGE:  [Ti,Td,error]=Franken_like_model(datafile,signal)
%
% datafile: a sleep data file from Jonathan Wisor where sleep
%           state is in the first column, lactate in the second column and
%           EEG data in the columns after that.  
%
% signal: either 'delta' or 'lactate' 
%
% filename: the name of the .txt data file so I can use it in figure titles, etc.
%
% epoch_length: length of the epoch in seconds 
%
% window_length: length of moving window (in hours) used to compute UA and LA if signal is lactate
%
% OUTPUT:
% Ti: the optimum value for tau_i, the rate of increase, using a
% two-process-like model, similar to Franken et al 2001
%
% Td: the optimum value for tau_d, the decay rate. 
% 
% error: the mean square error for the best fit
tic

%window_length=4; 

% make a frequency plot, and use it to figure out upper and lower
% bounds for the model (like Franken et al. 2001 Figure 1)
[LA,UA]=make_frequency_plot(datafile,window_length,signal,epoch_length);


% -- if using delta power normalize UA and LA to mean SWS delta 
% -- power in last 4 hours of baseline light period and find 
% -- all SWS episodes of longer than 5 minutes (like 
% -- Franken et al)
%if strcmp(signal,'delta1') || strcmp(signal,'delta2')
  baseline_start_hours = 17;
  baseline_end_hours = 21;
  ind_start = baseline_start_hours*(60*60/epoch_length);
  ind_end = baseline_end_hours*(60*60/epoch_length);
  
  %[t_mdpt_SWS,data_at_SWS_midpoints,t_mdpt_indices]=find_all_SWS_episodes2(datafile);

  locs = find(datafile(ind_start:ind_end,1)==1); % find SWS epochs in last 4 hr of baseline
  mn   = mean(datafile(locs+ind_start-1,2));     % mean delta power during SWS in last 4hr of baseline

  LAnormalized = (LA/mn)*100;   % lower asymptote normalized to mean delta power during SWS in last 4hr of baseline
  UAnormalized = (UA/mn)*100;   % upper asymptote normalized to mean delta power during SWS in last 4hr of baseline

%end


% if using delta power as a signal, prepare the data we will compare 
% to by finding all SWS episodes of longer than 5 minutes (like 
% Franken et al)
if strcmp(signal,'delta1') || strcmp(signal,'delta2') || strcmp(signal,'EEG1') || strcmp(signal,'EEG2')
  [t_mdpt_SWS,data_at_SWS_midpoints,t_mdpt_indices]=find_all_SWS_episodes2(datafile,epoch_length);
end

% if using a moving window for the upper and lower assymptotes, S
% will have 720 fewer elements than the number of rows of datafile,
% so set up a new index for S
% mask=find(t_mdpt_indices>(60*60/epoch_length) & t_mdpt_indices<(size(datafile,1)-(60*60/epoch_length)));
% t_mdpt_SWS_moving_window=t_mdpt_SWS(mask);
% data_at_SWS_midpoints_moving_window=data_at_SWS_midpoints(mask);
% t_mdpt_indices_moving_window=t_mdpt_indices(mask);
%mask=361:size(datafile,1)-(60*60/epoch_length)';
mask=(window_length/2)*(60*60/epoch_length)+1:size(datafile,1)-(window_length/2)*(60*60/epoch_length);


dt=1/(60*60/epoch_length);  % assuming data points are every 10 seconds and t is in hours 
if strcmp(signal,'delta1') || strcmp(signal,'delta2') || strcmp(signal,'EEG1') || strcmp(signal,'EEG2')
tau_i=0.05:.1:5; %1:.12:25; %following Franken (was 0.05:.1:5)
tau_d=0.05:0.05:5; %0.1:.025:5; %following Franken (was 0.05:0.025:5)
elseif strcmp(signal,'lactate')
tau_i=linspace(0.01,2,201); %0.01:.005:2;  %1:.12:25  % make sure these vectors have the same length as tau_i and tau_d for delta
tau_d=linspace(0.01,2,197); %0.01:0.005:2; %0.1:.025:5
end

error=zeros(length(tau_i),length(tau_d));

% COMPUTING LOOP
% run the model and compute error for all combinations of tau_i and tau_d
for i=1:length(tau_i)
  for j=1:length(tau_d)
   
    S=run_S_model(datafile,dt,(LA(1)+UA(1))/2,LA,UA,tau_i(i),tau_d(j),window_length,0,epoch_length,filename); % run model
   
    % compute error (depending on if delta power or lactate was used)
    if strcmp(signal,'delta1') || strcmp(signal,'delta2') || strcmp(signal,'EEG1') || strcmp(signal,'EEG2')
      error(i,j)=sqrt((sum((S([t_mdpt_indices])-data_at_SWS_midpoints).^2))/length(t_mdpt_indices)); %RMSE
    elseif strcmp(signal,'lactate')
      error(i,j)=sqrt((sum((S'-datafile([mask],2)).^2))/(size(datafile,1)-(window_length*(60*60/epoch_length)))); %RMSE
    end
      
    % display progress only at intervals of .25*total 
    display_progress(length(tau_d)*(i-1)+j,length(tau_i)*length(tau_d));

  end
end

best_error=min(min(error));
[r,c]=find(error==min(min(error)));
Ti=tau_i(r);
Td=tau_d(c);

% run one more time with best fit and plot it (add a plot with circles)
if  strcmp(signal,'lactate')
  best_S=run_S_model(datafile,dt,(LA(1)+UA(1))/2,LA,UA,tau_i(r),tau_d(c),window_length,1,epoch_length,filename);
%error_instant=run_instant_model(datafile,LA,UA,window_length);
  error_instant = 0;
end

if strcmp(signal,'delta1') || strcmp(signal,'delta2') || strcmp(signal,'EEG1') || strcmp(signal,'EEG2')
  best_S=run_S_model(datafile,dt,(LA(1)+UA(1))/2,LA,UA,Ti,Td,window_length,0,epoch_length,filename);
error_instant = 0;
end


figure
t=0:dt:dt*(size(datafile,1)-1);
if strcmp(signal,'delta1') | strcmp(signal,'delta2') || strcmp(signal,'EEG1') || strcmp(signal,'EEG2')
  plot(t_mdpt_SWS,data_at_SWS_midpoints,'ro')
  hold on
  plot(t,best_S)
  ylabel('Delta power')
  title(['Best fit of model to delta power data for file ' filename])
elseif strcmp(signal,'lactate')
  plot(t,datafile(:,2),'ro')
  hold on
  tS=t((window_length/2)*(60*60/epoch_length)+1:end-(window_length/2)*(60*60/epoch_length));
  plot(tS,best_S)
  plot(tS,LA,'--')
  plot(tS,UA,'--')
  ylabel('lactate')
  title('Best fit of model to lactate data')
end
xlabel('Time (hours)')
hold off

ElapsedTime=toc

% hold on

% if strcmp(signal,'delta') 
%   plot(t_mdpt_SWS,data_at_SWS_midpoints,'ro')
% ylabel('Delta power')
% elseif strcmp(signal,'lactate')
%   plot(t,datafile(:,2),'ro')
% ylabel('lactate')
% end

% hold off
% if strcmp(signal,'delta')
%   title('Best fit of model to delta power data')
% elseif strcmp(signal,'lactate')
%   title('Best fit of model to lactate data')
% end  
% xlabel('Time (hours)')



% make a contour plot of the errors
figure
[X,Y]=meshgrid(tau_d,tau_i);
contour(X,Y,error,100)
ylabel('\tau_i')
xlabel('\tau_d')
%colorbar

