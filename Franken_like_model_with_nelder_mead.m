function [Ti,Td,LA,UA,best_error,error_instant,best_S,ElapsedTime]=Franken_like_model_with_nelder_mead(datafile,signal,filename,epoch_length,window_length)
% USAGE:  [Ti,Td,LA,UA,error]=Franken_like_model_with_nelder_mead(datafile,signal)
%
% datafile: a sleep data file from Jonathan Wisor where sleep
%           state is in the first column, lactate or EEG data in the second column 
%
% signal: 'delta1' or 'delta2' or 'EEG1' or 'EEG2' or 'lactate' 
%
% filename: the name of the .txt data file so I can use it in figure titles, etc.
%
% epoch_length: length of the epoch in seconds 
%
% window_length: length of moving window (in hours) to use when computing LA and UA for lactate signal
%
% OUTPUT:
% Ti: the optimum value for tau_i, the rate of increase, using a
% two-process-like model, similar to Franken et al 2001
%
% Td: the optimum value for tau_d, the decay rate. 
% 
% LA: the lower asymptote found from make_frequency_plot.m (not nomralized)
% 
% UA: the upper asymptote found from make_frequency_plot.m (not normalized)
%
% best_error: the mean square error for the best fit
%
% error_instant: the error of the model run with lactate as the 
% signal where instead of following the exponential model it 
% just instantly jumps up or down to the UA or LA if the state
% switches.  
%
% best_S: the vector S of the best-fit homeostatic model
% NOTE: in PROCESSLBATCHMODE.m LA and UA get normalized by the mean SWS delta power 
% in the last 4 hours of the baseline light period

tic


%window_length=4;  % size of moving window (in hours) used to compute
                  % the upper and lower asymptotes for the L model.  


% make a frequency plot, and use it to figure out upper and lower
% bounds for the model (like Franken et al. 2001 Figure 1)

[LA,UA]=make_frequency_plot(datafile,window_length,signal,epoch_length);


% -- if using delta power normalize UA and LA to mean SWS delta 
% -- power in last 4 hours of baseline light period and find 
% -- all SWS episodes of longer than 5 minutes (like 
% -- Franken et al)
  % baseline_start_hours = 17;   % These values are only good if the dataset starts at 8:00 PM
  % baseline_end_hours = 21;
  % ind_start = baseline_start_hours*(60*60/epoch_length)
  % ind_end = baseline_end_hours*(60*60/epoch_length)
  
  % locs = find(datafile(ind_start:ind_end,1)==1); % find SWS epochs in last 4 hr of baseline
  % mn   = mean(datafile(locs+ind_start-1,2));     % mean delta power during SWS in last 4hr of baseline

  % LAnormalized = (LA/mn)*100;   % lower asymptote normalized to mean delta power during SWS in last 4hr of baseline
  % UAnormalized = (UA/mn)*100;   % upper asymptote normalized to mean delta power during SWS in last 4hr of baseline


if strcmp(signal,'delta1') || strcmp(signal,'delta2') || strcmp(signal,'EEG1') || strcmp(signal,'EEG2')
  [t_mdpt_SWS,data_at_SWS_midpoints,t_mdpt_indices]=find_all_SWS_episodes2(datafile,epoch_length);
disp(['Average delta power: ' num2str(mean(data_at_SWS_midpoints))])
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


dt=1/(60*60/epoch_length);  % assuming  t is in hours 
% tau_i=[0.05:0.01:1 1.1:.5:5];  %1:.12:25
% tau_d=[0.05:0.01:1 1.1:.5:5]; %0.1:.025:5
% error=zeros(length(tau_i),length(tau_d));

% COMPUTING LOOP
% use the nelder_mead algorithm to find the global minimum error
% fminsearch uses Nelder-Mead
initial_guess_delta = [2 2];     % one starting guess
initial_guess_lactate = [0.3 0.3];

if strcmp(signal,'delta1') || strcmp(signal,'delta2') || strcmp(signal,'EEG1') || strcmp(signal,'EEG2')
  [bestparams,best_error] = fminsearch(@(p) myobjectivefunction(signal,t_mdpt_indices,data_at_SWS_midpoints, ...
								datafile,dt,LA,UA,window_length,epoch_length,mask,p),initial_guess_delta,optimset('TolX',1e-3));
end

if strcmp(signal,'lactate')
[bestparams,best_error] = fminsearch(@(p) myobjectivefunction(signal,0,0,datafile,dt,LA,UA, ...
								window_length,epoch_length,mask,p),initial_guess_lactate,optimset('TolX',1e-3));
end
best_tau_i=bestparams(1);
best_tau_d=bestparams(2);

Ti=best_tau_i;    %output the best taus
Td=best_tau_d;


% run one more time with best fit and plot it (add a plot with circles)
if  strcmp(signal,'lactate')
  best_S=run_S_model(datafile,dt,(LA(1)+UA(1))/2,LA,UA,Ti,Td,window_length,0,epoch_length,filename);
  %error_instant=run_instant_model(datafile,LA,UA,window_length);
error_instant = 0;
end
if strcmp(signal,'delta1') || strcmp(signal,'delta2') || strcmp(signal,'EEG1') || strcmp(signal,'EEG2')
 best_S=run_S_model(datafile,dt,(LA(1)+UA(1))/2,LA,UA,Ti,Td,window_length,0,epoch_length,filename);
end


ElapsedTime=toc


% plot the best fit
t=0:dt:dt*(size(datafile,1)-1);


if strcmp(signal,'delta1') || strcmp(signal,'delta2') || strcmp(signal,'EEG1') || strcmp(signal,'EEG2')
  error_instant=0;  % this won't get set if signal is delta, but the function returns it
  figure
  %only_sleep_indices=find(datafile(:,1)==1);  
  %sleep_eeg1=datafile(only_sleep_indices,3);
  %sleep_eeg2=datafile(only_sleep_indices,4);
  %scatter(t(only_sleep_indices),sleep_eeg2,25,'r')
  plot(t_mdpt_SWS,data_at_SWS_midpoints,'go')
  hold on
  plot(t,best_S)
  ylabel('Delta power')
  xlabel('Time (hours)')
  title(['Best fit of model to delta power data for file ' filename ' using ' num2str(epoch_length) '-second epochs'])
    hold off


  elseif strcmp(signal,'lactate')
    figure
    %plot(t,datafile(:,2),'ro')
    only_sleep_indices = find(datafile(:,1)==1);
    only_wake_indices  = find(datafile(:,1)==0);
    only_rem_indices   = find(datafile(:,1)==2);
    sleep_lactate=datafile(only_sleep_indices,2);
    wake_lactate=datafile(only_wake_indices,2);
    rem_lactate=datafile(only_rem_indices,2);
    scatter(t(only_wake_indices),wake_lactate,25,'r')
      
    hold on
    scatter(t(only_sleep_indices),sleep_lactate,25,'k')
    scatter(t(only_rem_indices),rem_lactate,25,'c')
      
    
   %tS=t(361:end-(60*60/epoch_length));
    tS=t((window_length/2)*(60*60/epoch_length)+1:end-(window_length/2)*(60*60/epoch_length));
    
    plot(tS,best_S,'b')
    plot(tS,LA,'k--')
    plot(tS,UA,'k--')
    ylabel('lactate')
    xlabel('Time (hours)')
    title(['Best fit of model to lactate data for file ' filename 'using ' num2str(epoch_length) '-second epochs'])
    hold off
    

    figure
    L_indices = (window_length/2)*(60*60/epoch_length)+1:length(t)-(window_length/2)*(60*60/epoch_length);
    scaled_lactate_data = ((UA-LA)-(UA-datafile(L_indices,2)'))./(UA-LA);
    only_sleep_indices_L = find(datafile(L_indices,1)==1);
    only_wake_indices_L  = find(datafile(L_indices,1)==0);
    only_rem_indices_L   = find(datafile(L_indices,1)==2);


    sleep_lactate_scaled=scaled_lactate_data(only_sleep_indices_L);
    wake_lactate_scaled=scaled_lactate_data(only_wake_indices_L);
    rem_lactate_scaled=scaled_lactate_data(only_rem_indices_L);
    scatter(tS(only_wake_indices_L),wake_lactate_scaled,25,'r')
    hold on
    scatter(tS(only_sleep_indices_L),sleep_lactate_scaled,25,'k')
    scatter(tS(only_rem_indices_L),rem_lactate_scaled,25,'c')


    scaled_L = ((UA-LA)-(UA-best_S))./(UA-LA);
    %plot(tS,scaled_lactate_data,'ro')
    
    plot(tS,scaled_L,'b')
    hold off
    ylabel('Lactate (scaled)')
    xlabel('Time (hours)')
    title(['Scaled lactate data for file ' filename 'using ' num2str(epoch_length) '-second epochs'])

     best_S = scaled_L;   % Set best_S to the scaled lactate output so the function returns scaled_L

  end

% make a contour plot of the errors
% figure
% [X,Y]=meshgrid(tau_d,tau_i);
% contour(X,Y,error,100)
% ylabel('\tau_i')
% xlabel('\tau_d')
%colorbar


