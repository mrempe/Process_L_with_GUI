function [Ti,Td,LA,UA,best_error,error_instant,best_S,ElapsedTime]=Franken_like_model_with_nelder_mead(datafile,timestampvec,signal,filename,model,epoch_length,window_length)
% USAGE:  [Ti,Td,LA,UA,error]=Franken_like_model_with_nelder_mead(datafile,signal)
%
% datafile: a sleep data file from Jonathan Wisor where sleep
%           state is in the first column, lactate or EEG data in the second column 
%
% timestampvec: a vector of timestamps made using datetime.  With epochs removed for artifacts and negative lactate values.  
% signal: 'delta1' or 'delta2' or 'EEG1' or 'EEG2' or 'lactate' 
%
% filename: the name of the .txt data file so I can use it in figure titles, etc.
%
% model: 5state or 3state depending on whether the model includes 3 states (W,SWS,REMS) or 5 (AW,QW,W,SWS,REMS)
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



%% ---- Set up the time vector for lactate since you ignore data ------ %%
%% ---- at the beginning and the end (due to the moving average) ------ %%
% being careful to think about artifacts and allowing the sensor to settle down.  
if strcmp(signal,'lactate')
  tL_start_time  = timestampvec(1) + hours(window_length/2);
  tL_end_time    = timestampvec(end) - hours(window_length/2);
  indices_start = find(timestampvec>=tL_start_time);  % this handles the case where the start time has been left out due to artifact
  tL_start_index = indices_start(1);
  indices_end = find(timestampvec>=tL_end_time);
  tL_end_index   = indices_end(1);
  tL = timestampvec(tL_start_index:tL_end_index);
else 
  tL =0;
end
% make a frequency plot, and use it to figure out upper and lower
% bounds for the model (like Franken et al. 2001 Figure 1)
[LA,UA]=make_frequency_plot(datafile,window_length,signal,timestampvec,tL,epoch_length,0,0);


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
  [t_mdpt_SWS,data_at_SWS_midpoints,t_mdpt_indices]=find_all_SWS_episodes5(datafile,timestampvec,epoch_length);
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



%mask=(window_length/2)*(60*60/epoch_length)+1:size(datafile,1)-(window_length/2)*(60*60/epoch_length);
if strcmp(signal,'lactate')
  size(timestampvec)
  tL(1)
  tL(end)
  mask=find(timestampvec==tL(1)):find(timestampvec==tL(end));
else
  mask=0;
end

%dt=1/(60*60/epoch_length);  % assuming  t is in hours 
dt=timestampvec(2:end)-timestampvec(1:end-1);
dt=hours(dt);  % convert dt into hours to use in run_S_model
dt=min(dt);        % only want one value, not a vector. Take the minimum in case there is missing data
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
								datafile,dt,(LA(1)+UA(1))/2,LA,UA,window_length,timestampvec,0,epoch_length,mask,p),initial_guess_delta,optimset('TolX',1e-3));
end

if strcmp(signal,'lactate')
% set up the time vector for lactate since you ignore data at the beginning and the end (due to the moving average)
% being careful to think about artifacts and allowing the sensor to settle down.  
  % tL_start_time  = timestampvec(1) + hours(window_length/2);
  % tL_end_time    = timestampvec(end) - hours(window_length/2);
  % tL_start_index = find(timestampvec==tL_start_time);
  % tL_end_index   = find(timestampvec==tL_end_time);
  % tL = timestampvec(tL_start_index:tL_end_index);
  S_start_index = find(timestampvec==tL(1));
  S0 = datafile(S_start_index,2);  % lactate data 
  [bestparams,best_error] = fminsearch(@(p) myobjectivefunction(signal,0,0,datafile,dt,S0,LA,UA, ...
								window_length,timestampvec,tL,epoch_length,mask,p),initial_guess_lactate,optimset('TolX',1e-3));
end
best_tau_i=bestparams(1);
best_tau_d=bestparams(2);

Ti=best_tau_i;    %output the best taus
Td=best_tau_d;


% run one more time with best fit and plot it (with data points)
if  strcmp(signal,'lactate')
  best_S=run_S_model(datafile,dt,(LA(1)+UA(1))/2,LA,UA,Ti,Td,window_length,0,timestampvec,tL,epoch_length,filename);
  %error_instant=run_instant_model(datafile,LA,UA,window_length);
error_instant = 0;
end
if strcmp(signal,'delta1') || strcmp(signal,'delta2') || strcmp(signal,'EEG1') || strcmp(signal,'EEG2')
 best_S=run_S_model(datafile,dt,(LA(1)+UA(1))/2,LA,UA,Ti,Td,window_length,0,timestampvec,tL,epoch_length,filename);
end


ElapsedTime=toc;
disp(['ElapsedTime = ', num2str(ElapsedTime), ' seconds.'])

% plot the best fit
%t=0:dt:dt*(size(datafile,1)-1);


if strcmp(signal,'delta1') || strcmp(signal,'delta2') || strcmp(signal,'EEG1') || strcmp(signal,'EEG2')
  error_instant=0;  % this won't get set if signal is delta, but the function returns it
  figure
  %only_sleep_indices=find(datafile(:,1)==1);  
  %sleep_eeg1=datafile(only_sleep_indices,3);
  %sleep_eeg2=datafile(only_sleep_indices,4);
  %scatter(t(only_sleep_indices),sleep_eeg2,25,'r')
 
  plot(t_mdpt_SWS,data_at_SWS_midpoints,'go')
  hold on
  %plot(t,best_S)
  plot(timestampvec,best_S)
  xmin = datenum(timestampvec(1));
  xmax = datenum(timestampvec(end));
  xlim([xmin xmax])
  ylabel('Delta power')
  xlabel('Time')
  title(['Best fit of model to delta power data for file ' filename ' using ' num2str(epoch_length) '-second epochs and a ' model(1) '-state model' ])
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
    
    %scatter(timestampvec(only_wake_indices),wake_lactate,25,'r','filled')
    plot(timestampvec(only_wake_indices),wake_lactate,'r.','MarkerSize',16)
    hold on
    plot(timestampvec(only_sleep_indices),sleep_lactate,'k.','MarkerSize',16)
    plot(timestampvec(only_rem_indices),rem_lactate,'g.','MarkerSize',16)
    % scatter(timestampvec(only_sleep_indices),sleep_lactate,25,'k','filled')
    % scatter(timestampvec(only_rem_indices),rem_lactate,25,'g','filled')
   

    if strcmp(model,'5state')
      only_quiet_wake_indices  = find(datafile(:,1)==3);
      only_active_wake_indices = find(datafile(:,1)==4);
      quiet_wake_lactate  = datafile(only_quiet_wake_indices,2);
      active_wake_lactate = datafile(only_active_wake_indices,2);
      % scatter(timestampvec(only_quiet_wake_indices),quiet_wake_lactate,25,[1 0.5 0],'filled')  % orange
      % scatter(timestampvec(only_active_wake_indices),active_wake_lactate,25,[0.67 0.45 0.2],'filled') % brown
      plot(timestampvec(only_quiet_wake_indices),quiet_wake_lactate,'.','MarkerFaceColor',[1 0.5 0],'MarkerSize',16)  % orange
      plot(timestampvec(only_active_wake_indices),active_wake_lactate,'.','MarkerFaceColor',[0.67 0.45 0.2],'MarkerSize',16) % brown
    end
    
    % Set up the vector to plot the simulation output using lactate as the signal
    % You must remove a chunk of time at the beginning and ending of timestampvec
    % since we're using a moving average of the lactate data
    %tS=t(361:end-(60*60/epoch_length));
    %tS=t((window_length/2)*(60*60/epoch_length)+1:end-(window_length/2)*(60*60/epoch_length));
    % tL_start_time  = timestampvec(1) + hours(window_length/2);
    % tL_end_time    = timestampvec(end) - hours(window_length/2);
    % tL_start_index = find(timestampvec==tL_start_time);
    % tL_end_index   = find(timestampvec==tL_end_time);
    % tL = timestampvec(tL_start_index:tL_end_index);
    % timestampvec(end-10:end)
    % length(best_S)
    % length(tL)
    plot(tL,best_S,'b','LineWidth',1.5)
    plot(tL,LA,'k--')
    plot(tL,UA,'k--')
    xmin = datenum(timestampvec(1));
    xmax = datenum(timestampvec(end));
    xlim([xmin xmax])
    ylabel('lactate')
    xlabel('Time (hours)')
    tick_locations=datenum(timestampvec(1:60/epoch_length*60*4:end)); % ticks every 4 hours
    set(gca,'XTick',tick_locations)
    title(['Best fit of model to lactate data for file ' filename 'using ' num2str(epoch_length) '-second epochs'])
    hold off
    
    if strcmp(model,'3state')
      legend('Wake','SWS','REMS','Model')
    end

    if strcmp(model,'5state')
      legend('Wake','SWS','REMS','QW','AW','Model')
    end
    

    figure
    %L_indices = (window_length/2)*(60*60/epoch_length)+1:length(t)-(window_length/2)*(60*60/epoch_length);
    %L_indices = find(timestampvec>=tL(1) & timestampvec<=tL(end));
    L_indices = tL_start_index:tL_end_index;
    size(L_indices)
    size(UA)
    size(LA)
    scaled_lactate_data = ((UA-LA)-(UA-datafile(L_indices,2)'))./(UA-LA);
    only_sleep_indices_L = find(datafile(L_indices,1)==1);
    only_wake_indices_L  = find(datafile(L_indices,1)==0);
    only_rem_indices_L   = find(datafile(L_indices,1)==2);


    sleep_lactate_scaled=scaled_lactate_data(only_sleep_indices_L);
    wake_lactate_scaled=scaled_lactate_data(only_wake_indices_L);
    rem_lactate_scaled=scaled_lactate_data(only_rem_indices_L);
    plot(tL(only_wake_indices_L),wake_lactate_scaled,'rd','MarkerSize',4)
    hold on
    plot(tL(only_sleep_indices_L),sleep_lactate_scaled,'b+','MarkerSize',4)
    plot(tL(only_rem_indices_L),rem_lactate_scaled,'gx','MarkerSize',4)

  if strcmp(model,'5state')
    only_quiet_wake_indices_L  = find(datafile(L_indices,1)==3);
    only_active_wake_indices_L = find(datafile(L_indices,1)==4);
    quiet_wake_lactate_scaled=scaled_lactate_data(only_quiet_wake_indices_L);
    active_wake_lactate_scaled=scaled_lactate_data(only_active_wake_indices_L);
    plot(tL(only_active_wake_indices_L),active_wake_lactate_scaled,'o','MarkerEdgeColor',[0.67 0.45 0.2])
    plot(tL(only_quiet_wake_indices_L),quiet_wake_lactate_scaled,'s','MarkerEdgeColor',[1 0.5 0])
  end



    scaled_L = ((UA-LA)-(UA-best_S))./(UA-LA);
    %plot(tS,scaled_lactate_data,'ro')
    
    plot(tL,scaled_L,'b','LineWidth',1.5)
    hold off
    ylabel('Lactate (scaled)')
    xlabel('Time (hours)')
    title(['Scaled lactate data for file ' filename 'using ' num2str(epoch_length) '-second epochs'])

    if strcmp(model,'3state')
      legend('Wake','SWS','REMS','Model')
    end

    if strcmp(model,'5state')
      legend('Wake','SWS','REMS','QW','AW','Model')
    end



     best_S = scaled_L;   % Set best_S to the scaled lactate output so the function returns scaled_L

end

% make a contour plot of the errors
% figure
% [X,Y]=meshgrid(tau_d,tau_i);
% contour(X,Y,error,100)
% ylabel('\tau_i')
% xlabel('\tau_d')
%colorbar


