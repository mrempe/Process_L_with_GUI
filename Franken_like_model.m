function [Ti,Td,LA,UA,best_error,error_instant,best_S,ElapsedTime]=Franken_like_model(datafile,timestampvec,signal,filename,model,epoch_length,window_length)
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
%[LA,UA]=make_frequency_plot(datafile,window_length,signal,epoch_length);
[LA,UA]=make_frequency_plot(datafile,window_length,signal,timestampvec,tL,epoch_length,0,0);

% -- if using delta power normalize UA and LA to mean SWS delta 
% -- power in last 4 hours of baseline light period and find 
% -- all SWS episodes of longer than 5 minutes (like 
% -- Franken et al)
%if strcmp(signal,'delta1') || strcmp(signal,'delta2')
  % baseline_start_hours = 17;
  % baseline_end_hours = 21;
  % ind_start = baseline_start_hours*(60*60/epoch_length);
  % ind_end = baseline_end_hours*(60*60/epoch_length);
  
  % %[t_mdpt_SWS,data_at_SWS_midpoints,t_mdpt_indices]=find_all_SWS_episodes2(datafile);

  % locs = find(datafile(ind_start:ind_end,1)==1); % find SWS epochs in last 4 hr of baseline
  % mn   = mean(datafile(locs+ind_start-1,2));     % mean delta power during SWS in last 4hr of baseline

  % LAnormalized = (LA/mn)*100;   % lower asymptote normalized to mean delta power during SWS in last 4hr of baseline
  % UAnormalized = (UA/mn)*100;   % upper asymptote normalized to mean delta power during SWS in last 4hr of baseline

%end


% if using delta power as a signal, prepare the data we will compare 
% to by finding all SWS episodes of longer than 5 minutes (like 
% Franken et al)
if strcmp(signal,'delta1') || strcmp(signal,'delta2') || strcmp(signal,'EEG1') || strcmp(signal,'EEG2')
 % [t_mdpt_SWS,data_at_SWS_midpoints,t_mdpt_indices]=find_all_SWS_episodes2(datafile,epoch_length);
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
   
    S=run_S_model(datafile,dt,(LA(1)+UA(1))/2,LA,UA,tau_i(i),tau_d(j),window_length,0,timestampvec,tL,epoch_length); % run model
    %S=run_S_model(datafile,dt,(LA(1)+UA(1))/2,LA,UA,p(1),p(2),0,0,timestampvec,tL,epoch_length);
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
    
   
    scatter(t(only_wake_indices),wake_lactate,25,'r','filled')
    hold on
    scatter(t(only_sleep_indices),sleep_lactate,25,'k','filled')
    scatter(t(only_rem_indices),rem_lactate,25,'g','filled')
   

    if strcmp(model,'5state')
      only_quiet_wake_indices  = find(datafile(:,1)==3);
      only_active_wake_indices = find(datafile(:,1)==4);
      quiet_wake_lactate  = datafile(only_quiet_wake_indices,2);
      active_wake_lactate = datafile(only_active_wake_indices,2);
      scatter(t(only_quiet_wake_indices),quiet_wake_lactate,25,[1 0.5 0],'filled')  % orange
      scatter(t(only_active_wake_indices),active_wake_lactate,25,[0.67 0.45 0.2],'filled') % brown
    end
    
   %tS=t(361:end-(60*60/epoch_length));
    tS=t((window_length/2)*(60*60/epoch_length)+1:end-(window_length/2)*(60*60/epoch_length));
    
    plot(tS,best_S,'b','LineWidth',1.5)
    plot(tS,LA,'k--')
    plot(tS,UA,'k--')
    ylabel('lactate')
    xlabel('Time (hours)')
    title(['Best fit of model to lactate data for file ' filename 'using ' num2str(epoch_length) '-second epochs'])
    hold off
    
    if strcmp(model,'3state')
      legend('Wake','SWS','REMS','Model')
    end

    if strcmp(model,'5state')
      legend('Wake','SWS','REMS','QW','AW','Model')
    end
    

    figure
    L_indices = (window_length/2)*(60*60/epoch_length)+1:length(t)-(window_length/2)*(60*60/epoch_length);
    scaled_lactate_data = ((UA-LA)-(UA-datafile(L_indices,2)'))./(UA-LA);
    only_sleep_indices_L = find(datafile(L_indices,1)==1);
    only_wake_indices_L  = find(datafile(L_indices,1)==0);
    only_rem_indices_L   = find(datafile(L_indices,1)==2);


    sleep_lactate_scaled=scaled_lactate_data(only_sleep_indices_L);
    wake_lactate_scaled=scaled_lactate_data(only_wake_indices_L);
    rem_lactate_scaled=scaled_lactate_data(only_rem_indices_L);
    scatter(tS(only_wake_indices_L),wake_lactate_scaled,25,'rd')
    hold on
    scatter(tS(only_sleep_indices_L),sleep_lactate_scaled,25,'b+')
    scatter(tS(only_rem_indices_L),rem_lactate_scaled,25,'gx')

  if strcmp(model,'5state')
    only_quiet_wake_indices_L  = find(datafile(L_indices,1)==3);
    only_active_wake_indices_L = find(datafile(L_indices,1)==4);
    quiet_wake_lactate_scaled=scaled_lactate_data(only_quiet_wake_indices_L);
    active_wake_lactate_scaled=scaled_lactate_data(only_active_wake_indices_L);
    scatter(tS(only_active_wake_indices_L),active_wake_lactate_scaled,25,[0.67 0.45 0.2],'o')
    scatter(tS(only_quiet_wake_indices_L),quiet_wake_lactate_scaled,25,[1 0.5 0],'s')
  end



    scaled_L = ((UA-LA)-(UA-best_S))./(UA-LA);
    %plot(tS,scaled_lactate_data,'ro')
    
    plot(tS,scaled_L,'b','LineWidth',1.5)
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
figure
[X,Y]=meshgrid(tau_d,tau_i);
contour(X,Y,error,100)
ylabel('\tau_i')
xlabel('\tau_d')
colorbar
hold on
plot(Td,Ti,'rx')
hold off


