function S=run_S_model(dataset,dt,S0,LA,UA,ti,td,window_length,makeplot,timestampvec,tL,epoch_length,filename)

%usage: S=run_S_model(dataset,S0,LA,UA,ti,td)
% dataset contains 2 columns. sleep state in the
% first colum (0 for wake, 1 for sleep, 2 for REM), lactate or 
% delta power in the second column.
%
% dt: is the epoch length (converted to hours)
%
% S0: the starting value for S
%
% LA: the value of the lower asymptote (found from
%     make_frequency_plot.m)
%
% UA: the value of the upper asymptote (found from
%     make_frequency_plot.m)
%
% ti: the time constant for the increasing exponential for wake, active wake, and REM
%
% td: the time constant for the decreasing exponential
%
% window_length: is using lactate as the signal, use a moving
% window of length window_length (in hours) to compute LA and UA.
% they are computed as the 90th and 10th percentile for a 2 hour
% moving window centered at each data point. 
%
% makeplot: 1 if you want to make a plot, 0 if you don't
%
% timestampvec: a datetime vector that contains the times at which we have data (lactate or delta).
%
% tL: a datetime vector that contains only those times at which the lactate model can be plotted since it uses a moving window for the lactate signal.
%     this handles things like artifacts and letting the lactate sensor settle down.  
% 
% epoch_length: length of 1 epoch in seconds (typically 2,4,8, or 10)



% this function calculates a simple exponential model, like the one
% in Franken, et al J Neurosci 2001.  S goes up if the mouse is
% awake or in REM sleep and goes down if it is asleep.  
%
%
% If lactate is used as the signal, a moving window of length
% window_length (in hours) is used to compute LA and UA.  This is
% to account for drifting of the signal. 

%makeplot=1;  % flag 



  if nargin==10
    filename='';
  end

 
  exp_rise=exp(-dt/ti);  %compute these once to save time
  exp_fall=exp(-dt/td);

d=diff(timestampvec);   % array of time differences, so I know when an epoch has been dropped due to artifact
d=hours(d);   % convert d to hours

dL=diff(tL);
dL=hours(dL);


% CASE 1: using delta power histogram to choose upper and lower
% assymptotes for the model
if length(LA)==1

  % preallocate for speed
  S=zeros(1,size(dataset,1));

    S(1)=S0;
  

  % first run it for 24 hours, and use the ending value as the
  % starting value (like Franken et al 2000 Fig 1c)
    if size(dataset,1) >= 8640 
      iters=8640;
    else
      iters=size(dataset,1);
    end
    
  
    for i=1:iters-1                 % 8640 10-second intervals=24 hours
      if d(i) ~= dt       % if there is missing data (artifact) don't update S for the next time step
        S(i+1) = S(i);
      else
        if dataset(i,1)==0 || dataset(i,1)==2 || dataset(i,1)==4 %wake, REMS, or active wake 
          S(i+1)=UA-(UA-S(i))*exp_rise;
        elseif dataset(i,1)==1 || dataset(i,1)==3 %SWS or quiet wake
          S(i+1)=LA+(S(i)-LA)*exp_fall;
        else
        error('I found a sleepstate value that was not 0,1,2,3, or 4')
        end 
      end
    end
  

 
  % Now start the simulation over, using the last value out of S to be the new 
  % starting value of S
    if size(dataset,1) >= 8640
      S(1)=S(8640);
    else
      S(1) = S(size(dataset,1));
    end

  %start_value = S(end);

  % Testing: use MATLAB's ode45 for updating S
  % tspan=0:dt:dt*(size(dataset,1)-1);
  % issleep=(dataset(:,1)==1);
  % [T,S]=ode45(@(t,S) homeostatode(t,S,issleep,ti,td,LA,UA),tspan,start_value); 
  
    for i=1:size(dataset,1)-1
      if d(i) ~= dt       % if there is missing data (artifact) don't update S for the next time step
        S(i+1) = S(i);
      else 
        if dataset(i,1)==0 || dataset(i,1)==2 || dataset(i,1)==4 %wake, REMS, or active wake 
          S(i+1)=UA-(UA-S(i))*exp_rise;
        elseif(dataset(i,1)==1) || dataset(i,1)==3 %SWS or quiet wake
        S(i+1)=LA+(S(i)-LA)*exp_fall;
        else
          error('I found a sleepstate value that was not 0,1,2,3, or 4')
        end
      end
    end
  



% CASE 2: using lactate to choose upper and lower asymptotes. In
% this case we use a moving window, so S is not the same size at
% dataset(:,1). The moving window is 2 hours long, so an hour of
% the data at the beginning is discarded, as well as an hour of
% data at the end.  So S has length length(dataset(:,1))-720. 
% UA(1) does not correspond to the beginning of the recording, 
% but rather half of one window length into the recording.  
elseif length(LA) ~=1

%preallocate for speed
%S = zeros(1,size(dataset,1)-(window_length*(60*60/epoch_length)));
S = zeros(size(tL));

% initialize S to lactate value 
    %S(1)=dataset((window_length/2)*(60*60/epoch_length),2) 
S(1) = S0;
  
% determine a way to index through the data excluding the half-windows at the beginning and the end
l_start_index = find(timestampvec==tL(1));

 
    %for i=1:size(dataset,1)-(window_length*(60*60/epoch_length)+1)
    for i=1:length(tL)-1
      if dL(i) ~= dt       % if there is missing data (artifact) don't update S for the next time step
        S(i+1) = S(i);
      else
      %if dataset(i+(window_length/2)*(60*60/epoch_length),1)==0 || dataset(i+(window_length/2)*(60*60/epoch_length),1)==2 || dataset(i+(window_length/2)*(60*60/epoch_length),1)==4    %wake,REMS, or active wake 
      %UA_LA_index = find(timestampvec==tL(i));  % UA and LA are length of timestampvec, not tL.  
        if dataset(l_start_index+i-1,1)==0 || dataset(l_start_index+i-1,1)==2 || dataset(l_start_index+i-1,1)==4  %wake, REMS, or active wake 
          S(i+1)=UA(i)-(UA(i)-S(i))*exp_rise;
        % if S(i+1) > UA(UA_LA_index+1)
        %   S(i+1) = UA(UA_LA_index+1);
        %end
      
      %elseif(dataset(i+(window_length/2)*(60*60/epoch_length),1)==1)  || dataset(i+(window_length/2)*(60*60/epoch_length),1)==3  %SWS or quiet wake
        elseif dataset(l_start_index+i-1,1)==1  || dataset(l_start_index+i-1,1)==3  %SWS or quiet wake
       
          S(i+1)=LA(i)+(S(i)-LA(i))*exp_fall;
        % if S(i+1) < LA(UA_LA_index+1)
        %   S(i+1) = LA(UA_LA_index+1);
        %end
        else
          error('I found a sleepstate value that was not 0,1,2,3 or 4')
        end
      end
    end
  
  end %cases: delta power or lactate used for asymptotes

 
  if makeplot==1
    
    figure
    %t=0:dt:dt*(size(dataset,1)-1);
    %tS=t((window_length/2)*(60*60/epoch_length)+1:end-(window_length/2)*(60*60/epoch_length));

  
  
    if length(LA) ~= 1  
      %t=0:dt:dt*(size(dataset,1)-1);
      
      only_sleep_indices=find(dataset(:,1)==1);
      only_wake_indices=find(dataset(:,1)==0);
      only_rem_indices=find(dataset(:,1)==2);
      sleep_lactate=dataset(only_sleep_indices,2);
      wake_lactate=dataset(only_wake_indices,2);
      rem_lactate=dataset(only_rem_indices,2);

      sleep_state = dataset(:,1);
      scatter(timestampvec(only_wake_indices),wake_lactate,25,'r')
      
      hold on
      scatter(timestampvec(only_sleep_indices),sleep_lactate,25,'k')
      scatter(timestampvec(only_rem_indices),rem_lactate,25,'c')
      plot(tL,S) 
      
    else 
      plot(timestampvec,S)
    end
    title(['best fit for file ' filename])
      xlabel('TIME [h]')
      legend('wake','sleep','rem','best fit model')
      legend BOXOFF
    
      hold off
    end


    