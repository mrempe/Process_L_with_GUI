function [LA,UA,h]=make_frequency_plot(datafile,window_length,signal,timestampvec,tL,epoch_length,makefig,newfig)
% USAGE: [LA,UA]=make_frequency_plot(datafile,signal)
%
% this function reads in a sleep data file with two columns: sleep state, signal (delta or lactate)
% given to me by J. Wisor 2011. 
% signal is 'lactate' or 'delta'
% 
% timestampvec: datetime vector with the epochs that correspond to all the data (delta or lactate),
%               after removing artifacts and bad data.  
% 
% tL: datetime vector with the epochs that correspond to the lactate signal, after removing 
%     the first window_length/2 and last window_length/2 due to the moving window averaging. 
%     This was set up taking into account artifacts and so on. tL is used in this function 
%     for plotting and for setting up UA and LA if the signal is lactate
% optional argument is 1 if make the frequency histogram at all, newfig is if you
% want the histogram in a new figure, if not, use the current figure. a new figure
if nargin==5 makefig=0; newfig=0; end

if nargin==6 newfig=0; end 
%
% window_length is the length (in hours) of the moving window used
% in the calculation of UA and LA.
% 

data=datafile(:,2);


% flag for plotting, if 1 you get a plot of the histogram for the
% moving window of the data, showing UA and LA as they change with
% time
animation=0;




if strcmp(signal,'delta1') || strcmp(signal,'delta2') || strcmp(signal,'EEG1') || strcmp(signal,'EEG2')
  % If using delta power like Franken 2001 does, 
  % to find LA (lower assymptote), the intersection of the REM
  % histogram and SWS histogram

  % Find all the rows where sleep state is SWS (denoted by a 1)
  %rowsleep=find(datafile(:,1)==1);
  sleepdata=data(datafile(:,1)==1);
 

  % Find all the rows where sleep state is REM (denoted by a 2)
  %rowrem=find(datafile(:,1)==2);
  remdata=data(datafile(:,1)==2);
  
  % Find all the rows where sleep state is wake (denoted by 0)
  %rowwake=find(datafile(:,1)==0);
  wakedata=data(datafile(:,1)==0);

% size(datafile)
% size(sleepdata)
% max(sleepdata)  
  xbins=linspace(0,max(sleepdata),30);
  
  % compute the histograms
  [ns,xs]=hist(sleepdata,xbins);  %sleep data
  [nr,xr]=hist(remdata,xbins);   %REM data 
  difference=ns-nr;  % need to find where this is 0
  id=find(diff(difference >= 0)); % finds crossings of the difference
                                  % vector, (keep only last one)
  max_REM_loc = xr(find(nr==max(nr)));
  max_SWS_loc = xs(find(ns==max(ns)));

  % New method to find crossing of REM and SWS histograms.
  % First find all of the crossing (some of these may be way off in the tails)
  % Then choose the one that is between the maximums of each curve
  LA = xs(id)-(difference(id).*(xs(id+1)-xs(id))./(difference(id+1)-difference(id)));
  index_crossing_between_maxREM_andmaxSWS = find(LA>max_REM_loc & LA<max_SWS_loc);
  if ~isempty(index_crossing_between_maxREM_andmaxSWS) 
    LA = LA(index_crossing_between_maxREM_andmaxSWS);
  else                        % if the two histograms don't cross between their maximums
   LA=xr(find(nr==max(nr)));  % set to max value in REMS histogram 
  end

  % if(isempty(id)) % if they don't cross
  %   loc=find(nr>0);  % find non-zero bins of REM
  %   loc2=find(ns>0); % find non-zero bins of SWS
  %   if isempty(loc) 
  %     LA=xs(loc2(1))
  %   %else LA=xr(loc(1)); % set to first non-zero bin  
  %   else 
  %     LA=xr(find(nr==max(nr)));  % set to max value in REMS histogram
  %   end
  % else
  %   id=id(end);
  %   LA = xs(id)-(difference(id)*(xs(id+1)-xs(id))/(difference(id+1)-difference(id)));
  % end
% UA (upper assymptote)
  UA=quantile(sleepdata,.9);
LA
UA


% plot the histograms
  if newfig==1 
    figure
    bar(xs,ns) 
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',[0.5 0.5 0.5],'EdgeColor','k')
    hold on
    bar(xr,nr)
    %h = findobj(gca,'Type','patch');
    %set(h,'FaceColor',[0 0.5 0.5],'EdgeColor','w')
  
  %[nw,xw]=hist(wakedata,xbins)  % wake data
  % h = findobj(gca,'Type','patch');
  % set(h,'FaceColor',[0 1 0],'EdgeColor','w')
    line([LA LA],[0 max(ns)])
    line([UA UA],[0 max(ns)])
    %plot(LA,0:max(ns)) % plot vertical line at LA
    %plot(UA,0:max(ns))
    
    hold off
    xlabel('DELTA POWER [\mu V^2]')
    ylabel('FREQUENCY')
    
  end
end


% --- Case where lactate is the signal --------
%
% For this case use a sliding window to compute UA and LA.
% Compute the frequency plot of all states (sleep, wake, REM) for
% the window centered at that time. Compute UA and LA for
% these data and use them in the next step of S.  Compute UA and LA
% again for the next window and use this in the next update of S.

if strcmp(signal,'lactate') || strcmp(signal,'Lactate')
  
 % make a histogram of the data initially
  if newfig==1 
    figure 
    xbins=linspace(0,max(data(1:(window_length)*(60*60/epoch_length)+1)),30);
    [nall,xall]=hist(data(1:(window_length)*(60*60/epoch_length)+1),xbins);
    h=bar(xall,nall);
    axis([0 19.5 0 500])
  end


  % animation stuff
  if(animation)
    figure
    xbins=linspace(0,max(data(1:(window_length)*(60*60/epoch_length)+1)),30);
    [nall,xall]=hist(data(1:(window_length)*(60*60/epoch_length)+1),xbins);
    h=bar(xall,nall)
    axis([0 19.5 0 500]) 
  end
  
  UA=zeros(1,length(tL));
  LA=zeros(1,length(tL));

  %l_start_index = find(timestampvec==tL(1))
  
  total_epochs_no_missing_data = seconds(timestampvec(end)-timestampvec(1))/epoch_length; %assuming no missing data
  epochs_in_window = (window_length*60*60)/epoch_length;


  % start_time=timestampvec(1);
  % for i=1:total_epochs_no_missing_data-epochs_in_window+1
  %   %start_time = timestampvec(i);
  %   end_time = timestampvec(i)+hours(window_length);
  %   indices_in_window = find(timestampvec>start_time & timestampvec<end_time);

  %   LA(i) = quantile(data(indices_in_window),0.01);
  %   UA(i) = quantile(data(indices_in_window),0.99);
  
  %   start_time = start_time + epoch_length;
  % end

% TRY AGAIN. UA and LA should be the same size as tL
half_window = hours(window_length/2);
start_times = tL-half_window;   % vectorized?
end_times = tL+half_window;

for i=1:length(tL)
  %start_time = tL(i) - half_window;
  %end_time   = tL(i) + half_window;
  indices_in_window = find(timestampvec>=start_times(i) & timestampvec<=end_times(i));

  LA(i) = quantile(data(indices_in_window),0.01);
  UA(i) = quantile(data(indices_in_window),0.99);
  
end


  % pause

  % for i=1:length(tL)-1
  % start_window_index = l_start_index+i-1-(epochs_in_window/2)
  % end_window_index   = l_start_index+i-1+(epochs_in_window/2)

  % LA(i) = quantile(data(start_window_index:end_window_index),0.01);
  % UA(i) = quantile(data(start_window_index:end_window_index),0.99);

  % end


  % shift=(window_length/2)*(60*60/epoch_length);
  % %for i=361:(length(data)-(60*60/epoch_length))
  % for i=shift+1:(length(data)-shift)
  %      % if using lactate, the histograms for SWS and REM will overlap,
  %   % so just use the 1st percentile for SWS
  %   % LA(i-(60*60/epoch_length))=quantile(data(i-(60*60/epoch_length):i+(60*60/epoch_length)),.01);
  %   % UA(i-(60*60/epoch_length))=quantile(data(i-(60*60/epoch_length):i+(60*60/epoch_length)),.99);
  %      LA(i-shift)=quantile(data(i-shift:i+shift),.01);
  %      UA(i-shift)=quantile(data(i-shift:i+shift),.99);
    
    
  %   if(animation)
  %     xbins=linspace(0,max(data(i-shift:i+shift)),30);
  %     [nall,xall]=hist(data(i-shift:i+shift),xbins);
  %     set(h,'XData',xall,'YData',nall)
  %     l1=line([LA(i-shift) LA(i-shift)],[0 max(nall)]);
  %     l2=line([UA(i-shift) UA(i-shift)],[0 max(nall)]);
  %     drawnow
  %     delete(l1)
  %     delete(l2)
  %   end
    
  % end
h=gcf;
end

