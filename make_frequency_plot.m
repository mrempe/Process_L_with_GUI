function [LA,UA,h]=make_frequency_plot(datafile,window_length,signal,epoch_length,makefig,newfig)
% USAGE: [LA,UA]=make_frequency_plot(datafile,signal)
%
% this function reads in a sleep data file with two columns: sleep state, signal (delta or lactate)
% given to me by J. Wisor 2011. 
% signal is 'lactate' or 'delta'
%
% optional argument is 1 if make the frequency histogram at all, newfig is if you
% want the histogram in a new figure, if not, use the current figure. a new figure
if nargin==4 makefig=0; newfig=0; end
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
  
  if(isempty(id)) % if they don't cross
    loc=find(nr>0);  % find non-zero bins of REM
    loc2=find(ns>0); % find non-zero bins of SWS
    if isempty(loc) LA=xs(loc2(1))
    else LA=xr(loc(1)); % set to first non-zero bin  
    end
  else
    id=id(end);
    LA = xs(id)-(difference(id)*(xs(id+1)-xs(id))/(difference(id+1)-difference(id)));
  end
% UA (upper assymptote)
  UA=quantile(sleepdata,.9);

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
    plot(LA,0:max(ns)) % plot vertical line at LA
    plot(UA,0:max(ns))
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
  
  
  shift=(window_length/2)*(60*60/epoch_length);
  %for i=361:(length(data)-(60*60/epoch_length))
  for i=shift+1:(length(data)-shift)
       % if using lactate, the histograms for SWS and REM will overlap,
    % so just use the 1st percentile for SWS
    % LA(i-(60*60/epoch_length))=quantile(data(i-(60*60/epoch_length):i+(60*60/epoch_length)),.01);
    % UA(i-(60*60/epoch_length))=quantile(data(i-(60*60/epoch_length):i+(60*60/epoch_length)),.99);
       LA(i-shift)=quantile(data(i-shift:i+shift),.01);
       UA(i-shift)=quantile(data(i-shift:i+shift),.99);
    
    
    if(animation)
      xbins=linspace(0,max(data(i-shift:i+shift)),30);
      [nall,xall]=hist(data(i-shift:i+shift),xbins);
      set(h,'XData',xall,'YData',nall)
      l1=line([LA(i-shift) LA(i-shift)],[0 max(nall)]);
      l2=line([UA(i-shift) UA(i-shift)],[0 max(nall)]);
      drawnow
      delete(l1)
      delete(l2)
    end
    
  end
h=gcf;
end

