function[newstate]=RescoreQuietVsActiveWakeTEST(statechars,Emgdata,QThreshold,AThreshold,FileCounter,files);

EmgWakeOnly=Emgdata(statechars==0);
newstate=statechars;
newstate(Emgdata<=quantile(EmgWakeOnly,QThreshold) & statechars==0)=3;    %Quiet wake values are relabeled as state 3.
newstate(Emgdata>=quantile(EmgWakeOnly,AThreshold) & statechars==0)=4;    %Active wake values are relabeled as state 4.

disp(['33rd percentile: ', quantile(EmgWakeOnly,QThreshold)])

% NOt finished.  Keep working on this.  
disp(['1SD left of mean of log(EMG):', exp(mean(log(EmgWakeOnly))-std(log(EmgWakeOnly)))])


if quantile(EmgWakeOnly,AThreshold)/quantile(EmgWakeOnly,QThreshold)<1.5
    warning('Poor EMG dynamics in this file:');
    display( files(FileCounter) );
end

end