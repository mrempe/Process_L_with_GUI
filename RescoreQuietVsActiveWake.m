function[newstate]=RescoreQuietVsActiveWake(statechars,Emgdata,QThreshold,AThreshold,FileCounter,files);

EmgWakeOnly=Emgdata(statechars==0);
newstate=statechars;
newstate(Emgdata<=round(quantile(EmgWakeOnly,QThreshold),4) & statechars==0)=3;    %Quiet wake values are relabeled as state 3.
newstate(Emgdata>=round(quantile(EmgWakeOnly,AThreshold),4) & statechars==0)=4;    %Active wake values are relabeled as state 4.



if quantile(EmgWakeOnly,AThreshold)/quantile(EmgWakeOnly,QThreshold)<1.5
    warning('Poor EMG dynamics in this file:');
    display( files(FileCounter) );
end

end