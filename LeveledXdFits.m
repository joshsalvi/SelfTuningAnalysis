function [fitparams, gof]=  LeveledXdFits(timeseries,window,PulseStart,PulseEnd)

LevelXd = zeros(1,length(timeseries)-20);
TimeLen = length(timeseries);

LevelXd(1:PulseStart)=timeseries(1:PulseStart)'-mean(timeseries(1:PulseStart));

j=0;
startval=PulseStart+10;
endval=PulseStart+150;
while j<1
        
    [fitparams,gof]=fit([0:endval-startval]',timeseries(startval:endval)-timeseries(startval),'exp2')
    coeffs=coeffvalues(fitparams);
    confints=confint(fitparams);

    figure
    plot(timeseries);
    hold on
    plot([startval-20:endval],fitparams(0-20:endval-startval)+timeseries(startval),'r');
    
    if input('Keep? Y or N (enter as string) ')=='Y'
        j=1;
%        LevelXd(PulseStart+1:endval-9)=timeseries(startval:endval)-(fitparams(0:endval-startval)+timeseries(startval));
    else
        startval=input('Start value?');
        endval=input('End value?');
    end
    
end


% MovMean3 = zeros(1,TimeLen-(PulseEnd+150));
% endnum3 = (TimeLen-(PulseEnd+150))-ceil(window/2);
% 
% for k = 1:endnum3
%     MovMean3(k) = mean(timeseries(k+PulseEnd+150-ceil(window/2):k+PulseEnd+150+ceil(window/2)));
% end
% %MovMean2(1:1+ceil(window/2))=MovMean2(1+ceil(window/2));
% MovMean3(endnum3:length(MovMean3))=MovMean3(endnum3);


end
