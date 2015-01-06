close all;
nperc=0.7;
ind = [2 8 12 15 19 20 27 35:45];

for j = ind
      if isempty(pk{1,j}) == 0 && isempty(tr{1,j}) == 0
       loopmaxtr = length(tr{1,j});loopmaxpk = length(pk{1,j});
       if pk{1,j}(1) < tr{1,j}(1) && pk{1,j}(end) > tr{1,j}(end)
           for l = 1:min([loopmaxpk loopmaxtr])-1
               if l==1
                   clear trthresh
                   trthresh = nperc*(mean([Xds{1,j}(pk{1,j}(l)) Xds{1,j}(pk{1,j}(l+1))])-Xds{1,j}(tr{1,j}(l))) + Xds{1,j}(tr{1,j}(l));
                   trtime{1,j}(l) = sum(Xds{1,j}(pk{1,j}(l):pk{1,j}(l+1))<=trthresh)/Fs;
               else
               clear trthresh pkthresh
               trthresh = nperc*(mean([Xds{1,j}(pk{1,j}(l)) Xds{1,j}(pk{1,j}(l+1))])-Xds{1,j}(tr{1,j}(l))) + Xds{1,j}(tr{1,j}(l));
               trtime{1,j}(l) = sum(Xds{1,j}(pk{1,j}(l):pk{1,j}(l+1))<=trthresh)/Fs;
               pkthresh = nperc*(mean([Xds{1,j}(tr{1,j}(l)) Xds{1,j}(tr{1,j}(l+1))])-Xds{1,j}(pk{1,j}(l))) + Xds{1,j}(pk{1,j}(l));
               pktime{1,j}(l) = sum((Xds{1,j}(tr{1,j}(l):tr{1,j}(l+1))>=pkthresh))/Fs; 
               end
               if l==loopmaxpk
               pkthresh = nperc*(mean([Xds{1,j}(tr{1,j}(l)) Xds{1,j}(tr{1,j}(l+1))])-Xds{1,j}(pk{1,j}(l))) + Xds{1,j}(pk{1,j}(l));
               pktime{1,j}(l) = sum((Xds{1,j}(tr{1,j}(l):tr{1,j}(l+1))>=pkthresh))/Fs; 
               end
           end
       elseif pk{1,j}(1) < tr{1,j}(1) && pk{1,j}(end) < tr{1,j}(end)
           for l = 1:min([loopmaxpk loopmaxtr])-1
               if l==1
                   clear trthresh
                   trthresh = nperc*(mean([Xds{1,j}(pk{1,j}(l)) Xds{1,j}(pk{1,j}(l+1))])-Xds{1,j}(tr{1,j}(l))) + Xds{1,j}(tr{1,j}(l));
                   trtime{1,j}(l) = sum(Xds{1,j}(pk{1,j}(l):pk{1,j}(l+1))<=trthresh)/Fs;
               else
               clear trthresh pkthresh
               trthresh = nperc*(mean([Xds{1,j}(pk{1,j}(l)) Xds{1,j}(pk{1,j}(l+1))])-Xds{1,j}(tr{1,j}(l))) + Xds{1,j}(tr{1,j}(l));
               trtime{1,j}(l) = sum(Xds{1,j}(pk{1,j}(l):pk{1,j}(l+1))<=trthresh)/Fs;
               pkthresh = nperc*(mean([Xds{1,j}(tr{1,j}(l)) Xds{1,j}(tr{1,j}(l+1))])-Xds{1,j}(pk{1,j}(l))) + Xds{1,j}(pk{1,j}(l));
               pktime{1,j}(l) = sum((Xds{1,j}(tr{1,j}(l):tr{1,j}(l+1))>=pkthresh))/Fs; 
               end
               if l==loopmaxpk
               trthresh = nperc*(mean([Xds{1,j}(pk{1,j}(l)) Xds{1,j}(pk{1,j}(l+1))])-Xds{1,j}(tr{1,j}(l))) + Xds{1,j}(tr{1,j}(l));
               trtime{1,j}(l) = sum(Xds{1,j}(pk{1,j}(l):pk{1,j}(l+1))<=trthresh)/Fs;
               end
           end
       elseif pk{1,j}(1) > tr{1,j}(1) && pk{1,j}(end) < tr{1,j}(end)
           for l = 1:min([loopmaxpk loopmaxtr])-1
               if l==1
                   clear trthresh
                   pkthresh = nperc*(mean([Xds{1,j}(tr{1,j}(l)) Xds{1,j}(tr{1,j}(l+1))])-Xds{1,j}(pk{1,j}(l))) + Xds{1,j}(pk{1,j}(l));
                   pktime{1,j}(l) = sum((Xds{1,j}(tr{1,j}(l):tr{1,j}(l+1))>=pkthresh))/Fs; 
               else
               clear trthresh pkthresh
               trthresh = nperc*(mean([Xds{1,j}(pk{1,j}(l)) Xds{1,j}(pk{1,j}(l+1))])-Xds{1,j}(tr{1,j}(l))) + Xds{1,j}(tr{1,j}(l));
               trtime{1,j}(l) = sum(Xds{1,j}(pk{1,j}(l):pk{1,j}(l+1))<=trthresh)/Fs;
               pkthresh = nperc*(mean([Xds{1,j}(tr{1,j}(l)) Xds{1,j}(tr{1,j}(l+1))])-Xds{1,j}(pk{1,j}(l))) + Xds{1,j}(pk{1,j}(l));
               pktime{1,j}(l) = sum((Xds{1,j}(tr{1,j}(l):tr{1,j}(l+1))>=pkthresh))/Fs;
               end
               if l==loopmaxpk
               trthresh = nperc*(mean([Xds{1,j}(pk{1,j}(l)) Xds{1,j}(pk{1,j}(l+1))])-Xds{1,j}(tr{1,j}(l))) + Xds{1,j}(tr{1,j}(l));
               trtime{1,j}(l) = sum(Xds{1,j}(pk{1,j}(l):pk{1,j}(l+1))<=trthresh)/Fs;
               end
           end
       elseif pk{1,j}(1) > tr{1,j}(1) && pk{1,j}(end) > tr{1,j}(end)
           for l = 1:min([loopmaxpk loopmaxtr])-1
               if l==1
                   clear trthresh
                   pkthresh = nperc*(mean([Xds{1,j}(tr{1,j}(l)) Xds{1,j}(tr{1,j}(l+1))])-Xds{1,j}(pk{1,j}(l))) + Xds{1,j}(pk{1,j}(l));
                   pktime{1,j}(l) = sum((Xds{1,j}(tr{1,j}(l):tr{1,j}(l+1))>=pkthresh))/Fs;
               else
               clear trthresh pkthresh
               trthresh = nperc*(mean([Xds{1,j}(pk{1,j}(l)) Xds{1,j}(pk{1,j}(l+1))])-Xds{1,j}(tr{1,j}(l))) + Xds{1,j}(tr{1,j}(l));
               trtime{1,j}(l) = sum(Xds{1,j}(pk{1,j}(l):pk{1,j}(l+1))<=trthresh)/Fs;
               pkthresh = nperc*(mean([Xds{1,j}(tr{1,j}(l)) Xds{1,j}(tr{1,j}(l+1))])-Xds{1,j}(pk{1,j}(l))) + Xds{1,j}(pk{1,j}(l));
               pktime{1,j}(l) = sum((Xds{1,j}(tr{1,j}(l):tr{1,j}(l+1))>=pkthresh))/Fs; 
               end
               if l==min([loopmaxpk loopmaxtr])-1
               trthresh = nperc*(mean([Xds{1,j}(pk{1,j}(l)) Xds{1,j}(pk{1,j}(l+1))])-Xds{1,j}(tr{1,j}(l))) + Xds{1,j}(tr{1,j}(l));
               trtime{1,j}(l) = sum(Xds{1,j}(pk{1,j}(l):pk{1,j}(l+1))<=trthresh)/Fs;
               end
           end
       end
      else
          trtime{1,j}=0;pktime{1,j}=0;
      end
end

for j = ind
    figure(1);
    hold on;scatter(j,mean(abs(pktime{j})),'k'); hold on;
    hold on;scatter(j,mean(abs(trtime{j})),'r'); hold on;
end

for j = ind
    figure(2);
    hold on;scatter(j,mean(abs(pktime{j}))/mean(abs(trtime{j})),'b');
end

