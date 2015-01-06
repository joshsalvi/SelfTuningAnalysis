cell=5;
N1=10;
N2=13;
measuredist= [10 50 100 200 270];
ksf=120e-6;
forcein = [0 8.3 -8.3 16.7 -16.7 25 -25 33.3 -33.3 41.7 -41.7 50 -50];
for j = 1:15
    for k = 1:13
        ind2{j,k}=100+(k-1)*400+(j-1)*5200;
    end
end

for j = 1:N1
    for k = 1:N2
        for l = 1:length(measuredist)
            respsize{l}(j,k) = mean(Xd{1,cell}(ind2{j,k}+measuredist(l):ind2{j,k}+measuredist(l)+20)) - mean(Xd{1,cell}(ind2{j,k}-50:ind2{j,k}-30));
            inputsize{l}(j,k) = mean(Xo{1,cell}(ind2{j,k}+measuredist(l):ind2{j,k}+measuredist(l)+20)) - mean(Xo{1,cell}(ind2{j,k}-50:ind2{j,k}-30));
        end
    end
end

for k = 1:N2
    for l = 1:length(measuredist)
        respsizeavg{l}(k)=mean(respsize{l}(:,k));
        respsizesem{l}(k)=std(respsize{l}(:,k))/sqrt(N1);
        inputsizeavg{l}(k)=mean(inputsize{l}(:,k));
        inputsizesem{l}(k)=std(inputsize{l}(:,k))/sqrt(N1);
        forceavg{l} = (inputsizeavg{l}-respsizeavg{l})*ksf*1e12*1e-9;
        [fitcalc{l},fitgofcalc{l},fitoutputcalc{l}] = fit(sort(respsizeavg{l})',sort(forceavg{l})','poly1');
        [fitinput{l},fitgofinput{l},fitoutputinput{l}] = fit(sort(respsizeavg{l})',sort(forcein)','poly1');
    end
end

figure(cell);
set(0,'DefaultAxesColorOrder',cool(5));
for l = 1:5
    plot(sort(respsizeavg{l}),sort(forceavg{l}));hold all;
end
title('Calculated Force');

figure(cell+1);
set(0,'DefaultAxesColorOrder',cool(5));
for l = 1:5
    plot(sort(respsizeavg{l}),sort(forcein));hold all;
end
title('Input Force');

%%
cell=6;
N2=13;
measuredist= [10 50 100 200 270];
ksf=120e-6;
forcein = [0 8.3 -8.3 16.7 -16.7 25 -25 33.3 -33.3 41.7 -41.7 50 -50];

  for l = 1:5
    for k = 1:N2
  
    avgresp{l}(k) = mean(Xd{1,cell}(100+(k-1)*400+measuredist(l):100+(k-1)*400+measuredist(l)+10))-mean(Xd{1,cell}(100+(k-1)*400-70:100+(k-1)*400-20));
    avginp{l}(k) = mean(Xo{1,cell}(100+(k-1)*400+measuredist(l):100+(k-1)*400+measuredist(l)+10))-mean(Xo{1,cell}(100+(k-1)*400-70:100+(k-1)*400-20));
    avgforcecalc{l} = (avginp{l}-avgresp{l})*ksf*1e12*1e-9;
    
    end
    fitmodel='kinf*x-30*((1+exp(-0.7*(x-x0)/(4.1)))^-1)*0.7+F0';
    [fitcalc{l},fitgofcalc{l},fitoutputcalc{l}] = fit(sort(avgresp{l})',sort(avgforcecalc{l})',fitmodel);
    [fitinput{l},fitgofinput{l},fitoutputinput{l}] = fit(sort(avginp{l})',sort(forcein)',fitmodel);
  end

  figure(cell);
set(0,'DefaultAxesColorOrder',cool(5));
for l = 1:5
    plot(sort(avgresp{l}),sort(avgforcecalc{l}));hold all;
end
title('Calculated Force');
legend(num2str(fitcalc{1}.kinf),num2str(fitcalc{2}.kinf),num2str(fitcalc{3}.kinf),num2str(fitcalc{4}.kinf),num2str(fitcalc{5}.kinf))
%}
 
figure(cell+1);
set(0,'DefaultAxesColorOrder',cool(5));
for l = 1:5
    plot(sort(avgresp{l}),sort(forcein));hold all;
end
title('Input Force');
legend(num2str(fitinput{1}.kinf),num2str(fitinput{2}.kinf),num2str(fitinput{3}.kinf),num2str(fitinput{4}.kinf),num2str(fitinput{5}.kinf))
%}


%% PLOT (from RAW)
l = 3;
figure(cell);
scatter(respsizeavg{l},forceavg{l},'k.');

%% PLOT (from AVERAGED)

l = 3;
figure(cell);
scatter(avgresp{l},forceavg{l},'k.');


