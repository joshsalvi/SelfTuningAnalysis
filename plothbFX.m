function plothbFX(filename)
%
% This function plots the force-displacement relation for a hair bundle.
% Use importVLCdata3() to save a MAT-file that can be imported for this
% function.
% 
% plothbFX(filename)
%
% where 'filename' is the path+name of your MAT file.
%
% Example:
% plothbFX('/Users/.../Extracted Data.mat')
%
% Joshua D. Salvi
% jsalvi@rockefeller.edu
%
%

clc;
setfiguredefaults
close all

% Load the data file
load(filename);

% Cell to analyze?
disp(['Raw indices: ' num2str(raw)]);
disp(['Averaged indices: ' num2str(nonraw)]);
disp('Select an index for analysis...');
ind = input('Index: ');
if isempty(intersect(ind,raw)) == 0
    disp('Selected RAW index.');
    mm = 1;
elseif isempty(intersect(ind,nonraw)) == 0
    disp('Selected AVERAGED index.');
    mm = 2;
else
    disp('Did not selected appropriate index.');
    return;
end

% Fiber properties
ksf = input('Fiber stiffness (µN/m): ');
kv = input('Virtual stiffness (µN/m): ');


if mm == 2
    
    % Time vector
    dt = 1/Fs;
    tvec = 0:dt:length(Xd{1,ind})*dt-dt;
    
% Plot the traces
hrt(1)=figure(1);
for j = 1:length(Xd)
    lengthX(j) = isempty(Xd{j,ind})==0;
end
lengthX = sum(lengthX);
for j = 1:lengthX
    plot(Xd{j,ind},'k');hold on;
end
set(1,'WindowStyle','docked')
grid on;

% Choose range for averaging
tr=1;
while tr==1
    disp('Analysis parameters. Select the range over which you will average...');
    range1 = input('Pre-pulse start: ');
    range2 = input('Pre-pulse end: ');
    range3 = input('Pulse start: ');
    range4 = input('Pulse end: ');
    close;
    for j = 1:lengthX
        Xd{j,ind} = Xd{j,ind}-mean(Xd{j,ind}(range1:range2));
        xdmin(j) = min(Xd{j,ind});xdmax(j) = max(Xd{j,ind});
        plot(tvec,Xd{j,ind},'k');
        hold on;plot(tvec(range1:range2),Xd{j,ind}(range1:range2),'g');
        plot(tvec(range3:range4),Xd{j,ind}(range3:range4),'r');
    end
    title('Traces with Selected Averaging Range');
    grid on;axis([tvec(1) tvec(end) 1.1*min(xdmin) 1.1*max(xdmax)]);

% Plot the traces to be averaged
close all;
hrt(1)=figure(1); subplot(2,2,1);
set(1,'WindowStyle','docked')
setfiguredefaults(length(Xd))
for j = 1:lengthX
    Xd{j,ind} = Xd{j,ind}-mean(Xd{j,ind}(range1:range2));
    xdmin(j) = min(Xd{j,ind});xdmax(j) = max(Xd{j,ind});
    plot(tvec,Xd{j,ind},'k');
    hold on;plot(tvec(range1:range2),Xd{j,ind}(range1:range2),'g');
    plot(tvec(range3:range4),Xd{j,ind}(range3:range4),'r');
    xdpremean(j) = mean(Xd{j,ind}(range1:range2));
    xdpresem(j) = std(Xd{j,ind}(range1:range2))/sqrt(length(Xd{j,ind}(range1:range2)));
    xdpulsemean(j) = mean(Xd{j,ind}(range3:range4));
    xdpulsesem(j) = std(Xd{j,ind}(range3:range4))/sqrt(length(Xd{j,ind}(range3:range4)));
    xddiffmean(j) = xdpulsemean(j) - xdpremean(j);
    xddiffsem(j) = xdpulsesem(j) + xdpresem(j);
end
title('Photodiode output');xlabel('Time (s)');ylabel('Position (nm)');
grid on;axis([tvec(1) tvec(end) 1.1*min(xdmin) 1.1*max(xdmax)]);
hrt(1)=figure(1); subplot(2,2,3);
set(1,'WindowStyle','docked')
setfiguredefaults(length(Xd))
for j = 1:lengthX
    Xo{j,ind} = Xo{j,ind}-mean(Xo{j,ind}(range1:range2));
    xdmin(j) = min(Xo{j,ind});xdmax(j) = max(Xo{j,ind});
    plot(tvec,Xo{j,ind},'k');
    hold on;plot(tvec(range1:range2),Xo{j,ind}(range1:range2),'g');
    plot(tvec(range3:range4),Xo{j,ind}(range3:range4),'r');
    xopremean(j) = mean(Xo{j,ind}(range1:range2));
    xopresem(j) = std(Xo{j,ind}(range1:range2))/sqrt(length(Xo{j,ind}(range1:range2)));
    xopulsemean(j) = mean(Xo{j,ind}(range3:range4));
    xopulsesem(j) = std(Xo{j,ind}(range3:range4))/sqrt(length(Xo{j,ind}(range3:range4)));
    xodiffmean(j) = xopulsemean(j) - xopremean(j);
    xodiffsem(j) = xdpulsesem(j) + xdpresem(j);
    Ffibermean(j) = kv*1e-6*(xodiffmean(j)*1e-9 - xddiffmean(j)*1e-9);Ffibermean(j) = Ffibermean(j)*1e12;
    Ffibersem(j) = kv*1e-6*(xodiffsem(j)*1e-9 + xddiffsem(j)*1e-9);Ffibersem(j) = Ffibersem(j)*1e12;
end
title('Fiber displacement');xlabel('Time (s)');ylabel('Position (nm)');
grid on;axis([tvec(1) tvec(end) 1.1*min(xdmin) 1.1*max(xdmax)]);


forces = input('Input forces manually (1=yes)? ');
if forces == 1
    forcerange = 0;
    while length(forcerange) ~= length(xddiffmean)
        disp(['Please input ' num2str(length(xddiffmean)) ' forces.']);
        forcerange = input('Forces (pN): ');
        if length(forcerange) ~= length(xddiffmean)
            disp('Incorrect number of forces.')
        else
            %return;
        end
    end
    subplot(2,2,2);
    plot(xddiffmean,forcerange,'r:','LineWidth',2); hold on;
    herrorbar(xddiffmean,forcerange,xddiffsem,'b.');
    title('Force versus displacement');xlabel('Bundle displacement (nm)');ylabel('Force (pN)');
    grid on;axis([1.1*min(xddiffmean) 1.1*max(xddiffmean) 1.1*min(forcerange) 1.1*max(forcerange)]);
    subplot(2,2,4);
    diffF = gradient(forcerange*1e-12);
    diffX = gradient(xddiffmean*1e-9);
    stiffFx = diffF./diffX * 1e6;
    plot(xddiffmean,stiffFx,'LineWidth',2);
    hold on;plot(linspace(-100,100,20),zeros(1,20),'g:');
    stiffmedrange=floor(length(stiffFx)/4):floor(3*length(stiffFx)/4);stdstiff=std(stiffFx(stiffmedrange));
    grid on;grid on;axis([1.1*min(xddiffmean) 1.1*max(xddiffmean) -stdstiff stdstiff]);
    title('Stiffness versus displacement');xlabel('Bundle displacement (nm)');ylabel('Stiffness (µN/m)');
    
    % Fit the data and plot the open probability
    numevals = input('Number of fitting iterations: ');mm=1;
    if numevals<1
        numevals=1;
    end 
    disp('Fitting...');fitobject1 = fittype('kinf.*x-N*z*(1./(1+exp(-z.*(x-x0)./4.11)))+F0');
    fitopts=fitoptions('Method','NonlinearLeastSquares','TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',6000,'MaxIter',4000,'Robust','Off');
    [fitobj gof]=fit(xddiffmean',forcerange',fitobject1,fitopts);
    figure(6);set(6,'WindowStyle','docked');scatter(xddiffmean,forcerange,'bo');hold on;plot(xddiffmean,feval(fitobj,xddiffmean));hold all;
    disp(['R2 = ' num2str(gof.rsquare)]);
    while mm < numevals
        fdata=feval(fitobj,xddiffmean);I = abs(fdata - forcerange') > 1.5*std(forcerange');
        disp('Removing outliers...');outliers = excludedata(xddiffmean,forcerange,'indices',I);
        fitopts=fitoptions('Method','NonlinearLeastSquares','TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',6000,'MaxIter',4000,'Robust','Off','Exclude',outliers);
        [fitobj gof]=fit(xddiffmean',forcerange',fitobject1,fitopts);
        plot(xddiffmean,feval(fitobj,xddiffmean));hold all;
        disp(['R2 = ' num2str(gof.rsquare)]);
        mm=mm+1;
    end
    close(6);fitobj.z=abs(fitobj.z);Po=1./(1+exp(-fitobj.z.*(xddiffmean-fitobj.x0)./4.11));disp(['z = ' num2str(fitobj.z)]);disp(['sqrt(4*k*T*Kinf/N) = ' num2str(sqrt(4*4.11*fitobj.kinf/fitobj.N))]);disp(['Kinf = ' num2str(fitobj.kinf)]);
    Ff=fitobj.kinf.*xddiffmean-fitobj.N*fitobj.z.*Po+fitobj.F0;
    hrt(4)=figure(4);set(4,'WindowStyle','docked');
    plot(xddiffmean,Po,'LineWidth',2,'Color','b');
    grid on;axis([1.1*min(xddiffmean) 1.1*max(xddiffmean) 0 1.5]);
    xlabel('Bundle position (nm)');ylabel('Open probability (from fit)');
    axes('Position',[.7 .7 .2 .2]); box on;
    scatter(xddiffmean,forcerange,'bo');hold on;
    plot(xddiffmean,Ff,'LineWidth',2,'Color','r');
    grid on;axis([1.1*min(xddiffmean) 1.1*max(xddiffmean) -2*std(forcerange) 2*std(forcerange)]);
    set(gca,'xticklabel','','yticklabel','');hold on;
    xlabel('X (nm)');ylabel('F (pN)');disp('Fitting complete.');
else
    subplot(2,2,2);
    plot(xddiffmean,Ffibermean,'r:','LineWidth',2); hold on;
    errorbar(xddiffmean,Ffibermean,Ffibersem,'b.');
    hold on;herrorbar(xddiffmean,Ffibermean,xddiffsem,'b.');
    title('Force versus displacement');xlabel('Bundle displacement (nm)');ylabel('Force (pN)');
    grid on;axis([1.1*min(xddiffmean) 1.1*max(xddiffmean) -2*std(Ffibermean) 2*std(Ffibermean)]);
    subplot(2,2,4);
    diffF = gradient(Ffibermean*1e-12);
    diffX = gradient(xddiffmean*1e-9);
    stiffFx = diffF./diffX * 1e6;
    plot(xddiffmean,stiffFx,'LineWidth',2);
    hold on;plot(linspace(-100,100,20),zeros(1,20),'g--');
    stiffmedrange=floor(length(stiffFx)/4):floor(3*length(stiffFx)/4);stdstiff=std(stiffFx(stiffmedrange));
    grid on;axis([1.1*min(xddiffmean) 1.1*max(xddiffmean) -stdstiff stdstiff]);
    title('Stiffness versus displacement');xlabel('Bundle displacement (nm)');ylabel('Stiffness (µN/m)');
    
    % Fit the data and plot the open probability
    numevals = input('Number of fitting iterations: ');mm=1;
    if numevals<1
        numevals=1;
    end 
    disp('Fitting...');fitobject1 = fittype('kinf.*x-N*z*(1./(1+exp(-z.*(x-x0)./4.11)))+F0');
    fitopts=fitoptions('Method','NonlinearLeastSquares','TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',6000,'MaxIter',4000,'Robust','Off');
    [fitobj gof]=fit(xddiffmean',Ffibermean',fitobject1,fitopts);
    figure(6);set(6,'WindowStyle','docked');scatter(xddiffmean,Ffibermean,'bo');hold on;plot(xddiffmean,feval(fitobj,xddiffmean));hold all;
    disp(['R2 = ' num2str(gof.rsquare)]);
    while mm < numevals
        fdata=feval(fitobj,xddiffmean);I = abs(fdata - Ffibermean') > 1.5*std(Ffibermean');
        disp('Removing outliers...');outliers = excludedata(xddiffmean,Ffibermean,'indices',I);
        fitopts=fitoptions('Method','NonlinearLeastSquares','TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',6000,'MaxIter',4000,'Robust','Off','Exclude',outliers);
        [fitobj gof]=fit(xddiffmean',Ffibermean',fitobject1,fitopts);
        plot(xddiffmean,feval(fitobj,xddiffmean));hold all;
        disp(['R2 = ' num2str(gof.rsquare)]);
        mm=mm+1;
    end
    close(6);fitobj.z=abs(fitobj.z);Po=1./(1+exp(-fitobj.z.*(xddiffmean-fitobj.x0)./4.11));disp(['z = ' num2str(fitobj.z)]);disp(['sqrt(4*k*T*Kinf/N) = ' num2str(sqrt(4*4.11*fitobj.kinf/fitobj.N))]);disp(['Kinf = ' num2str(fitobj.kinf)]);
    Ff=fitobj.kinf.*xddiffmean-fitobj.N*fitobj.z.*Po+fitobj.F0;
    hrt(4)=figure(4);set(4,'WindowStyle','docked');
    plot(xddiffmean,Po,'LineWidth',2,'Color','b');
    grid on;axis([1.1*min(xddiffmean) 1.1*max(xddiffmean) 0 1.5]);
    xlabel('Bundle position (nm)');ylabel('Open probability (from fit)');
    axes('Position',[.7 .7 .2 .2]); box on;
    scatter(xddiffmean,Ffibermean,'bo');hold on;
    plot(xddiffmean,Ff,'LineWidth',2,'Color','r');
    grid on;axis([1.1*min(xddiffmean) 1.1*max(xddiffmean) -2*std(Ffibermean) 2*std(Ffibermean)]);
    set(gca,'xticklabel','','yticklabel','');hold on;
    xlabel('X (nm)');ylabel('F (pN)');disp('Fitting complete.');
end
    tr = input('Repeat (1=yes)? ');

end

elseif mm == 1
    
    % Time vector
    clear dt tvec
    dt = 1/Fs;
    tvec = 0:dt:length(Xd_pulse{1,ind})*dt-dt;
    
% Plot the traces
for j = 1:length(Xd_pulse)
    lengthX(j) = isempty(Xd_pulse{j,ind})==0;
end
lengthX = sum(lengthX);

hrt(1)=figure(1);
for j = 1:1
    %subplot(ceil(sqrt(lengthX)),ceil(sqrt(lengthX)),j);
    sizeX = size(Xd_pulse{j,ind});
    for i = 1:sizeX(2)
        plot(Xd_pulse{j,ind}(:,i)-mean(Xd_pulse{j,ind}(1:pre(ind),i)),'k');hold all;grid on;axis tight;
    end
    set(1,'WindowStyle','docked')
    grid on;
end


% Choose range for averaging
tr=1;
while tr==1
    disp('Analysis parameters. Select the range over which you will average...');
    range1 = input('Pre-pulse start: ');
    range2 = input('Pre-pulse end: ');
    range3 = input('Pulse start: ');
    range4 = input('Pulse end: ');
    close;
    hrt(1)=figure(1);
    for j = 1:lengthX
        subplot(ceil(sqrt(lengthX)),ceil(sqrt(lengthX)),j);
        sizeX = size(Xd_pulse{j,ind});
        for i = 1:sizeX(2)
            Xd_pulse{j,ind}(:,i) = Xd_pulse{j,ind}(:,i)-mean(Xd_pulse{j,ind}(range1:range2,i));
            plot(tvec,Xd_pulse{j,ind}(:,i),'k');hold all;grid on;axis tight;
            plot(tvec(range1:range2),Xd_pulse{j,ind}(range1:range2,i),'g');
            plot(tvec(range3:range4),Xd_pulse{j,ind}(range3:range4,i),'r');
            xlabel('Time (s)');ylabel('Position (nm)');title('Bundle displacement');
            xdpremean(j,i) = mean(Xd_pulse{j,ind}(range1:range2,i));
            xdpresem(j,i) = std(Xd_pulse{j,ind}(range1:range2,i))/length(Xd_pulse{j,ind}(range1:range2,i));
            xdpulsemean(j,i) = mean(Xd_pulse{j,ind}(range3:range4,i));
            xdpulsesem(j,i) = std(Xd_pulse{j,ind}(range3:range4,i))/length(Xd_pulse{j,ind}(range3:range4,i));
            xddiffmean(j,i) = xdpulsemean(j,i) - xdpremean(j,i);
            xddiffsem(j,i) = xdpulsesem(j,i) + xdpresem(j,i);
            xddiffmeanavg(i) = mean(xddiffmean(:,i));xddiffsemerr(i) = sum(xddiffsem(:,i));
        end
    end
    hrt(2)=figure(2);
    for j = 1:lengthX
        subplot(ceil(sqrt(lengthX)),ceil(sqrt(lengthX)),j);
        sizeX = size(Xo_pulse{j,ind});
        for i = 1:sizeX(2)
            Xo_pulse{j,ind}(:,i) = Xo_pulse{j,ind}(:,i)-mean(Xo_pulse{j,ind}(range1:range2,i));
            plot(tvec,Xo_pulse{j,ind}(:,i),'k');hold all;grid on;axis tight;
            plot(tvec(range1:range2),Xo_pulse{j,ind}(range1:range2,i),'g');
            plot(tvec(range3:range4),Xo_pulse{j,ind}(range3:range4,i),'r');
            xlabel('Time (s)');ylabel('Position (nm)');title('Fiber displacement');
            xopremean(j,i) = mean(Xo_pulse{j,ind}(range1:range2,i));
            xopresem(j,i) = std(Xo_pulse{j,ind}(range1:range2,i))/length(Xo_pulse{j,ind}(range1:range2,i));
            xopulsemean(j,i) = mean(Xo_pulse{j,ind}(range3:range4,i));
            xopulsesem(j,i) = std(Xo_pulse{j,ind}(range3:range4,i))/length(Xo_pulse{j,ind}(range3:range4,i));
            xodiffmean(j,i) = xopulsemean(j,i) - xopremean(j,i);
            xodiffsem(j,i) = xopulsesem(j,i) + xopresem(j,i);
            xodiffmeanavg(i) = mean(xodiffmean(:,i));xodiffsemerr(i) = sum(xodiffsem(:,i));
            Ffibermeanavg(i) = kv*1e-6*(xodiffmeanavg(i)*1e-9 - xddiffmeanavg(i)*1e-9)*1e12;
            Ffibermeanerr(i) = kv*1e-6*(xodiffsemerr(i)*1e-9 + xddiffsemerr(i)*1e-9)*1e12;
        end
    end
    
    
    
    forces = input('Input forces manually (1=yes)? ');
    
    hrt(3)=figure(3);
    set(3,'WindowStyle','docked')
    j=1;
    subplot(2,2,1);
    for i = 1:sizeX(2)
            Xd_pulse{j,ind}(:,i) = Xd_pulse{j,ind}(:,i)-mean(Xd_pulse{j,ind}(range1:range2,i));
            plot(tvec,Xd_pulse{j,ind}(:,i),'k');hold all;grid on;axis tight;
            plot(tvec(range1:range2),Xd_pulse{j,ind}(range1:range2,i),'g');
            plot(tvec(range3:range4),Xd_pulse{j,ind}(range3:range4,i),'r');
            xlabel('Time (s)');ylabel('Position (nm)');title('Bundle displacement');
    end
    subplot(2,2,3);
    for i = 1:sizeX(2)
            Xo_pulse{j,ind}(:,i) = Xo_pulse{j,ind}(:,i)-mean(Xo_pulse{j,ind}(range1:range2,i));
            plot(tvec,Xo_pulse{j,ind}(:,i),'k');hold all;grid on;axis tight;
            plot(tvec(range1:range2),Xo_pulse{j,ind}(range1:range2,i),'g');
            plot(tvec(range3:range4),Xo_pulse{j,ind}(range3:range4,i),'r');
            xlabel('Time (s)');ylabel('Position (nm)');title('Fiber displacement');
    end
    
    
    
if forces == 1
    forcerange = 0;
    while length(forcerange) ~= length(xddiffmeanavg)
        disp(['Please input ' num2str(length(xddiffmeanavg)) ' forces.']);
        forcerange = input('Forces (pN): ');
        if length(forcerange) ~= length(xddiffmeanavg)
            disp('Incorrect number of forces.')
        end
    end
    subplot(2,2,2);
    plot(xddiffmeanavg,forcerange,'r:','LineWidth',2); hold on;
    herrorbar(xddiffmeanavg,forcerange,xddiffsemerr,'b.');
    title('Force versus displacement');xlabel('Bundle displacement (nm)');ylabel('Force (pN)');
    grid on;axis([1.1*min(xddiffmeanavg) 1.1*max(xddiffmeanavg) 1.1*min(forcerange) 1.1*max(forcerange)]);
    subplot(2,2,4);
    diffF = gradient(forcerange*1e-12);
    diffX = gradient(xddiffmeanavg*1e-9);
    stiffFx = diffF./diffX * 1e6;
    plot(xddiffmeanavg,stiffFx,'LineWidth',2);
    hold on;plot(linspace(-100,100,20),zeros(1,20),'g:');
    stiffmedrange=floor(length(stiffFx)/4):floor(3*length(stiffFx)/4);stdstiff=std(stiffFx(stiffmedrange));
    grid on;axis([1.1*min(xddiffmeanavg) 1.1*max(xddiffmeanavg) -stdstiff stdstiff]);
    title('Stiffness versus displacement');xlabel('Bundle displacement (nm)');ylabel('Stiffness (µN/m)');
    
    
    % Fit the data and plot the open probability
    numevals = input('Number of fitting iterations: ');mm=1;
    if numevals<1
        numevals=1;
    end 
    disp('Fitting...');fitobject1 = fittype('kinf.*x-N*z*(1./(1+exp(-z.*(x-x0)./4.11)))+F0');
    fitopts=fitoptions('Method','NonlinearLeastSquares','TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',6000,'MaxIter',4000,'Robust','Off');
    [fitobj gof]=fit(xddiffmeanavg',forcerange',fitobject1,fitopts);
    figure(6);set(6,'WindowStyle','docked');scatter(xddiffmeanavg,forcerange,'bo');hold on;plot(xddiffmeanavg,feval(fitobj,xddiffmeanavg));hold all;
    disp(['R2 = ' num2str(gof.rsquare)]);
    while mm < numevals
        fdata=feval(fitobj,xddiffmeanavg);I = abs(fdata - forcerange') > 1.5*std(forcerange');
        disp('Removing outliers...');outliers = excludedata(xddiffmeanavg,forcerange,'indices',I);
        fitopts=fitoptions('Method','NonlinearLeastSquares','TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',6000,'MaxIter',4000,'Robust','Off','Exclude',outliers);
        [fitobj gof]=fit(xddiffmeanavg',forcerange',fitobject1,fitopts);
        plot(xddiffmeanavg,feval(fitobj,xddiffmeanavg));hold all;
        disp(['R2 = ' num2str(gof.rsquare)]);
        mm=mm+1;
    end
    close(6);fitobj.z=abs(fitobj.z);Po=1./(1+exp(-fitobj.z.*(xddiffmeanavg-fitobj.x0)./4.11));disp(['z = ' num2str(fitobj.z)]);disp(['sqrt(4*k*T*Kinf/N) = ' num2str(sqrt(4*4.11*fitobj.kinf/fitobj.N))]);disp(['Kinf = ' num2str(fitobj.kinf)]);
    Ff=fitobj.kinf.*xddiffmeanavg-fitobj.N*fitobj.z.*Po+fitobj.F0;
    hrt(4)=figure(4);set(4,'WindowStyle','docked');
    plot(xddiffmeanavg,Po,'LineWidth',2,'Color','b');
    grid on;axis([1.1*min(xddiffmeanavg) 1.1*max(xddiffmeanavg) 0 1.5]);
    xlabel('Bundle position (nm)');ylabel('Open probability (from fit)');
    axes('Position',[.7 .7 .2 .2]); box on;
    scatter(xddiffmeanavg,forcerange,'bo');hold on;
    plot(xddiffmeanavg,Ff,'LineWidth',2,'Color','r');
    grid on;axis([1.1*min(xddiffmeanavg) 1.1*max(xddiffmeanavg) -2*std(forcerange) 2*std(forcerange)]);
    set(gca,'xticklabel','','yticklabel','');hold on;
    xlabel('X (nm)');ylabel('F (pN)');disp('Fitting complete.');
    
else
    subplot(2,2,2);
    plot(xddiffmeanavg,Ffibermeanavg,'r:','LineWidth',2); hold on;
    errorbar(xddiffmeanavg,Ffibermeanavg,Ffibermeanerr,'b.');
    hold on;herrorbar(xddiffmeanavg,Ffibermeanavg,xddiffsemerr,'b.');
    title('Force versus displacement');xlabel('Bundle displacement (nm)');ylabel('Force (pN)');
    grid on;axis([1.1*min(xddiffmeanavg) 1.1*max(xddiffmeanavg) -2*std(Ffibermeanavg) 2*std(Ffibermeanavg)]);
    subplot(2,2,4);
    diffF = gradient(Ffibermeanavg*1e-12);
    diffX = gradient(xddiffmeanavg*1e-9);
    stiffFx = diffF./diffX * 1e6;
    plot(xddiffmeanavg,stiffFx,'LineWidth',2);
    hold on;plot(linspace(-100,100,20),zeros(1,20),'g--');
    stiffmedrange=floor(length(stiffFx)/4):floor(3*length(stiffFx)/4);stdstiff=std(stiffFx(stiffmedrange));
    grid on;axis([1.1*min(xddiffmeanavg) 1.1*max(xddiffmeanavg) -stdstiff stdstiff]);
    title('Stiffness versus displacement');xlabel('Bundle displacement (nm)');ylabel('Stiffness (µN/m)');
    
    
    % Fit the data and plot the open probability
    numevals = input('Number of fitting iterations: ');mm=1;
    if numevals<1
        numevals=1;
    end 
    disp('Fitting...');fitobject1 = fittype('kinf.*x-N*z*(1./(1+exp(-z.*(x-x0)./4.11)))+F0');
    fitopts=fitoptions('Method','NonlinearLeastSquares','TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',6000,'MaxIter',4000,'Robust','Off');
    [fitobj gof]=fit(xddiffmeanavg',Ffibermeanavg',fitobject1,fitopts);
    figure(6);set(6,'WindowStyle','docked');scatter(xddiffmeanavg,Ffibermeanavg,'bo');hold on;plot(xddiffmeanavg,feval(fitobj,xddiffmeanavg));hold all;
    disp(['R2 = ' num2str(gof.rsquare)]);
    while mm < numevals
        fdata=feval(fitobj,xddiffmeanavg);I = abs(fdata - Ffibermeanavg') > 1.5*std(Ffibermeanavg');
        disp('Removing outliers...');outliers = excludedata(xddiffmeanavg,Ffibermeanavg,'indices',I);
        fitopts=fitoptions('Method','NonlinearLeastSquares','TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',6000,'MaxIter',4000,'Robust','Off','Exclude',outliers);
        [fitobj gof]=fit(xddiffmeanavg',Ffibermeanavg',fitobject1,fitopts);
        plot(xddiffmeanavg,feval(fitobj,xddiffmeanavg));hold all;
        disp(['R2 = ' num2str(gof.rsquare)]);
        mm=mm+1;
    end
    close(6);fitobj.z=abs(fitobj.z);Po=1./(1+exp(-fitobj.z.*(xddiffmeanavg-fitobj.x0)./4.11));disp(['z = ' num2str(fitobj.z)]);disp(['sqrt(4*k*T*Kinf/N) = ' num2str(sqrt(4*4.11*fitobj.kinf/fitobj.N))]);disp(['Kinf = ' num2str(fitobj.kinf)]);
    Ff=fitobj.kinf.*xddiffmeanavg-fitobj.N*fitobj.z.*Po+fitobj.F0;
    hrt(4)=figure(4);set(4,'WindowStyle','docked');
    plot(xddiffmeanavg,Po,'LineWidth',2,'Color','b');
    grid on;axis([1.1*min(xddiffmeanavg) 1.1*max(xddiffmeanavg) 0 1.5]);
    xlabel('Bundle position (nm)');ylabel('Open probability (from fit)');
    axes('Position',[.7 .7 .2 .2]); box on;
    scatter(xddiffmeanavg,Ffibermeanavg,'bo');hold on;
    plot(xddiffmeanavg,Ff,'LineWidth',2,'Color','r');
    grid on;axis([1.1*min(xddiffmeanavg) 1.1*max(xddiffmeanavg) -2*std(Ffibermeanavg) 2*std(Ffibermeanavg)]);
    set(gca,'xticklabel','','yticklabel','');hold on;
    xlabel('X (nm)');ylabel('F (pN)');disp('Fitting complete.');

    
end
    
    
    tr = input('Repeat (1=yes)? ');

end
    
end

saveyn = input('Save the figures? (1=yes): ');
hrt(hrt==0)=[];
if saveyn==1
    disp('Saving...');
    filename2 = filename(1:end-18);
    savefig(hrt,sprintf('%s%s%s%s',filename2,'FX-Figures-ind',num2str(ind),'.fig'));
    disp('Finished.');
else
    disp('Not saved.');
end
saveyn = input('Save the data? (1=yes): ');
if saveyn==1
    disp('Saving...');
    filename2 = filename(1:end-18);
    save([filename2 'FXdata-ind' num2str(ind) '.mat']);
    disp('Finished.');
else
    disp('Not saved.');
end
end

function setfiguredefaults(N)
% Set the figure defaults, with axes order for N axes
%
% setfiguredefaults(N,dockyn)
%
% where N is the number of axes for your colormap
% If N is not included, the default line color will be black.
%
%
% Joshua Salvi
% jsalvi@rockefeller.edu
%
%
% NOTE: Requires the 'ametrine' colormap
%

set(0,'defaultFigureColormap',jet)
if exist('N')==1
    set(0,'defaultAxesColorOrder',ametrine(N))
else
    set(0,'defaultAxesColorOrder',[0 0 0])
end
set(0,'DefaultTextFontName','Times')
set(0,'DefaultTextFontUnits','Points')
set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontName','Times')
set(0,'DefaultAxesFontSize',16)
set(0,'DefaultAxesGridLineStyle',':')
set(0,'defaultFigureColor','White')
set(0,'DefaultFigurePaperType','a4letter')
set(0,'defaultFigurePosition',[1200 100 700 700])
set(0,'defaultAxesColor',[1 1 1])
set(0,'defaultLineMarker','None')
end
