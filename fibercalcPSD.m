function [ksf,lambdaxx,lambdadx,R2,CI] = fibercalcPSD(varargin)
%
% ----------------------
% [ksf,lambdaxx,lambdadx,R2,CI] = fibercalcPSD(file,minf,maxf)
% ----------------------
%
% Calculates a fiber's stiffness and drag coefficients using two methods.
%
% file: MAT file path (e.g. '/Users/username/Desktop/Extracted Data.mat');
%   -> Use importPCdata() to generate this file.
% minf,maxf: minimum and maximum frequencies over which to evaluate
%
% ksf: fiber's stiffness, from two methods
% lambdaxx,lambdadx: fiber's drag coefficients
% R2: R2-values for fits
% CI: 95% confidence intervals for fits
%
% ----------------------
% Joshua D. Salvi
% jsalvi@rockefeller.edu
%


datafile = varargin{1};

% load dat
load(datafile)

lP = length(PSDAll);

% define fits
fit1 = fittype('1/(a*(b^2+x^2))');
fit2 = fittype('a/(b^2+x^2)');

if numel(varargin) < 2
    minfreq = 20;
else
    minfreq = varargin{2};
end
if numel(varargin) < 3
    maxfreq = 2e3;
else
    maxfreq = varargin{3};
end

for j = 1:lP
    
    % convert to rad/s
    freqrad{j} = FreqAll{j}.*(2*pi);
    
    % indices of start- and end-points
    qmin{j} = findnearest(FreqAll{j},minfreq);qmin{j}=qmin{j}(1);
    qmax{j} = findnearest(FreqAll{j},maxfreq);qmax{j}=qmax{j}(1);
    
    % complete fit A, from Bormuth et al. (2014)
    disp('Fit A')
    rt1=0;rt2=0;
    clear outliers
    warning off
    fitopts=fitoptions('Method','NonlinearLeastSquares','TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',4000,'MaxIter',4000,'Robust','Off');
    [fitA{j},gofA{j}] = fit(freqrad{j}(qmin{j}:qmax{j}),PSDAll{j}(qmin{j}:qmax{j}),fit1,'Algorithm','Levenberg-Marquardt');
    Aconfint = confint(fitA{j});
    rt1 = sign(Aconfint(1,1)) + sign(Aconfint(2,1));
    rt2 = sign(Aconfint(1,2)) + sign(Aconfint(2,2));
    m = 1;
    Acoeffs{j} = coeffvalues(fitA{j});
    % remove outliers if necessary
    while rt1 == 0 || rt2 == 0 || sign(Acoeffs{j}(1)) == -1 || sign(Acoeffs{j}(2)) == -1
        m = m +1;
        fdata = feval(fitA{j},freqrad{j});
        I = abs(fdata(qmin{j}:qmax{j}) - PSDAll{j}(qmin{j}:qmax{j})) > 1.5*std(PSDAll{j}(qmin{j}:qmax{j}));
        outliers = excludedata(freqrad{j}(qmin{j}:qmax{j}),PSDAll{j}(qmin{j}:qmax{j}),'indices',I);
        fitopts=fitoptions('Method','NonlinearLeastSquares','TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',4000,'MaxIter',4000,'Robust','Off','Exclude',outliers);
        [fitA{j},gofA{j}] = fit(freqrad{j}(qmin{j}:qmax{j}),PSDAll{j}(qmin{j}:qmax{j}),fit1,fitopts);
        Aconfint = confint(fitA{j});
        rt1 = sign(Aconfint(1,1)) + sign(Aconfint(2,1));
        rt2 = sign(Aconfint(1,2)) + sign(Aconfint(2,2));
        if m == 10
            break;
        end
    end
    disp(['Completed fit A in ' num2str(m) ' iterations.']);
    
    % complete fit B, using previously reported Lorentzian fit
    disp('Fit B');
    clear outliers
    fitopts=fitoptions('Method','NonlinearLeastSquares','TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',4000,'MaxIter',4000,'Robust','Off');
    [fitB{j},gofB{j}] = fit(FreqAll{j}(qmin{j}:qmax{j}),PSDAll{j}(qmin{j}:qmax{j}),fit2,'Algorithm','Levenberg-Marquardt');
    Aconfint = confint(fitB{j});
    rt1 = sign(Aconfint(1,1)) + sign(Aconfint(2,1));
    rt2 = sign(Aconfint(1,2)) + sign(Aconfint(2,2));
    Bcoeffs{j} = coeffvalues(fitB{j});
    m=1;
    % remove outliers as necessary
    while rt1 == 0 || rt2 == 0 || sign(Bcoeffs{j}(1)) == -1 || sign(Bcoeffs{j}(2)) == -1
        m=m+1;
        fdata = feval(fitB{j},FreqAll{j});
        I = abs(fdata(qmin{j}:qmax{j}) - PSDAll{j}(qmin{j}:qmax{j})) > 1.5*std(PSDAll{j}(qmin{j}:qmax{j}));
        outliers = excludedata(FreqAll{j}(qmin{j}:qmax{j}),PSDAll{j}(qmin{j}:qmax{j}),'indices',I);
        fitopts=fitoptions('Method','NonlinearLeastSquares','TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',4000,'MaxIter',4000,'Robust','Off','Exclude',outliers);
        [fitB{j},gofB{j}] = fit(FreqAll{j}(qmin{j}:qmax{j}),PSDAll{j}(qmin{j}:qmax{j}),fit2,'Algorithm','Levenberg-Marquardt','Exclude',outliers);
        Aconfint = confint(fitB{j});
        rt1 = sign(Aconfint(1,1)) + sign(Aconfint(2,1));
        rt2 = sign(Aconfint(1,2)) + sign(Aconfint(2,2));
        if m == 10
            break;
        end
    end
    disp(['Completed fit B in ' num2str(m) ' iterations.']);
    
    % output coefficient values
    Acoeffs{j} = coeffvalues(fitA{j});
    Bcoeffs{j} = coeffvalues(fitB{j});
    
    % calculate fiber's stiffness and drag coefficients for each fit
    
    lambdasfA(j) = Acoeffs{j}(1)*2*4.114;  % pN�s/nm
    ksfA(j) = lambdasfA(j)*Acoeffs{j}(2);   %pN/nm
    lambdasfA(j) = lambdasfA(j)*1e6; % nN�s/m
    ksfA(j) = ksfA(j)*1e3; % �N/m
    
    lambdasfB(j) = 4.114/(pi^2*Bcoeffs{j}(1)); % pN�s/nm    
    ksfB(j) = 2*pi*lambdasfB(j)*Bcoeffs{j}(2); % pN/nm
    lambdasfB(j) = lambdasfB(j)*1e6; % nN�s/m
    ksfB(j) = ksfB(j)*1e3; % �N/m
    
    kf(j,:) = [ksfA(j) ksfB(j)];
    lambdasf(j,:) = [lambdasfA(j) lambdasfB(j)];
    R2(j,:) = [gofA{j}.rsquare gofB{j}.rsquare];
    CI{j} = [confint(fitA{j}) confint(fitB{j})];
end

% cleanup
lambdaxx = lambdasf.* 0.94;
lambdadx = lambdasf.*0.57;
ksf = kf.*0.97;

% plot the data
plotyn = 1;
if plotyn == 1
    close all
    figure
    setfiguredefaults();
for j = 1:lP
    subplot(1,lP,j);
    loglog(FreqAll{j},PSDAll{j},'k');hold on; plot(FreqAll{j},feval(fitB{j},FreqAll{j}),'b');
    loglog(freqrad{j},PSDAll{j},'k--');plot(freqrad{j},feval(fitA{j},freqrad{j}),'r');
    legend({'PSD(Hz)';'Fit B';'PSD(rad/s)';'Fit A'})
    title(['Fit A R^2: ' num2str(R2(j,1)) '; Fit B R^2: ' num2str(R2(j,2))]);
end
end

saveyn = 0;
if saveyn == 1
    save([datafile(1:end-4) '-fibercalcPSD.mat'],'ksf','lambdaxx','lambdadx','R2','CI','fitA','fitB','FreqAll','PSDAll');
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

set(0,'defaultFigureColormap',jet)
if exist('N')==1
    set(0,'defaultAxesColorOrder',jet(N))
else
    set(0,'defaultAxesColorOrder',[0 0 0])
end
set(0,'DefaultTextFontName','Helvetica Neue')
set(0,'DefaultTextFontUnits','Points')
set(0,'defaultlinelinewidth',2)
set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontName','Helvetica Neue')
set(0,'DefaultAxesFontSize',16)
set(0,'DefaultAxesGridLineStyle',':')
set(0,'defaultFigureColor','White')
set(0,'DefaultFigurePaperType','a4letter')
set(0,'defaultFigurePosition',[1200 100 700 700])
set(0,'defaultAxesColor',[1 1 1])
set(0,'defaultLineMarker','None')
end
