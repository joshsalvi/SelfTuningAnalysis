function [ksf,lambdaxx,lambdadx,kf,lambdasf,R2,CI,fitA,fitB,iterations,tel] = fibercalcPSDLV(FreqAll,PSDAll,minfreq,maxfreq,maxiter,freqguess)
%
% ----------------------
% [ksf,lambdaxx,lambdadx,kf,lambdasf,R2,CI,fitA,fitB,iterations] = fibercalcPSDLV(FreqAll,PSDAll,minf,maxf,maxiter)
% ----------------------
%
% Calculates a fiber's stiffness and drag coefficients using two methods.
% Fit A from Bormuth et al. 2014 PNAS
% Fit B from J. Howard, used in Salvi et al 2015 PNAS
%
% INPUTS:
%   -> FreqAll:  1D array of frequencies (Hz)
%   -> PSDAll:   1D array of the fiber's PSD (nm^2/Hz)
%   -> minf,maxf: minimum/maximum frequencies over which to fit
%   -> maxiter: maximum number of iterations per fit (number of times to
%   exclude outliers)
%   -> freqguess: guess of initial frequency (in Hz)
%
% OUTPUTS:
%   -> ksf: fiber's stiffness (0.97*kf) (1x2,fits A/B) [煮/m]
%   -> lambdaxx: fiber's drag coefficient owing to motion at the tip (1x2)
%   [nN新/m]
%   -> lambdadx: fiber's drag coefficient owing to motion at the base (1x2)
%   [nN新/m]
%   -> kf:  fiber's stiffness directly from fit [煮/m]
%   -> lambdasf: fiber's drag directly from fit [nN新/m]
%   -> R2:  R-squared for each fit (1x2)
%   -> CI:  95% confidence intervals for each model parameter in each fit
%   -> fitA/fitB:   fitted models for fits A and B
%   -> iterations: number of iterations for each of the two fits (1x2)
%   -> tel: total time for fitting algorithm to complete
%
% ----------------------
% Joshua D. Salvi
% jsalvi@rockefeller.edu
% ----------------------
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Initialization               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fit1 = fittype('1/(4*pi^2*a*(b^2+x^2))');                   % Bormuth et al. 2014, PNAS (fit A)
fit2 = fittype('a/((b^2+x^2))');                            % J. Howard, 2001 (fit B)
qmin = findnearest(FreqAll,minfreq);qmin=qmin(1);           % starting index
qmax = findnearest(FreqAll,maxfreq);qmax=qmax(1);           % ending index
qguess = findnearest(FreqAll,freqguess);qguess=qguess(1);   % starting frequency index
startptA = 1/(4*pi^2*(2*freqguess^2)*PSDAll(qguess));       % starting a parameter (fitA)
startptB = 2*freqguess^2*PSDAll(qguess);                    % starting a parameter (fitB)
tic;                                                        % start clock


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Fit A, from Bormuth et al. (2014)     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rt1=0;rt2=0;clear outliers; warning off;            % initialize
fitopts = fitoptions('Algorithm','Levenberg-Marquardt','Method','NonlinearLeastSquares','TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',4000,'MaxIter',4000,'Robust','Off','StartPoint',[startptA freqguess]);
[fitA,gofA] = fit(FreqAll(qmin:qmax),PSDAll(qmin:qmax),fit1,fitopts);
Aconfint = confint(fitA);                           % check 95% confidence intervals
rt1 = sign(Aconfint(1,1)) + sign(Aconfint(2,1));    % do CI95 cross zero?
rt2 = sign(Aconfint(1,2)) + sign(Aconfint(2,2));    % do CI95 cross zero?
m = 0;Acoeffs = coeffvalues(fitA);                  % save model parameters
% Remove outliers if necessary
% Repeat until CI95 do not cross zero and model parameters are positive, or
% until the maximum number of iterations is achieved.
while rt1 == 0 || rt2 == 0 || sign(Acoeffs(1)) == -1 || sign(Acoeffs(2)) == -1
    if m >= maxiter                                  % only loop maxiter times; otherwise, break
        break;
    end
    m = m + 1;                                       % loop index
    fdata = feval(fitA,FreqAll);                     % evaluate the fit
    I = abs(fdata(qmin:qmax) - PSDAll(qmin:qmax)) > 1*std(PSDAll(qmin:qmax));   % find outliers
    outliers = excludedata(FreqAll(qmin:qmax),PSDAll(qmin:qmax),'indices',I);   % remove outliers
    fitopts = fitoptions('Algorithm','Levenberg-Marquardt','Method','NonlinearLeastSquares','TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',4000,'MaxIter',4000,'Robust','Off','Exclude',outliers,'StartPoint',[startptA freqguess]);
    [fitA,gofA] = fit(FreqAll(qmin:qmax),PSDAll(qmin:qmax),fit1,fitopts);
    Aconfint = confint(fitA);                        % check 95% confidence interval
    rt1 = sign(Aconfint(1,1)) + sign(Aconfint(2,1)); % do CI95 cross zero?
    rt2 = sign(Aconfint(1,2)) + sign(Aconfint(2,2)); % do CI95 cross zero?
end
iterA = m;                                           % save number of iterations
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Fit B, from J. Howard           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rt1=0;rt2=0;clear outliers; warning off;                % initialize
fitopts = fitoptions('Algorithm','Levenberg-Marquardt','Method','NonlinearLeastSquares','TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',4000,'MaxIter',4000,'Robust','Off','StartPoint',[startptA freqguess]);
[fitB,gofB] = fit(FreqAll(qmin:qmax),PSDAll(qmin:qmax),fit2,fitopts);
Aconfint = confint(fitB);                               % check 95% confidence intervals
rt1 = sign(Aconfint(1,1)) + sign(Aconfint(2,1));        % do CI95 cross zero?
rt2 = sign(Aconfint(1,2)) + sign(Aconfint(2,2));        % do CI95 cross zero?
m=0;Bcoeffs = coeffvalues(fitB);                        % save model parameters
% Remove outliers if necessary
% Repeat until CI95 do not cross zero and model parameters are positive, or
% until the maximum number of iterations is achieved.
while rt1 == 0 || rt2 == 0 || sign(Bcoeffs(1)) == -1 || sign(Bcoeffs(2)) == -1
    if m >= maxiter  % Only loop maxiter times; otherwise, break
        break;
    end
    m = m + 1;                                          % loop index
    fdata = feval(fitB,FreqAll);                        % evaluate the fit
    I = abs(fdata(qmin:qmax) - PSDAll(qmin:qmax)) > 1*std(PSDAll(qmin:qmax));   % find outliers
    outliers = excludedata(FreqAll(qmin:qmax),PSDAll(qmin:qmax),'indices',I);   % remove outliers
    fitopts = fitoptions('Algorithm','Levenberg-Marquardt','Method','NonlinearLeastSquares','TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',4000,'MaxIter',4000,'Robust','Off','Exclude',outliers,'StartPoint',[startptA freqguess]);
    [fitB,gofB] = fit(FreqAll(qmin:qmax),PSDAll(qmin:qmax),fit2,'Algorithm','Levenberg-Marquardt','Exclude',outliers);
    Aconfint = confint(fitB);                           % check 95% confidence intervals
    rt1 = sign(Aconfint(1,1)) + sign(Aconfint(2,1));    % do CI95 cross zero?
    rt2 = sign(Aconfint(1,2)) + sign(Aconfint(2,2));    % do CI95 cross zero?
end
iterB = m;                                              % save number of iterations
Acoeffs = coeffvalues(fitA);                            % save model parameters for fit A
Bcoeffs = coeffvalues(fitB);                            % save model parameters for fit B
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Fiber parameters              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambdasfA = Acoeffs(1)*2*4.114;         % pN新/nm
ksfA = 2*pi*lambdasfA*Acoeffs(2);       % pN/nm
lambdasfA = lambdasfA*1e6;              % nN新/m
ksfA = ksfA*1e3;                        % 煮/m
    
lambdasfB = 4.114/(2*pi^2*Bcoeffs(1));    % pN新/nm    
ksfB = 2*pi*lambdasfB*Bcoeffs(2);       % pN/nm
lambdasfB = lambdasfB*1e6;              % nN新/m
ksfB = ksfB*1e3;                        % 煮/m

kf = [ksfA ksfB];                       % stiffness
lambdasf = [lambdasfA lambdasfB];       % drag coefficients
R2 = [gofA.rsquare gofB.rsquare];       % R-squared values
CI = [confint(fitA) confint(fitB)];     % 95% confidence intervals
lambdaxx = lambdasf.* 0.94;             % drag owing to motion at tip
lambdadx = lambdasf.*0.57;              % drag owing to motion at base
ksf = kf.*0.97;                         % fiber stiffness
iterations = [iterA iterB];             % number of iterations for each fit
tel = toc;                              % time elapsed
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              findnearest()               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [r,c,V] = findnearest(srchvalue,srcharray,bias)
if nargin<2
    error('Need two inputs: Search value and search array')
elseif nargin<3
    bias = 0;
end
srcharray = srcharray-srchvalue;
if bias == -1  
    srcharray(srcharray>0) =inf;        
elseif bias == 1   
    srcharray(srcharray<0) =inf;        
end

if nargout==1 | nargout==0    
    if all(isinf(srcharray(:)))
        r = [];
    else
        r = find(abs(srcharray)==min(abs(srcharray(:))));
    end         
elseif nargout>1
    if all(isinf(srcharray(:)))
        r = [];c=[];
    else
        [r,c] = find(abs(srcharray)==min(abs(srcharray(:))));
    end    
    if nargout==3
        V = srcharray(r,c)+srchvalue;
    end
end
end
