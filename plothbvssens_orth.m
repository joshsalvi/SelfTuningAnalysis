function plothbvssens_orth(datapath)
% This function calculates the sensitivity and vector strength across
% various sets of control parameters.
%
% plothbvssens_orth(datapath)
%
% Joshua D. Salvi
% jsalvi@rockefeller.edu

clc;
setfiguredefaults
close all

% Search for "Extracted Data.mat"
filename = dir(sprintf('%s%s',datapath,'Extracted Data.mat'));

% Load the data file
% If data file not present, import it (with user prompt)
if isempty(filename) == 0
	disp('MAT file found! Loading...');
    load([datapath filename.name]);
    ramp = logdata.data(:,41)*Fs;
    disp('Done.');
else
    disp('No MAT file found!');
    importdatayn = input('Import (1=yes)?  ');
    if importdatayn == 1
        importVLCdata3(datapath,1);
        filename = dir(sprintf('%s%s',datapath,'Extracted Data.mat'));
        load([datapath filename.name]);
        ramp = logdata.data(:,41)*Fs;
        disp('Done.');
    else
        disp('Exiting...');
        return;
    end
end

analyzeyn=1;
close all;

% Select indices for analysis
ksf = input('Fiber stiffness (µN/m): ');
disp(['Indices: ' num2str(nonraw)]);
disp('Select a range of indices for analysis...');
ind = input('Indices (e.g. 1:5 or [1 2 4 5]: ');
if isempty(intersect(ind,nonraw)) == 0
    disp(['Selected indices ' num2str(ind) '.']);
else
    disp('Did not selected appropriate index.');
    return;
end

% Ask user to input stimulus data
% Assume linear spacing
numfreq = input('Number of frequencies: ');
minfreq = input('Minimum frequency: ');
maxfreq = input('Maximum frequency: ');
freqrange = linspace(minfreq,maxfreq,numfreq);

% Time vector
clear dt tvec
dt = 1/Fs;
tvec = 0:dt:length(Xd_pulse{1,ind(1)})*dt-dt;

% Plot the traces
setfiguredefaults(length(ind));    
numpoints = size(Xd_pulse,1);
numtraces = length(ind);

XsegL = floor(length(tvec));

% Split the time series into individual waveforms
% NOTE: This assumes that there is no delay 
disp('Splitting waveforms.');
clear Xd_split Xo_split
stimlength = length(Xd_pulse{1,ind(1)}) / numfreq;
for j = 1:numpoints
    for k = 1:numtraces
        stimlength = length(Xd_pulse{j,ind(k)}) / numfreq;
        for l = 1:numfreq
            warning off
            Xd_split{j,k,l} = Xd_pulse{j,ind(k)}(1+(l-1)*stimlength:l*stimlength);
            Xo_split{j,k,l} = Xo_pulse{j,ind(k)}(1+(l-1)*stimlength:l*stimlength);
            Xd_split{j,k,l} = Xd_split{j,k,l} - mean(Xd_split{j,k,l});
            Xo_split{j,k,l} = Xo_split{j,k,l} - mean(Xo_split{j,k,l});
        end
    end
end


% FFT Preliminaries
NFFT = (2^1)*2^nextpow2(numel(tvec));   % CHOOSE number of FFT points (2^n)*(n_p)
nw = 1;     % one window

welchwin = round(XsegL);
NPSD = floor(NFFT/nw);
noverlap = 0;       % zero overlap
winfunc = rectwin(welchwin);    % window/taper
f = Fs/2*linspace(0,1,NFFT/2+1);
% Test function
freq = 1;
Xsine = sin(2*pi*freq.*tvec);
Xsinefft = fft(Xsine,NFFT)./XsegL; Xsinefft = Xsinefft(1:NFFT/2+1);
winpeaknorm = sqrt(max(abs(Xsinefft)).*(2.*Fs.*XsegL.*(sum(abs(winfunc).^2)./XsegL)))./XsegL;

% Calculate FFT for Xo and Xd
disp(['Calculating FFT with ' num2str(NFFT) ' points.']);
m = 0;
disp([num2str(round(m/(numpoints*numtraces))) '% complete']);
for j = 1:numpoints
    for k = 1:numtraces
        for l = 1:numfreq
            Xd_fft{j,k,l} = fft(Xd_split{j,k,l},NFFT)./stimlength; Xd_fft{j,k,l} = Xd_fft{j,k,l}(1:NFFT/2+1);
            Xo_fft{j,k,l} = fft(Xo_split{j,k,l},NFFT)./stimlength; Xo_fft{j,k,l} = Xo_fft{j,k,l}(1:NFFT/2+1);
            freq0 = findnearest(f,freqrange(l)); freq0 = freq0(1);
            %Xdpk(j,k,l) = Xd_fft{j,k,l}(abs(Xd_fft{j,k,l}(freq0-5:freq0+5)) == max(abs(Xd_fft{j,k,l}(freq0-5:freq0+5))));
            %Xopk(j,k,l) = Xo_fft{j,k,l}(abs(Xo_fft{j,k,l}(freq0-5:freq0+5)) == max(abs(Xo_fft{j,k,l}(freq0-5:freq0+5))));
            Xopk(j,k,l) = Xo_fft{j,k,l}(freq0);
            Xdpk(j,k,l) = Xd_fft{j,k,l}(freq0);
            %freqstim_check(j,k,l) = f(abs(Xd_fft{j,k,l}(freq0-5:freq0+5)) == max(abs(Xd_fft{j,k,l}(freq0-5:freq0+5))));
        end
        m = m+1;
        disp([num2str(round(m/(numpoints*numtraces)*100)) '% complete']);
    end
end
while analyzeyn == 1
    clear hrt
    close all
% Plot the individual results
disp('Plotting individual results.');
setfiguredefaults(numpoints);
hrt(1) = figure(1);
for k = 1:numtraces
    subplot(1,numtraces,k)
    for j = 1:numpoints
        plot(freqrange,abs(sqrt(squeeze(Xdpk(j,k,:)))));
        hold all;
        title(['Trial ' num2str(k)]);
        xlabel('Frequency (Hz)');ylabel('X_D (nm)');
        set(1,'WindowStyle','docked')
    end
end

plotfftyn = input('Plot spectra (1=yes)?  ');
if plotfftyn == 1
for m = 1:numfreq
hrt(m+5) = figure(m+5);
f_ind = freqrange(m);
for k = 1:numtraces
    subplot(1,numtraces,k)
    for j = 1:numpoints
        plot(f,abs(Xd_fft{j,k,m}));
        hold all;
        title(['freq ' num2str(f_ind)]);
        xlabel('Frequency (Hz)');ylabel('X_f_f_t');
        set(m+5,'WindowStyle','docked')
    end
end
end
end

analyzedind = input('Indices to analyze (e.g. 1:6, [2 3 5]):  ');

for j = 1:numpoints
    for l = 1:numfreq
        Xdpkmean(j,l) = abs(mean(Xdpk(j,analyzedind,l))) * 1e-9;
        Fepkmean(j,l) = abs(mean(Xopk(j,analyzedind,l))) * ksf * 1e-15;
        sensmean(j,l) = Xdpkmean(j,l) / Fepkmean(j,l);
        Xdpksem(j,l) = abs(std(Xdpk(j,analyzedind,l)))*1e-9/sqrt(length(analyzedind));
        Fepksem(j,l) = abs(std(Xopk(j,analyzedind,l)))*ksf*1e-15/sqrt(length(analyzedind));
        senssem(j,l) = sensmean(j,l) * sqrt( (Xdpksem(j,l)/Xdpkmean(j,l))^2 + (Fepksem(j,l)/Fepkmean(j,l))^2 );
    end
end

hrt(3) = figure(3);
for j = 1:numpoints
    errorbar(freqrange,Xdpkmean(j,:),Xdpksem(j,:)); hold all;
end
xlabel('Frequency (Hz)'); ylabel('Response amplitude (m)');
set(3,'WindowStyle','docked')

hrt(2) = figure(2);
for j = 1:numpoints
    errorbar(freqrange,sensmean(j,:),senssem(j,:)); hold all;
end
xlabel('Frequency (Hz)'); ylabel('Sensitivity (m/N)');
set(2,'WindowStyle','docked')



saveyn = input('Save the figures? (1=yes): ');
hrt(hrt==0)=[];
if saveyn==1
    disp('Saving...');
    savefig(hrt,[datapath 'Sens-ind' num2str(analyzedind) '-figures.fig']);
    disp('Finished.');
else
    disp('Not saved.');
end
saveyn = input('Save the data? (1=yes): ');
if saveyn==1
    disp('Saving...');
    save([datapath 'Sens-ind' num2str(analyzedind) '-data.mat'],'datapath','analyzedind','minfreq','maxfreq','numfreq','freqrange','Xdpk','Xopk','Xdpkmean','Fepkmean','Xdpksem','Fepksem','sensmean','senssem','Xd_fft','Xo_fft','f','Xd_split','Xo_split');
    disp('Finished.');
else
    disp('Not saved.');
end

analyzeyn = input('Analyze again? (1=yes):  ');
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%  setfiguredefaults()  %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function setfiguredefaults(N)

set(0,'defaultFigureColormap',jet)
if exist('N')==1
    set(0,'defaultAxesColorOrder',ametrine(N))
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  importVLCdata3()  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function importVLCdata3(datapath,MSVLC)

if MSVLC==1
    file = dir(sprintf('%s%s',datapath,'*VLC*.txt'));   % find all data files
    logfile = dir(sprintf('%s%s',datapath,'*.log'));    % find the logfile
    a = length(file);       % number of sessions in the directory
elseif MSVLC==2
    file = dir(sprintf('%s%s',datapath,'*MS*.txt'));   % find all data files
    logfile = dir(sprintf('%s%s',datapath,'*.log'));    % find the logfile
    a = length(file);       % number of sessions in the directory
end
logdata = importdata(sprintf('%s%s',datapath,logfile.name));  % logdata.data contains log data of interest
if isstruct(logdata) == 0
end
if isempty('logdata.textdata(isnan(logdata.data(:,8))==0,3)')==0
    comments = logdata.textdata(isnan(logdata.data(:,8))==0,3); % import comments
end
Fs = logdata.data(1,12);       % scan rate (Hz), CHECK THIS!
pre = logdata.data(:,22)*1e-3*Fs;   % delay before a stimulus, CHECK THIS!
pulse = logdata.data(:,23)*1e-3*Fs; % length of stimulus, CHECK THIS!
cyclesdur = logdata.data(:,24)*1e-3*Fs; % length of stimulus, CHECK THIS!
pulselength = pre + pulse;
cyclesdurlength = pre + cyclesdur;
post = pulselength - pre - pulse;
ntrace=logdata.data(isnan(logdata.data(:,8))==0,8);
numavg=logdata.data(isnan(logdata.data(:,8))==0,10);
for j = 1:a
    rawfiles(j) = isempty(findstr(file(j).name,'raw'));
end
nonraw = find(rawfiles==1);raw = find(rawfiles==0);
ntraceraw = ntrace(ntrace~=0);
numavgraw = numavg(numavg~=0);
for i = 1:a
    data2{i}=importdata(sprintf('%s%s',datapath,file(i).name));   % initial import
end

for i =1:a
    data0{i}=data2{i}.data;     % extract data from its structure format (not necessary, but easier)
end
clear data2 data 
k=1;
for j = nonraw
    for i = 1:ntrace(k)
        Xd{i,j} = data0{j}(:,(1+i));   % photodiode
        Xo{i,j} = data0{j}(:,(1+ntrace(k)+i));  % stimulus piezo
       
    end
    if size(data0{j},2)>1+2*ntrace(k)
        sizen = size(data0{j},2);
        for i = 1:size(data0{j},2)-(1+2*ntrace(k))
            Fe{i,j} = data0{j}(:,(1+2*ntrace(k)+i));
        end
    end
    k=k+1;
end
k=1;
for j = nonraw
    for i = 1:ntrace(k)
        Xd_pulse{i,j} = Xd{i,j};  % photodiode, stimulation only
        Xo_pulse{i,j} = Xo{i,j};  % stimulus piezo, stimulation only
        if size(data0{j},2)>1+2*ntrace(k)
            Fe_pulse{i,j} = Fe{i,j};
        end
    end
    k=k+1;
end
k=1;
for j = raw
    for i = 1:ntraceraw(k)
        Xd{i,j} = data0{j}(:,1);   % photodiode
        Xo{i,j} = data0{j}(:,2);  % stimulus piezo
    end
    if size(data0{j},2)>1+2*ntrace(k)
        sizen = size(data0{j},2);
        for i = 1:size(data0{1},2)-(1+2*ntrace(k))
            Fe{i,j} = data0{j}(:,(1+2*ntrace(k)+i));
        end
    end
    k=k+1;
end
if exist('Fe')
if size(Fe,1) > size(Xd,1)
    for j = 1:size(Xd,1)
        for k = nonraw
            kv{j,k} = Fe{j,k};
            Fe2{j,k} = Fe{j+size(Xd,1),k};
        end
    end
    Fe = Fe2; clear Fe2
end
end
mm=1;
for j = raw
for k=1:numavgraw(mm)
    lengthxdraw(mm) = length(Xd{1,j})/numavgraw(mm);
    Xd_split{k,j} = Xd{1,j}(1+(k-1)*lengthxdraw(mm):k*lengthxdraw(mm));
    Xo_split{k,j} = Xo{1,j}(1+(k-1)*lengthxdraw(mm):k*lengthxdraw(mm));
end
mm=mm+1;
end
mm=1;
for j = raw
    for k = 1:numavgraw(mm)
        for i = 1:ntraceraw(mm)
            trange=(1+(i-1)*(pre(mm)+pulse(mm)+post(mm))):(i*(pre(mm)+pulse(mm)+post(mm)));
            Xd_pulse{k,j}(:,i) = Xd_split{k,j}(trange);  % photodiode, stimulation only
            Xo_pulse{k,j}(:,i) = Xo_split{k,j}(trange);  % stimulus piezo, stimulation only    
        end
    end
    mm=mm+1;
end
dt = 1/Fs;
sizeX = size(Xd);
for j = 1:sizeX(1)
    for k = 1:sizeX(2)
        tvec{j,k} = 0:dt:length(Xd{j,k})*dt-dt;
    end
end
clear i j 
% Save extracted data
disp('Saving...');
savefile = sprintf('%s%s',datapath,'Extracted Data.mat');
save(savefile);
disp('Finished.');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  findnearest()  %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  ametrine()  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cmap=ametrine(n,varargin)

p=inputParser;
p.addParamValue('gamma',1.8, @(x)x>0);
p.addParamValue('minColor','none');
p.addParamValue('maxColor','none');
p.addParamValue('invert',0, @(x)x==0 || x==1);
if nargin==1
    p.addRequired('n', @(x)x>0 && mod(x,1)==0);
    p.parse(n);
elseif nargin>1
    p.addRequired('n', @(x)x>0 && mod(x,1)==0);
    p.parse(n, varargin{:});
else
    p.addParamValue('n',256, @(x)x>0 && mod(x,1)==0);
    p.parse();
end
config = p.Results;
n=config.n;
cP(:,1) = [30  60  150]./255; k(1)=1;  %cyan at index 1
cP(:,2) = [180 90  155]./255; k(3)=17; %purple at index 17
cP(:,3) = [230 85  65 ]./255; k(4)=32; %redish at index 32
cP(:,4) = [220 220 0  ]./255; k(5)=64; %yellow at index 64
for i=1:3
    f{i}   = linspace(0,1,(k(i+1)-k(i)+1))';  % linear space between these controlpoints
    ind{i} = linspace(k(i),k(i+1),(k(i+1)-k(i)+1))';
end
cmap = interp1((1:4),cP',linspace(1,4,64)); % for non-iso points, a normal interpolation gives better results
cmap = abs(interp1(linspace(0,1,size(cmap,1)),cmap,linspace(0,1,n)));
if config.invert
    cmap = flipud(cmap);
end
if ischar(config.minColor)
    if ~strcmp(config.minColor,'none')
        switch config.minColor
            case 'white'
                cmap(1,:) = [1 1 1];
            case 'black'
                cmap(1,:) = [0 0 0];
            case 'lightgray'
                cmap(1,:) = [0.8 0.8 0.8];
            case 'darkgray'
                cmap(1,:) = [0.2 0.2 0.2];
        end
    end
else
    cmap(1,:) = config.minColor;
end
if ischar(config.maxColor)
    if ~strcmp(config.maxColor,'none')
        switch config.maxColor
            case 'white'
                cmap(end,:) = [1 1 1];
            case 'black'
                cmap(end,:) = [0 0 0];
            case 'lightgray'
                cmap(end,:) = [0.8 0.8 0.8];
            case 'darkgray'
                cmap(end,:) = [0.2 0.2 0.2];
        end
    end
else
    cmap(end,:) = config.maxColor;
end
end
