function importVLCdataclass(datapath,MSVLC)
% This script imports data, without spacefiles, LabVIEW's
% variableloadclamp.vi to MATLAB. You will need to input the full directory
% name containing all files with your data, including the logfile.
%
% importVLCdata(datapath)
%
% datapath : string of the path containing your data
% MSVLC : Mech Stim or VLC? (1=VLC, 2=Mech Stim)
%
% Example:
% importVLCdata('/Users/username/Data/2014-09-30.01/Ear 1/Cell 2/')
%
% In some cases, the logfile will have incorrect indices. This becomes
% obvious for cases where "Fs" (the sampling rate) and "pulse" are not
% their expected values. If this occurs, modify the script and run it
% again.
%
% jsalvi@rockefeller.edu
%

VLCdata = DataImport;

% Find all appropriate files
if MSVLC==1
    file = dir(sprintf('%s%s',datapath,'*VLC*.txt'));   % find all data files
    logfile = dir(sprintf('%s%s',datapath,'*.log'));    % find the logfile
    a = length(file);       % number of sessions in the directory
elseif MSVLC==2
    file = dir(sprintf('%s%s',datapath,'*MS*.txt'));   % find all data files
    logfile = dir(sprintf('%s%s',datapath,'*.log'));    % find the logfile
    a = length(file);       % number of sessions in the directory
end
%{
% Find all appropriate files
file = dir(sprintf('%s%s',datapath,'*VLC*.txt'));   % find all data files
logfile = dir(sprintf('%s%s',datapath,'*.csv'));    % find the logfile
a = length(file);       % number of sessions in the directory
%}
% Import logdata
logdata = importdata(sprintf('%s%s',datapath,logfile.name));  % logdata.data contains log data of interest
if isstruct(logdata) == 0
%    logdata.data = logdata;
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


% Note that some logfiles will have formatting issues and the indices
% listed above may be off by �10 in certain cases. If this is true, open
% logdata.data and modify the script accordingly.



%Number of traces
ntrace=logdata.data(isnan(logdata.data(:,8))==0,8);
numavg=logdata.data(isnan(logdata.data(:,8))==0,10);
% Find all raw files
for j = 1:a
    rawfiles(j) = isempty(findstr(file(j).name,'raw'));
end
nonraw = find(rawfiles==1);raw = find(rawfiles==0);
ntraceraw = ntrace(ntrace~=0);
numavgraw = numavg(numavg~=0);

% Import the data, some of this may be redundant
for i = 1:a
    data2{i}=importdata(sprintf('%s%s',datapath,file(i).name));   % initial import
end

for i =1:a
    data0{i}=data2{i}.data;     % extract data from its structure format (not necessary, but easier)
end
clear data2 data 
%}
% Import non-raw time traces
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

% Extract pulses from full time traces
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
% Import non-raw time traces
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


% Extract pulses from full time traces

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

% Time vector
dt = 1/Fs;
sizeX = size(Xd);
for j = 1:sizeX(1)
    for k = 1:sizeX(2)
        tvec{j,k} = 0:dt:length(Xd{j,k})*dt-dt;
    end
end
clear i j

VLCdata = VLCdata.initialization(datapath,Fs,pre,pulse,cyclesdur,cyclesdurlength,post,ntrace,numavg,raw,nonraw);
VLCdata = VLCdata.logfileimport(logdata);
VLCdata = VLCdata.timeseries(tvec,Xd,Xo,Fe,Xd_pulse,Xo_pulse,Fe_pulse,Xd_split,Xo_split,kv);
VLCdata.savedata(VLCdata.datapath);

end
