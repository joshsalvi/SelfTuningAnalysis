function importVLCdata4(datapath,MSVLC)
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

% Find all appropriate files
file = dir(sprintf('%s%s',datapath,'*VLC*.txt'));   % find all data files
logfile = dir(sprintf('%s%s',datapath,'*.log'));    % find the logfile
a = length(file);       % number of sessions in the directory

% Import logdata
logdata = importdata(sprintf('%s%s',datapath,logfile.name));  % logdata.data contains log data of interest
%comments = logdata.textdata(isnan(logdata.data(:,8))==0,3); % import comments
Fs = logdata.data(1,12);       % scan rate (Hz), CHECK THIS!
pre = logdata.data(1,22)*1e-3*Fs;   % delay before a stimulus, CHECK THIS!
pulse = logdata.data(1,23)*1e-3*Fs; % length of stimulus, CHECK THIS!
% Note that some logfiles will have formatting issues and the indices
% listed above may be off by �10 in certain cases. If this is true, open
% logdata.data and modify the script accordingly.

%Number of traces
ntrace=logdata.data(isnan(logdata.data(:,8))==0,9);
% Find all raw files
for j = 1:a
    rawfiles(j) = isempty(findstr(file(j).name,'raw'));
end
nonraw = find(rawfiles==1);raw = find(rawfiles==0);

j=1;
    ntraceraw = ntrace;

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
    for i = 1:ntrace
        Xd{i,j} = data0{j}(:,(1+i));   % photodiode
        Xo{i,j} = data0{j}(:,(1+ntrace+i));  % stimulus piezo
    end
    k=k+1;
end
% Extract pulses from full time traces
k=1;
for j = nonraw
    for i = 1:ntrace
        Xd_pulse{i,j} = Xd{i,j}(1+pre:pre+pulse);  % photodiode, stimulation only
        Xo_pulse{i,j} = Xo{i,j}(1+pre:pre+pulse);  % stimulus piezo, stimulation only
    end
    k=k+1;
end
% Import raw time traces
k=1;
for j = raw
    for i = 1:size(data0{j},2)
        Xd{i,j} = data0{j}(:,(i));   % photodiode
        Xo{i,j} = data0{j}(:,(i));  % stimulus piezo
    end
    k=k+1;
end
% Extract pulses from full time traces

k=1;
for j = raw
    for i = 1:size(Xd,1)
        Xd_pulse{i,j} = Xd{i,j}(1+pre:pre+pulse);  % photodiode, stimulation only
        Xo_pulse{i,j} = Xo{i,j}(1+pre:pre+pulse);  % stimulus piezo, stimulation only
    end
    k=k+1;
end
%}

clear i j 

% Save extracted data
disp('Saving...');
savefile = sprintf('%s%s',datapath,'Extracted Data.mat');
save(savefile);
disp('Finished.');

end
