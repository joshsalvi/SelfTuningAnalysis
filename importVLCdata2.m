function importVLCdata2(datapath)
% This script imports data, without spacefiles, LabVIEW's
% variableloadclamp.vi to MATLAB. You will need to input the full directory
% name containing all files with your data, including the logfile.
%
% importVLCdata(datapath)
%
% datapath : string of the path containing your data
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
file = dir(sprintf('%s%s',datapath,'*VLC*.txt'));   % find all data files
logfile = dir(sprintf('%s%s',datapath,'*.log'));    % find the logfile
a = length(file);       % number of sessions in the directory

% Import logdata
logdata = importdata(sprintf('%s%s',datapath,logfile.name));  % logdata.data contains log data of interest
comments = logdata.textdata(isnan(logdata.data(:,8))==0,3); % import comments
Fs = logdata.data(1,12);       % scan rate (Hz), CHECK THIS!
for j = 1:a
pre(j) = logdata.data(j,22)*1e-3*Fs;   % delay before a stimulus, CHECK THIS!
pulse(j) = logdata.data(j,23)*1e-3*Fs; % length of stimulus, CHECK THIS!
end
% Note that some logfiles will have formatting issues and the indices
% listed above may be off by ±10 in certain cases. If this is true, open
% logdata.data and modify the script accordingly.

%Number of traces
ntrace=logdata.data(isnan(logdata.data(:,8))==0,8);
% Find all raw files
for j = 1:a
    rawfiles(j) = isempty(findstr(file(j).name,'raw'));
end
nonraw = find(rawfiles==1);raw = find(rawfiles==0);
for j = 1:length(raw)
    ntraceraw(j) = ntrace(nonraw==raw(j)+1);
end

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
    k=k+1;
end
% Extract pulses from full time traces
k=1;

for j = nonraw
    for i = 1:ntrace(k)
        Xd_pulse{i,j} = Xd{i,j}(1+pre(j):pre(j)+pulse(j));  % photodiode, stimulation only
        Xo_pulse{i,j} = Xo{i,j}(1+pre(j):pre(j)+pulse(j));  % stimulus piezo, stimulation only
    end
    k=k+1;
end

% Import non-raw time traces
k=1;
for j = raw
    for i = 1:ntraceraw(k)
        Xd{i,j} = data0{j}(:,(i));   % photodiode
        Xo{i,j} = data0{j}(:,(ntraceraw(k)+i));  % stimulus piezo
    end
    k=k+1;
end
% Extract pulses from full time traces
k=1;
for j = raw
    for i = 1:ntraceraw(k)
        Xd_pulse{i,j} = Xd{i,j}(1+pre:pre+pulse);  % photodiode, stimulation only
        Xo_pulse{i,j} = Xo{i,j}(1+pre:pre+pulse);  % stimulus piezo, stimulation only
    end
    k=k+1;
end


clear i j 

% Save extracted data
disp('Saving...');
savefile = sprintf('%s%s',datapath,'Extracted Data.mat');
save(savefile);
disp('Finished.');

end
