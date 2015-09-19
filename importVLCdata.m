function importVLCdata(datapath,MSVLC)
% This script imports data, without spacefiles, LabVIEW's
% variableloadclamp.vi to MATLAB. You will need to input the full directory
% name containing all files with your data, including the logfile.
%
% importVLCdata(datapath,MSVLC)
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

% Import logdata
logdata = importdata(sprintf('%s%s',datapath,logfile.name));  % logdata.data contains log data of interest
%comments = logdata.textdata{2}; % import comments
Fs = logdata.data(1,12);       % scan rate (Hz), CHECK THIS!
pre = logdata.data(1,22)*1e-3*Fs;   % delay before a stimulus, CHECK THIS!
pulse = logdata.data(1,23)*1e-3*Fs; % length of stimulus, CHECK THIS!
% Note that some logfiles will have formatting issues and the indices
% listed above may be off by ±10 in certain cases. If this is true, open
% logdata.data and modify the script accordingly.


% Import the data, some of this may be redundant
for i = 1:a
data=importdata(sprintf('%s%s',datapath,file(i).name));   % initial import
data2.data{i} = data.data;
%{
if a==1
    data2.data{i}=data.data; 
elseif i<a
    data2.data{i}=data.data;    % assign with proper indexing
else
    if length(data.data(:,1,1))<length(data2.data{1}(:,1))   % search for length discrepancies and initialize with zeros
        data.data(length(data2.data{i}(:,1)),:)=0;
        for j = (length(data.data(1,:))+1):length(data2.data(1,:,(i-1)))
           data.data(:,j)=zeros(1,length(data2.data(:,:)));
        end
        data2.data{i}=data.data(1:length(data2.data{i}(:,1)),:);
    else
        for j = (length(data.data(1,:))+1):length(data2.data{i-1}(1,:))
           data.data(:,j)=zeros(1,length(data2.data{1}(:,:)));
        end
        data2.data{i}=data.data(1:length(data2.data{1}(:,1)),:);
    end
end
%}
time{i}=data2.data{i}(:,1); 
end

clear data 

% Import time traces
if MSVLC==1
    for j = 1:a
    for i = 1:(logdata.data(1,8))
        Xd{i,j} = data2.data{j}(:,(1+i));   % photodiode
        Xo{i,j} = data2.data{j}(:,(1+logdata.data(1,8)+i));  % stimulus piezo
    end
    end
elseif MSVLC==2
    for j = 1:a
    for i = 1:(logdata.data(1,8))
        Xo{i,j} = data2.data{j}(:,(1+i));   % photodiode
        Xd{i,j} = data2.data{j}(:,(1+logdata.data(1,8)+i));  % stimulus piezo
    end
    end
end


% Extract pulses from full time traces
for j = 1:a
    for i = 1:(logdata.data(1,8))
        Xd_pulse{i,j} = Xd{i,j}(1+pre:pre+pulse);  % photodiode, stimulation only
        Xo_pulse{i,j} = Xo{i,j}(1+pre:pre+pulse);  % stimulus piezo, stimulation only
    end
end
%}

% Save extracted data
disp('Saving...');
savefile = sprintf('%s%s',datapath,'Extracted Data.mat');
save(savefile);
disp('Finished.');

end
