function importVLCdata(datapath)
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
comments = logdata.textdata{2}; % import comments
Fs = logdata.data(1,12);       % scan rate (Hz), CHECK THIS!
pre = logdata.data(1,22)*1e-3*Fs;   % delay before a stimulus, CHECK THIS!
pulse = logdata.data(1,23)*1e-3*Fs; % length of stimulus, CHECK THIS!
% Note that some logfiles will have formatting issues and the indices
% listed above may be off by �10 in certain cases. If this is true, open
% logdata.data and modify the script accordingly.


% Import the data, some of this may be redundant
for i = 1:a
data=importdata(sprintf('%s%s',datapath,file(i).name));   % initial import
 
if i<a
    data2.data(:,:,i)=data.data;    % assign with proper indexing
else
    if length(data.data(:,1,1))<length(data2.data(:,1,1))   % search for length discrepancies and initialize with zeros
        data.data(length(data2.data(:,1,1)),:)=0;
        for j = (length(data.data(1,:))+1):length(data2.data(1,:,(i-1)))
           data.data(:,j)=zeros(1,length(data2.data(:,:)));
        end
        data2.data(:,:,i)=data.data(1:length(data2.data(:,1,1)),:);
    else
        for j = (length(data.data(1,:))+1):length(data2.data(1,:,(i-1)))
           data.data(:,j)=zeros(1,length(data2.data(:,:)));
        end
        data2.data(:,:,i)=data.data(1:length(data2.data(:,1,1)),:);
    end
end
 
end

time = data2.data(:,1);    % time vector

for i =1:a
    data0(:,:,i)=data2.data(:,:,i);     % extract data from its structure format (not necessary, but easier)
end
clear data2 data 

% Import time traces
for j = 1:a
    for i = 1:(logdata.data(1,1))
        Xd(:,i,j) = data0(:,(1+i),j);   % photodiode
        Xo(:,i,j) = data0(:,(1+logdata.data(1,1)+i),j);  % stimulus piezo
    end
end

% Extract pulses from full time traces
for j = 1:a
    for i = 1:(logdata.data(1,1))
        Xd_pulse(:,i,j) = Xd(1+pre:pre+pulse,i,j);  % photodiode, stimulation only
        Xo_pulse(:,i,j) = Xo(1+pre:pre+pulse,i,j);  % stimulus piezo, stimulation only
    end
end

% Save extracted data
disp('Saving...');
savefile = sprintf('%s%s',datapath,'Extracted Data.mat');
save(savefile);
disp('Finished.');

end
