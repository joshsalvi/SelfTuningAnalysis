function importPCdata(datapath)
% Imports data from probe characterization
%
% importPCdata(datapath)
%
%
% Joshua D. Salvi
% jsalvi@rockefeller.edu
%
%

file = dir(sprintf('%s%s',datapath,'*PC*.txt'));   % find all data files
logfile = dir(sprintf('%s%s',datapath,'*.log'));    % find the logfile
a = length(file);       % number of sessions in the directory

% Import the data, some of this may be redundant
for i = 1:a
    data2{i}=importdata(sprintf('%s%s',datapath,file(i).name));   % initial import
end
for i =1:a
    data0{i}=data2{i}.data;     % extract data from its structure format (not necessary, but easier)
end

% Import power spectra
m=1;
for j = 1:2:a
    FreqAll{m} = data0{j}(:,1);
    PSDAll{m} = data0{j}(:,2);
    FreqNoO{m} = data0{j}(:,3);
    PSDNoO{m} = data0{j}(:,4);
    FreqFit{m} = data0{j}(:,5);
    PSDFit{m} = data0{j}(:,6);
    m = m+1;
end

% Import time-series data
m = 1;
for j = 2:2:a
    Timems{m} = data0{j}(:,1);
    RawPD{m} = data0{j}(:,2);
    SmoothedPD{m} = data0{j}(:,3);
    m = m+1;
end
clear data0 data2

% Save extracted data
disp('Saving...');
savefile = sprintf('%s%s',datapath,'Extracted Data.mat');
save(savefile);
disp('Finished.');

end
