function importVLCdata5(varargin)
% This script imports data, without spacefiles, LabVIEW's
% variableloadclamp.vi to MATLAB. You will need to input the full directory
% name containing all files with your data, including the logfile.
%
% importVLCdata(datapath,varsincluded)
%
% -> datapath : string of the path containing your data
% -> varsincluded : list of numbers (1-7) that states which variables were
% recorded. Use the following key:
%       (1):  Xd
%       (2):  Xo
%       (3):  Xc
%       (4):  kv
%       (5):  Fe
%       (6):  mv
%       (7):  gv
%
% EXAMPLES:
% importVLCdata('/Users/username/Data/2014-09-30.01/Ear 1/Cell 2/')
%   -> Datapath is a string. Defaults to all variables with only one arg.
% importVLCdata('/Users/username/Data/2014-09-30.01/Ear 1/Cell 2')
%   -> No '/' at end. OK. Defaults to all variables.
% importVLCdata('/Users/username/Data/2014-09-30.01/Ear 1/Cell 2', 1:7)
%   -> Includes all possible variables
% importVLCdata('/Users/username/Data/2014-09-30.01/Ear 1/Cell 2', [1 2 4])
%   -> Only includes variables Xd, Xo, and kv
%
% In some cases, the logfile will have incorrect indices. This becomes
% obvious for cases where "Fs" (the sampling rate) and "pulse" are not
% their expected values. If this occurs, modify the script and run it
% again.
%
% Coder:    Joshua D. Salvi
% Email:    jsalvi@rockefeller.edu
% Year:     2016
%

datapath = varargin{1};

if nargin == 1
    varsincluded = 1:7;
else
    varsincluded = varargin{2};
end

if datapath(end) ~= '/'
    datapath(end+1) = '/';
end

% Find all appropriate files
file = dir(sprintf('%s%s',datapath,'*VLC*.txt'));   % find all data files
logfile = dir(sprintf('%s%s',datapath,'*.log'));    % find the logfile
a = length(file);       % number of sessions in the directory
logfilekey = '/Users/joshsalvi/Documents/Lab/Lab/VLClogfilekey.xlsx';
if a == 0
    disp('No files found. Aborting.');
    return;
end

% Import logdata
try
    logdata = importdata(sprintf('%s%s',datapath,logfile.name));  % logdata.data contains log data of interest
    logstruc = getVLClog(datapath,logdata,logfilekey);
catch
    disp('No logfile found. Aborting.');
    return;
end

if isstruct(logdata) == 0
   logdata.data = logdata;
end
if isempty('logdata.textdata(isnan(logdata.data(:,8))==0,3)')==0
    try
        comments = logdata.textdata(:,3); % import comments
    catch
        disp('Cannot import comments.')
    end
end

% Label the comments
m=1;
for j = 1:length(comments)
    try
        if isempty(comments{j})
        else
            try
                comments{j} = ['(' num2str(m) '): ' comments{j}];
                m=m+1;
            catch
            end
        end
    catch
    end
end

try
    Fs = logdata.data(1,12);       % scan rate (Hz), CHECK THIS!
    pre = logdata.data(:,22)*1e-3*Fs;   % delay before a stimulus, CHECK THIS!
    pulse = logdata.data(:,23)*1e-3*Fs; % length of stimulus, CHECK THIS!
    cyclesdur = logdata.data(:,24)*1e-3*Fs; % length of stimulus, CHECK THIS!
    pulselength = pre + pulse;
    cyclesdurlength = pre + cyclesdur;
    post = pulselength - pre - pulse;
catch
    disp('Unable to get information from logfile. Aborting.');
    return;
end

% Note that some logfiles will have formatting issues and the indices
% listed above may be off by �10 in certain cases. If this is true, open
% logdata.data and modify the script accordingly.



%Number of traces
try
    ntrace=logdata.data(isnan(logdata.data(:,8))==0,8);
    numavg=logdata.data(isnan(logdata.data(:,8))==0,10);
catch
    disp('Unable to find number of traces in logfile. Aborting.');
    return;
end



% Find all raw files
rawfiles = zeros(1,a);
for j = 1:a
    try
        rawfiles(j) = isempty(findstr(file(j).name,'raw'));
    catch
        disp('Unable to find RAW files.');
    end
end

try
    nonraw = find(rawfiles==1);raw = find(rawfiles==0);
    ntraceraw = ntrace(ntrace~=0);
    numavgraw = numavg(numavg~=0);
catch
    disp('Unable to find raw or non-raw files. Aborting.');
    return;
end

% Import the data, some of this may be redundant
data2 = cell(1,a);
data0 = cell(1,a);
for i = 1:a
    try
        data2{i} = importdata(sprintf('%s%s',datapath,file(i).name));   % initial import
        data0{i} = data2{i}.data;
    catch
    end
end
clear data2 data 

% Import each average of the time traces
k=1;

while length(numavg) < a
    numavg(end+1) = numavg(end);
end

% Pre-allocate variables
Xd = cell(max(numavg),a);
Xo = cell(max(numavg),a);
Xc = cell(max(numavg),a);
kv = cell(max(numavg),a);
mv = cell(max(numavg),a);
gv = cell(max(numavg),a);
Fe = cell(max(numavg),a);

for j = 1:a
    xL = length(data0{j})/numavg(j);   
    for i = 1:numavg(j)
        m=1;
        try
            if isempty(intersect(varsincluded,1)) == 0
                Xd{i,j} = data0{j}(1+(i-1)*xL:i*xL,m);
                m = m + 1;
            end
        catch
        end
        
        try
            if isempty(intersect(varsincluded,2)) == 0
                Xo{i,j} = data0{j}(1+(i-1)*xL:i*xL,m);
                m = m + 1;
            end
        catch
        end
        try
            if isempty(intersect(varsincluded,3)) == 0
                Xc{i,j} = data0{j}(1+(i-1)*xL:i*xL,m);
                m = m + 1;
            end
        catch
        end
        try
            if isempty(intersect(varsincluded,4)) == 0
                kv{i,j} = data0{j}(1+(i-1)*xL:i*xL,m);
                m = m + 1;
            end
        catch
        end
        try
            if isempty(intersect(varsincluded,5)) == 0
                Fe{i,j} = data0{j}(1+(i-1)*xL:i*xL,m);
                m = m + 1;
            end
        catch
        end
        try
            if isempty(intersect(varsincluded,6)) == 0
                mv{i,j} = data0{j}(1+(i-1)*xL:i*xL,m);
                m = m + 1;
            end
        catch
        end
        try
            if isempty(intersect(varsincluded,7)) == 0
                gv{i,j} = data0{j}(1+(i-1)*xL:i*xL,m);
                m = m + 1;
            end
        catch
        end
        
        k=k+1;
    end
end

% Extract raw pulses from full time traces
try    
    while length(ntraceraw) < length(raw)
        ntraceraw(end+1) = ntraceraw(end);  % Check that ntraceraw is correct length
    end

    k=1;
    % Pre-allocate pulse variables
    Xd_pulse = cell(max(numavg),max(raw));
    Xo_pulse = cell(max(numavg),max(raw));
    Xc_pulse = cell(max(numavg),max(raw));
    kv_pulse = cell(max(numavg),max(raw));
    mv_pulse = cell(max(numavg),max(raw));
    Fe_pulse = cell(max(numavg),max(raw));
    gv_pulse = cell(max(numavg),max(raw));

    % Define pulse variables
    for j = raw
        for i = 1:numavg(k)
            xL2 = floor(length(Xd{j})/ntraceraw(k));
            for l = 1:ntraceraw(k)
                try
                    Xd_pulse{i,j}(:,l) = Xd{i,j}(1+(l-1)*xL2:l*xL2);
                catch
                end
                try
                    Xo_pulse{i,j}(:,l) = Xo{i,j}(1+(l-1)*xL2:l*xL2);
                catch
                end
                try
                    Xc_pulse{i,j}(:,l) = Xc{i,j}(1+(l-1)*xL2:l*xL2);
                catch
                end
                try
                    kv_pulse{i,j}(:,l) = kv{i,j}(1+(l-1)*xL2:l*xL2);
                catch
                end
                try
                    Fe_pulse{i,j}(:,l) = Fe{i,j}(1+(l-1)*xL2:l*xL2);
                catch
                end
                try
                    mv_pulse{i,j}(:,l) = mv{i,j}(1+(l-1)*xL2:l*xL2);
                catch
                end
                try
                    gv_pulse{i,j}(:,l) = gv{i,j}(1+(l-1)*xL2:l*xL2);
                catch
                end
            end
        end
        k=k+1;
    end
catch
    disp('Unable to import raw time series.');
end

% Extract non-raw pulses from full time traces
try
    k=1;
    while length(ntrace) < length(nonraw)
        ntrace(end+1) = ntrace(end);
    end

    % Pre-allocate averaged variables
    Xd_pulse_avg = cell(max(numavg(k),max(nonraw)));
    Xo_pulse_avg = cell(max(numavg(k),max(nonraw)));
    Xc_pulse_avg = cell(max(numavg(k),max(nonraw)));
    Fe_pulse_avg = cell(max(numavg(k),max(nonraw)));
    kv_pulse_avg = cell(max(numavg(k),max(nonraw)));
    mv_pulse_avg = cell(max(numavg(k),max(nonraw)));
    gv_pulse_avg = cell(max(numavg(k),max(nonraw)));

    % Define averaged variables
    for j = nonraw
        xL2 = length(data0{j})/ntrace(k);
        for i = 1:numavg(k)
            for l = 1:ntrace(k)
                try
                    Xd_pulse_avg{i,j}(:,l) = Xd{i,j}(1+(l-1)*xL2:l*xL2);
                catch
                end
                try
                    Xo_pulse_avg{i,j}(:,l) = Xo{i,j}(1+(l-1)*xL2:l*xL2);
                catch
                end
                try
                    Xc_pulse_avg{i,j}(:,l) = Xc{i,j}(1+(l-1)*xL2:l*xL2);
                catch
                end
                try
                    kv_pulse_avg{i,j}(:,l) = kv{i,j}(1+(l-1)*xL2:l*xL2);
                catch
                end
                try
                    Fe_pulse_avg{i,j}(:,l) = Fe{i,j}(1+(l-1)*xL2:l*xL2);
                catch
                end
                try
                    mv_pulse_avg{i,j}(:,l) = mv{i,j}(1+(l-1)*xL2:l*xL2);
                catch
                end
                try
                    gv_pulse_avg{i,j}(:,l) = gv{i,j}(1+(l-1)*xL2:l*xL2);
                catch
                end
            end
        end
        k=k+1;
    end
catch
    disp('Unable to import non-raw time series.');
end


% Time vector
dt = 1/Fs;
sizeX = size(Xd);
tvec = cell(sizeX(1),sizeX(2));
tvec_pulse = cell(sizeX(1),sizeX(2));
for j = 1:sizeX(1)
    for k = 1:sizeX(2)
        try
            tvec{j,k} = 0:dt:length(Xd{j,k})*dt-dt;
        catch
        end
        try
            tvec_pulse{j,k} = 0:dt:length(Xd_pulse{j,k})*dt-dt;
        catch
        end
    end
end

clear i j 

% Save extracted data
disp('Saving...');
savefile = sprintf('%s%s',datapath,'Extracted Data.mat');
save(savefile);
disp(['Saved as ' datapath 'Extracted Data.mat']);
disp('Finished.');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  getVLClog()  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = getVLClog(datapath,logdata,logfilepath)

% Log file path
if nargin < 3 || isempty(logfilepath)
    LFP = '/Users/joshsalvi/Documents/Lab/Lab/VLClogfilekey.xlsx';
else
    LFP = logfilepath;
end

% Import fieldnames
FNames = local_xls2fields(LFP,1,'A:EZ');
FNames = FNames(4:end);
FNames = unique(FNames,'stable');
Nfields = numel(FNames); % number of fields, some fields of arbitrary length
ilastvarlength = find(~cellfun(@isempty,strfind(FNames,'sizeOfStateSpaceFreqIncArray')), 1); % last index before varying length parameters

% Create file identifier for log file and open it
logfile = dir(sprintf('%s%s',datapath,'*.log')); % find the logfile
FN = fullfile(datapath,logfile.name);
fileID = fopen(FN);
if isempty(ilastvarlength)
    ilastvarlength = 5;
end
% No header line needs to be read in these log files.

% Generate an empty struct with field names as specified
S = cell2struct(repmat({[]},1,Nfields),FNames,2); S(:) = [];
fields = fieldnames(S);
for j = 1:size(logdata.data,1);
    for k = 1:numel(fields)
        S(j).(fields{k}) = logdata.data(j,k);
    end
end

% Close file
fclose(fileID);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  local_xls2fields()  %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function names = local_xls2fields(FN,sheet,range)
% LOCAL_XLS2FIELDS - read range of text from excel file and return fieldnames
%   names = local_xls2fields('foo.xls',sheet,range) reads excel file 'foo'
%   optional inputs sheet and range and converts any read text into valid
%   MATLAB fieldnames.
%
%   See also xlsread.

% Read column from excel file
[~, names] = xlsread(FN,sheet,range);

% Loop over all names and make them valid identifiers.
% Shamelessly copied and adapted from MATLAB 'built-in' function 'genvalidnames'.
if iscell(names), Nnames = numel(names); else Nnames = 1; end
for k = 1:Nnames
    name = names{k};
    if ~isvarname(name)
        
        % Replace # with 'number of'
        name = strrep(name, '#', 'number of');
        
        % Remove leading and trailing whitespace, and replace embedded
        % whitespace with camel/mixed casing.
        [~, afterSpace] = regexp(name,'\S\s+\S');
        if ~isempty(afterSpace)
            % Leave case alone except for chars after spaces
            name(afterSpace) = upper(name(afterSpace));
        end
        name = regexprep(name,'\s*','');
        if (isempty(name))
            name = 'x';
        end
        
        % Replace non-word character with underscore
        illegalChars = unique(name(regexp(name,'[^A-Za-z_0-9]')));
        for illegalChar = illegalChars
            replace = '_';
            name = strrep(name, illegalChar, replace);
        end
        
        % Remove double underscores
        name = regexprep(name,'__','_');
                
        % Insert x if the first column is non-letter.
        name = regexprep(name,'^\s*+([^A-Za-z])','x$1', 'once');
        
        % Prepend keyword with 'x' and camel case.
        if iskeyword(name)
            name = ['x' upper(name(1)) lower(name(2:end))];
        end
        
        % Truncate name to NAMLENGTHMAX
        name = name(1:min(length(name),namelengthmax));
        
        % Remove ending underscores
        name = fliplr(regexprep(fliplr(name),'_','', 'once'));
        names{k} = name;
    end
end

end

