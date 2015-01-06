function S = getVLClog(datapath,logfilepath)
% getVLClog - load VLC log parameters
%   S = getVLClog(datapath) imports the log file for VLC data from folder
%   datapath. It returns struct S with fieldnames based on the excel file
%   'log file key.xlsx', located in the 'Data' folder on the recording
%   computer in the green room.
%   S = getVLClog(datapath,logfile) is similar, but now the logfile can be
%   loaded from any location specified by 'logfile'.
%
%   Current log file key used: log files from 2014-12-03 or after (Column D)
%
% Programmer:   Corstiaen Versteegh
% Date:         3 December 2014

% Log file path
if nargin < 2 || isempty(logfilepath)
    LFP = fullfile('/Volumes/Data','log file key.xlsx');
else
    LFP = logfilepath;
end

% Import fieldnames
FNames = local_xls2fields(LFP,1,'D:D');
FNames = FNames(2:end);
Nfields = numel(FNames); % number of fields, some fields of arbitrary length
ilastvarlength = find(~cellfun(@isempty,strfind(FNames,'sizeOfStateSpaceFreqIncArray'))); % last index before varying length parameters

% Create file identifier for log file and open it
logfile = dir(sprintf('%s%s',datapath,'*.log')); % find the logfile
FN = fullfile(datapath,logfile.name);
fileID = fopen(FN);

% No header line needs to be read in these log files.

% Generate an empty struct with field names as specified
S = cell2struct(repmat({[]},1,Nfields),FNames,2); S(:) = [];

% Run while loop "forever", endpoint is not known
while true
    startFormatSpec = ['%s %s %s' repmat(' %f',[1,ilastvarlength-3])];
    Data = textscan(fileID,startFormatSpec,1,'Delimiter','\t'); % start of dataline
    if isempty(Data{1}), break, end % if empty values are found, time to quit loop
    for ii = 9:-1:0 % hardcoded, needs updating when log file (key) is changed
        Dval = Data{ilastvarlength-ii};
        if Dval > 0
            loopFormatSpec = repmat('%f ',[1,Dval]); % loop over data of variable length
            DataLoop = textscan(fileID,loopFormatSpec,1,'Delimiter','\t'); % looping on data line
            Data{end+1} = cell2mat(DataLoop); %#ok<AGROW> % add looped data to existing data line
        else
            Data{end+1} = NaN; %#ok<AGROW>
        end
    end
    S(end+1) = cell2struct(Data,FNames,2); %#ok<AGROW> % assign data to struct
end

% Close file
fclose(fileID);



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




