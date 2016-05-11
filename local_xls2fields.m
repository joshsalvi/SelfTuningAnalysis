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
