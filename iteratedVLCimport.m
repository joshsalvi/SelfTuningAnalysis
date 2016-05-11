function iteratedVLCimport(varargin)
% This function calls importVLCdata5() to iterate through all folders in a
% tree for that particular date. Each folder's data are be stored within
% that folder in a MAT file (Extracted Data.mat).
%
% iteratedVLCimport(datapath, variableorder)
%
% INPUTS:
%   datapath:      string containing the data path to be searched
%       -> Data path must be the DATE'S folder (see EXAMPLE)
%   variableorder: the order of variables (Xd,Xo,Xc,kv,Fe,mv,gv) for import
%       -> Type "help importVLCdata5" for more information
%       -> If no values are included, the order will default to "1:7"
%
% EXAMPLE:
%   iteratedVLCimport('/Users/joshsalvi/Documents/Lab/Lab/Clamp
%   Data/2016-05-02.01/')
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

folder1 = dir(datapath);
for j = 1:length(folder1)
    if folder1(j).name ~= '.'
        datapath1 = [datapath folder1(j).name '/'];
        folder2 = dir(datapath1);
        for k = 1:length(folder2)
            if folder2(k).name ~= '.'
                datapath2 = [datapath1 folder2(k).name '/'];
                FL = folder2(k).name;
                if FL(1) ~= 'F'
                    try
                        disp(['Importing ' datapath2 ' ...']);
                        importVLCdata5(datapath2,varsincluded);
                    catch
                        disp(['Unable to import ' datapath2 ' ...'])
                        disp('Skipped.');
                    end
                end
            end
        end
    end
end


end
