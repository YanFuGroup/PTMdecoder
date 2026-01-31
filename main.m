% SPDX-License-Identifier: BSD-3-Clause-Clear

function main(varargin)
% main function of PTMdecoder
% Input:
%   varargin (1 x N cell):
%       paths of parameter files

if nargin == 0
    error(sprintf(['Please input the parameter file.\n' ...
               'Syntax: PTMdecoder_core.exe <parameter files>\n' ...
               'Example: PTMdecoder_core.exe task1.param task2.param\n' ...
               'Please refer to the user manual for details on setting up the parameter files.']));
else
    % Although there is a member function in CPTMdecoder to check parameter files, 
    % it is necessary to check the parameter files before running the program with all parameter files.
    check_param_files(varargin);
    CPTMdecoder(varargin{:});
end

end



function check_param_files(param_files)
% Check whether the parameter files exist.
% Input:
%   param_files (1 x N cell)
%       the paths of parameter files.

for i = 1:length(param_files)
    if ~fileExists(param_files{i})
        error(['The parameter file ' param_files{i} ' does not exist!']);
    end
end
end



function exists = fileExists(fileName)
% Check whether the file exists.
% Input:
%   fileName (1 x N char/string)
% Output:
%   exists (1 x 1 logical)
fid = fopen(fileName, 'r');
if fid == -1
    exists = false;
else
    exists = true;
    fclose(fid);
end
end