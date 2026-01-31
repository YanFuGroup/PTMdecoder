function [result] = ReadDatResultFolder(path)
% Read all Mascot .dat files in a folder
% Input:
%   path (1 x 1 char/string)
%       folder path containing .dat files
% Output:
%   result (1 x N struct)
%       concatenated results

result = [];
datafiles = dir(fullfile(path,'*.dat'));
nfiles = length(datafiles);
for i=1:nfiles
    plfilepath = fullfile(path,datafiles(i).name);
    fprintf('Reading %s...',plfilepath);
    [resulti] = ReadDatResult(plfilepath);   
    result = [result resulti];
    fprintf('\n')   
end

end