function task_param_map = parse_file(~, param_file)
% Parse the task parameter file and return the dictionary of user settings.
% Input:
%   param_file
%       the file path of the task parameter file.
% Output:
%   task_param_map
%       the dictionary of user settings.

fid = fopen(param_file,'r');
if 0 >= fid
    error('Failed to open the task parameter file!');
end

task_param_map = containers.Map();
fileLine = 0;
while ~feof(fid)
    fileLine = fileLine + 1;
    strLine = fgetl(fid);

    % Remove comments and empty lines
    strLine = remove_comments(strLine);
    if isempty(strLine)
        continue;
    end

    % Record user settings to the dictionary
    task_param_map = record_user_settings(strLine, task_param_map, fileLine);
end
fclose(fid);
end



function strLine = remove_comments(strLine)
% Remove comments in the line.
% Input:
%   strLine
%       the line to be processed.
% Output:
%   strLine
%       the line without comments.

if contains(strLine,'#')
    iStr = strfind(strLine,'#');
    strLine(iStr(1):end) = [];
end
end



function task_param_map = record_user_settings(strLine, task_param_map, fileLine)
% Record user settings to the dictionary in current line.
% Input:
%   strLine
%       the line to be processed.
%   task_param_map
%       the dictionary of user settings.
%   fileLine
%       the line number in the file.
% Output:
%   task_param_map
%       the dictionary of user settings.

strSeg = split(strLine,'=');
if length(strSeg) > 2
    error(['Unexpected parameter format in line ',fileLine,'!']);
end
mKey = strtrim(strSeg(1));
mValue = strtrim(strSeg(2));
task_param_map(mKey{1}) = mValue{1};
end