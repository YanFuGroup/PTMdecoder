function task_param_map = parseFileToMap(param_file_path)
% Parse parameter file into a key-value map.
%   param_file_path (1 x 1 char/string)
%       path to the parameter file
fid = fopen(param_file_path, 'r');
if fid <= 0
    error('CWorkflowParamParser:OpenParamFileFailed', ...
        'Failed to open parameter file: %s', param_file_path);
end

task_param_map = containers.Map();
line_num = 0;
while ~feof(fid)
    line_num = line_num + 1;
    str_line = fgetl(fid);
    if ~ischar(str_line)
        continue;
    end

    str_line = CWorkflowParamParser.removeComments(str_line);
    str_line = strtrim(str_line);
    if isempty(str_line)
        continue;
    end

    str_seg = split(str_line, '=');
    if length(str_seg) ~= 2
        error('CWorkflowParamParser:InvalidParamFormat', ...
            'Unexpected parameter format in line %d: %s', line_num, str_line);
    end

    m_key = strtrim(str_seg(1));
    m_val = strtrim(str_seg(2));
    task_param_map(m_key{1}) = m_val{1};
end
fclose(fid);
end