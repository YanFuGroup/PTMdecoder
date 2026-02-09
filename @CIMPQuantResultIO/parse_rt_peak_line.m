function peak = parse_rt_peak_line(strline)
% Parse an RT peak line starting with '@'
% Line format: @\trt_start\trt_end\tratio\tcheck_label
% Exactly 5 tab-separated fields, starting with '@', all fields 2-5 must be numeric
parts = regexp(strline, '\t', 'split');
if numel(parts) ~= 5 || ~strcmp(parts{1}, '@')
    error(['The line: "', strline, '" representing RT peaks is in an unexpected format!']);
end

numbers = str2double(parts(2:5));
if any(isnan(numbers))
    error(['The line: "', strline, '" representing RT peaks is in an unexpected format!']);
end

peak = struct(...
    'rt_start', numbers(1), ...
    'rt_end', numbers(2), ...
    'ratio', numbers(3), ...
    'check_label', numbers(4));
end
