function rt_sorted_stats(action, value)
% Record or flush rt_sorted lengths for offline analysis.
% Input:
%   action (1 x 1 char/string)
%       'record': record a new rt_sorted length value
%       'flush': flush the recorded values to a .mat file and reset
%   value (1 x 1 double, optional)
%       the rt_sorted length value to record when action is 'record';
%       the file path to save the recorded values when action is 'flush'

persistent rt_sorted_counts
if isempty(rt_sorted_counts)
    rt_sorted_counts = zeros(0, 1);
end

if nargin < 1 || isempty(action)
    error('rt_sorted_stats requires an action.');
end
if (nargin < 2 || isempty(value)) && ~strcmpi(action, 'init')
    error('rt_sorted_stats requires a value.');
end

switch lower(action)
    case 'init'
        rt_sorted_counts = zeros(0, 1);
        
    case 'record'
        rt_sorted_counts(end+1, 1) = value;

    case 'flush'
        if isempty(rt_sorted_counts)
            return;
        end
        stats_path = value;
        if exist(stats_path, 'file')
            stats_data = load(stats_path, 'rt_sorted_counts');
            if isfield(stats_data, 'rt_sorted_counts')
                rt_sorted_counts = [stats_data.rt_sorted_counts; rt_sorted_counts];
            end
        end
        save(stats_path, 'rt_sorted_counts', '-v7');
        rt_sorted_counts = zeros(0, 1);

    otherwise
        error('Unknown action: %s', action);
end
end
