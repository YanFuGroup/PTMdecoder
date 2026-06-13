function stats = quant_group_stats(action, value)
% Record and summarize initial quantification group outcomes.

persistent counts
if isempty(counts)
    counts = make_empty_counts();
end

switch lower(action)
    case 'init'
        counts = make_empty_counts();
    case 'record'
        field_name = char(string(value));
        if ~isfield(counts, field_name)
            counts.(field_name) = 0;
        end
        counts.(field_name) = counts.(field_name) + 1;
    case 'get'
        % Return current counters without resetting them.
    otherwise
        error('CIMPQuantStats:UnknownQuantGroupStatsAction', ...
            'Unknown quant_group_stats action: %s', action);
end
stats = counts;
end

function counts = make_empty_counts()
counts = struct( ...
    'success', 0, ...
    'insufficient_psm_inputs', 0, ...
    'no_xic_peak', 0, ...
    'sparse_xic_peaks', 0, ...
    'zero_imp_area', 0, ...
    'unknown', 0);
end
