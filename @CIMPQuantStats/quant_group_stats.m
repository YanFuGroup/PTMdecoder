function stats = quant_group_stats(action, value)
% Record and summarize initial quantification group outcomes.

persistent counts
persistent warned_keys
if isempty(counts)
    counts = make_empty_counts(5);
end
if isempty(warned_keys)
    warned_keys = {};
end

switch lower(action)
    case 'init'
        counts = make_empty_counts(value);
        warned_keys = {};
    case 'record'
        diagnostics = value;
        counts.total_groups = counts.total_groups + 1;

        [reason, warning_key, warning_message] = parse_reason(diagnostics);
        if ~isempty(warning_key)
            warned_keys = warn_once(warned_keys, warning_key, warning_message);
        end
        counts.(reason) = counts.(reason) + 1;

        [candidate_count, filtered_count, warning_key, warning_message] = ...
            parse_sparse_metrics(diagnostics);
        if ~isempty(warning_key)
            warned_keys = warn_once(warned_keys, warning_key, warning_message);
        end
        counts.candidate_peaks = counts.candidate_peaks + candidate_count;
        counts.filtered_sparse_peaks = counts.filtered_sparse_peaks + filtered_count;
        if filtered_count > 0
            counts.sparse_filtered_groups = counts.sparse_filtered_groups + 1;
        end
    case 'get'
        % Return current counters without resetting them.
    otherwise
        error('CIMPQuantStats:UnknownQuantGroupStatsAction', ...
            'Unknown quant_group_stats action: %s', action);
end
stats = counts;
end

function counts = make_empty_counts(min_nonzero_points)
if ~isnumeric(min_nonzero_points) || ~isscalar(min_nonzero_points) || ...
        ~isfinite(min_nonzero_points) || min_nonzero_points < 1 || ...
        floor(min_nonzero_points) ~= min_nonzero_points
    error('CIMPQuantStats:InvalidMinXicNonzeroPoints', ...
        'min_xic_nonzero_points must be a positive integer.');
end
counts = struct( ...
    'min_xic_nonzero_points', min_nonzero_points, ...
    'total_groups', 0, ...
    'success', 0, ...
    'insufficient_psm_inputs', 0, ...
    'no_xic_peak', 0, ...
    'sparse_xic_peaks', 0, ...
    'zero_imp_area', 0, ...
    'unknown', 0, ...
    'candidate_peaks', 0, ...
    'sparse_filtered_groups', 0, ...
    'filtered_sparse_peaks', 0);
end

function [reason, warning_key, warning_message] = parse_reason(diagnostics)
known_reasons = {'success', 'insufficient_psm_inputs', 'no_xic_peak', ...
    'sparse_xic_peaks', 'zero_imp_area', 'unknown'};
reason = 'unknown';
warning_key = '';
warning_message = '';

if ~isstruct(diagnostics) || ~isscalar(diagnostics) || ~isfield(diagnostics, 'reason') || ...
        ~((ischar(diagnostics.reason) && isrow(diagnostics.reason)) || ...
        (isstring(diagnostics.reason) && isscalar(diagnostics.reason)))
    warning_key = 'invalid_reason';
    warning_message = ['[CIMPQuantStats:InvalidDiagnosticsReason] diagnostics.reason ', ...
        'must be a char or string scalar; recording as unknown.'];
    return;
end

raw_reason = char(string(diagnostics.reason));
if any(strcmp(raw_reason, known_reasons))
    reason = raw_reason;
else
    warning_key = ['unknown_reason:', raw_reason];
    warning_message = sprintf(['[CIMPQuantStats:UnknownDiagnosticsReason] Unknown quant ', ...
        'group reason "%s"; recording as unknown.'], raw_reason);
end
end

function [candidate_count, filtered_count, warning_key, warning_message] = parse_sparse_metrics(diagnostics)
candidate_count = 0;
filtered_count = 0;
warning_key = '';
warning_message = '';

if ~isstruct(diagnostics) || ~isscalar(diagnostics) || ...
        ~isfield(diagnostics, 'candidate_peak_count') || ...
        ~isfield(diagnostics, 'filtered_sparse_peak_count')
    warning_key = 'missing_sparse_metrics';
    warning_message = ['[CIMPQuantStats:MissingSparseMetrics] diagnostics must contain ', ...
        'candidate_peak_count and filtered_sparse_peak_count; recording sparse metrics as zero.'];
    return;
end

candidate_value = diagnostics.candidate_peak_count;
filtered_value = diagnostics.filtered_sparse_peak_count;
is_valid_count = @(x) isnumeric(x) && isscalar(x) && isfinite(x) && ...
    x >= 0 && floor(x) == x;
if ~is_valid_count(candidate_value) || ~is_valid_count(filtered_value) || ...
        filtered_value > candidate_value
    warning_key = 'invalid_sparse_metrics';
    warning_message = ['[CIMPQuantStats:InvalidSparseMetrics] Sparse diagnostic counts ', ...
        'must be non-negative integers with filtered <= candidate; recording sparse metrics as zero.'];
    return;
end

candidate_count = candidate_value;
filtered_count = filtered_value;
end

function warned_keys = warn_once(warned_keys, warning_key, warning_message)
if ~any(strcmp(warned_keys, warning_key))
    CLogger.warn('%s', warning_message);
    warned_keys{end + 1} = warning_key;
end
end
