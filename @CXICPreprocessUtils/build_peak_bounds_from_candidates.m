function [xic_peak_idx_bounds, xic_peak_rt_bounds, is_ok] = build_peak_bounds_from_candidates(xic_rt, candidate_rt_peaks)
% Map candidate RT intervals to index bounds on xic_rt.
% Input:
%   xic_rt (N x 1 double)
%       XIC retention-time grid
%   candidate_rt_peaks (struct array)
%       Candidate intervals with fields rt_start/rt_end
% Output:
%   xic_peak_idx_bounds (struct array)
%       Index bounds (idx_start/idx_end)
%   xic_peak_rt_bounds (struct array)
%       Sanitized RT bounds corresponding to index bounds
%   is_ok (logical)
%       True when at least one valid interval is projected

xic_peak_idx_bounds = repmat(struct('idx_start', 0, 'idx_end', 0), 0, 1);
xic_peak_rt_bounds = repmat(struct('rt_start', 0, 'rt_end', 0), 0, 1);

for idx_peak = 1:numel(candidate_rt_peaks)
    peak = candidate_rt_peaks(idx_peak);
    if ~isfield(peak, 'rt_start') || ~isfield(peak, 'rt_end')
        continue;
    end

    rt_start = peak.rt_start;
    rt_end = peak.rt_end;
    if isempty(rt_start) || isempty(rt_end) || rt_end < rt_start
        continue;
    end

    [~, idx_start] = min(abs(xic_rt - rt_start));
    [~, idx_end] = min(abs(xic_rt - rt_end));
    if idx_end < idx_start
        idx_end = idx_start;
    end

    xic_peak_idx_bounds(end + 1, 1) = struct('idx_start', idx_start, 'idx_end', idx_end); %#ok<AGROW>
    xic_peak_rt_bounds(end + 1, 1) = struct( ...
        'rt_start', xic_rt(idx_start), ...
        'rt_end', xic_rt(idx_end)); %#ok<AGROW>
end

is_ok = ~isempty(xic_peak_idx_bounds);
end