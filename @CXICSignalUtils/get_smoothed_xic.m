function [xic_rt, xic_intensity_smoothed, xic_intensity_raw] = get_smoothed_xic(ms12DatasetIO, raw_name, low_mz_bound, high_mz_bound, selected_charge)
% Extracts and smooths the Extracted Ion Chromatogram (XIC) from MS1 data.
%
% Inputs:
%   ms12DatasetIO (object)
%       Instance of CMS12DatasetIO containing parsed MS data.
%   raw_name (1 x 1 char/string)
%       The raw file name (e.g., 'sample.mgf').
%   low_mz_bound (1 x 1 double)
%       Lower m/z bound for the monoisotopic peak.
%   high_mz_bound (1 x 1 double)
%       Upper m/z bound for the monoisotopic peak.
%   selected_charge (1 x 1 double/int)
%       Charge state of the precursor ion.
%
% Output:
%   xic_rt (N x 1 double)
%       Retention time grid in minutes.
%   xic_intensity_smoothed (N x 1 double)
%       Smoothed XIC intensity of the monoisotopic peak.
%   xic_intensity_raw (N x 1 double)
%       Raw XIC intensity of the monoisotopic peak.

% -------------------------------------------------------------------------
% 1. Data Initialization & Retrieval
% -------------------------------------------------------------------------
mgf_stem = erase(raw_name, '.mgf');
ms1_stem = ms12DatasetIO.m_cMsFileMapper.get_ms1_stem(mgf_stem);

% MS1_index columns: [scan, retention time, peak number, baseline, injection time]
% MS1_peaks columns: [m/z, intensity]
MS1_index = ms12DatasetIO.m_mapNameMS1Index(ms1_stem);
MS1_peaks = ms12DatasetIO.m_mapNameMS1Peaks(ms1_stem);

xic_rt = MS1_index(:,2);
isotope_num = [-1, 0, 1, 2, 3, 4];
num_scans = size(MS1_index, 1);
num_iso = length(isotope_num);

% Initialize the raw XIC matrix (rows: scans, cols: isotopes)
xic_intensity_raw = zeros(num_scans, num_iso); 

% -------------------------------------------------------------------------
% 2. Pre-computation for Peak Extraction
% -------------------------------------------------------------------------
if isstruct(ms12DatasetIO)
    has_sorted_cache = isfield(ms12DatasetIO, 'm_mapNameMS1SortedMz') && ...
        isfield(ms12DatasetIO, 'm_mapNameMS1SortedInt') && ...
        isfield(ms12DatasetIO, 'm_mapNameMS1SortedScan') && ...
        isKey(ms12DatasetIO.m_mapNameMS1SortedMz, ms1_stem) && ...
        isKey(ms12DatasetIO.m_mapNameMS1SortedInt, ms1_stem) && ...
        isKey(ms12DatasetIO.m_mapNameMS1SortedScan, ms1_stem);
else
    has_sorted_cache = isprop(ms12DatasetIO, 'm_mapNameMS1SortedMz') && ...
        isprop(ms12DatasetIO, 'm_mapNameMS1SortedInt') && ...
        isprop(ms12DatasetIO, 'm_mapNameMS1SortedScan') && ...
        isKey(ms12DatasetIO.m_mapNameMS1SortedMz, ms1_stem) && ...
        isKey(ms12DatasetIO.m_mapNameMS1SortedInt, ms1_stem) && ...
        isKey(ms12DatasetIO.m_mapNameMS1SortedScan, ms1_stem);
end

if has_sorted_cache
    sorted_mz = ms12DatasetIO.m_mapNameMS1SortedMz(ms1_stem);
    sorted_intensity = ms12DatasetIO.m_mapNameMS1SortedInt(ms1_stem);
    sorted_scan = ms12DatasetIO.m_mapNameMS1SortedScan(ms1_stem);
else
    mz_vals = MS1_peaks(:, 1);
    intensity_vals = MS1_peaks(:, 2);
    scan_boundaries = [0; MS1_index(:, 3); inf];
    peak_global_idx = (1:size(MS1_peaks, 1))';
    scan_idx_all = discretize(peak_global_idx, scan_boundaries);

    [sorted_mz, sort_order] = sort(mz_vals, 'ascend');
    sorted_intensity = intensity_vals(sort_order);
    sorted_scan = scan_idx_all(sort_order);
end

% Pre-compute m/z extraction windows for all target isotopes
mz_shift = reshape(isotope_num * (CConstant.unitdiff / selected_charge), 1, []);
L_bounds = low_mz_bound + mz_shift;
H_bounds = high_mz_bound + mz_shift;

% -------------------------------------------------------------------------
% 3. Vectorized XIC Extraction
% -------------------------------------------------------------------------
% Run one binary-search pass per isotope window on globally sorted m/z.
for idx_iso = 1:num_iso
    idx_start = binary_search_first_gt(sorted_mz, L_bounds(idx_iso));
    idx_end = binary_search_last_lt(sorted_mz, H_bounds(idx_iso));
    if idx_start > idx_end
        continue;
    end

    target_scan = sorted_scan(idx_start:idx_end);
    target_intensity = sorted_intensity(idx_start:idx_end);

    % Keep the last hit in sorted-window order for each scan.
    [scan_unique, idx_last] = unique(target_scan, 'last');
    intensity_last = target_intensity(idx_last);
    xic_intensity_raw(scan_unique, idx_iso) = intensity_last;
end

% -------------------------------------------------------------------------
% 4. Quality Control & Filtering
% -------------------------------------------------------------------------
% Filter out scans that do not meet expected isotopic distribution criteria.

% Criterion 1: The M-1 peak intensity must not exceed the monoisotopic peak
invalid_mask = xic_intensity_raw(:, 1) > xic_intensity_raw(:, 2);

% Isolate the primary isotopic envelope (M0 to M4)
obs_matrix = xic_intensity_raw(:, 2:end);
max_obs = max(obs_matrix, [], 2);

% Identify valid rows to compute cosine similarity (non-zero and pass Criterion 1)
valid_rows = (max_obs > 0) & ~invalid_mask;

if any(valid_rows)
    % Normalize observed intensities by the maximum intensity in the envelope
    obs_norm = obs_matrix(valid_rows, :) ./ max_obs(valid_rows);
    
    % Retrieve and format the theoretical Isotope Pattern Vector (IPV)
    ipv_ref = double(CConstant.IPV(int64((high_mz_bound + low_mz_bound) / 2), :));
    
    % Criterion 2: Vectorized Cosine Similarity check against theoretical IPV
    % CosSim = (A dot B) / (||A|| * ||B||)
    dot_prod = obs_norm * ipv_ref';          
    norm_ipv = norm(ipv_ref);
    norm_obs = sqrt(sum(obs_norm.^2, 2));    
    
    cos_sim = dot_prod ./ (norm_ipv .* norm_obs);
    
    % Mark rows with cosine similarity < 0.6 as invalid
    failed_sim_idx = cos_sim < 0.6;
    valid_indices = find(valid_rows);
    invalid_mask(valid_indices(failed_sim_idx)) = true;
end

% Zero out all intensities for scans that failed the quality filters
xic_intensity_raw(invalid_mask, :) = 0;

% -------------------------------------------------------------------------
% 5. Output Formatting & Smoothing
% -------------------------------------------------------------------------
% Extract the monoisotopic trace (M0 is now the 2nd column)
xic_intensity_raw = xic_intensity_raw(:, 2);

% Apply moving average smoothing
xic_intensity_smoothed = smoothdata(xic_intensity_raw, 'movmean', 5);

end


function idx = binary_search_first_gt(sorted_vals, target)
% Locate the first index with value strictly greater than target.

low = 1;
high = length(sorted_vals);
idx = high + 1;

while low <= high
    mid = floor((low + high) / 2);
    if sorted_vals(mid) > target
        idx = mid;
        high = mid - 1;
    else
        low = mid + 1;
    end
end
end


function idx = binary_search_last_lt(sorted_vals, target)
% Locate the last index with value strictly smaller than target.

low = 1;
high = length(sorted_vals);
idx = 0;

while low <= high
    mid = floor((low + high) / 2);
    if sorted_vals(mid) < target
        idx = mid;
        low = mid + 1;
    else
        high = mid - 1;
    end
end
end