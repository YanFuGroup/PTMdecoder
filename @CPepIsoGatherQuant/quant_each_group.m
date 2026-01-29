function [bhave_non_zeros, idxNonZero, auxic, rt_bound, idx_selected, ratio_each_XIC_peak]...
    = quant_each_group(obj,raw_name,current_ratioMatrix,current_rts,...
    current_inten, low_mz_bound, high_mz_bound, selected_charge)
% Quantify each group
% input:    
%   raw_name
%       the name of the raw (mgf) file
%   current_iso_name
%       names of IMPs in current group
%   current_ratioMatrix
%       ratio matrix of quantification in current group
%   current_rts
%       retention time in current group
%   current_inten
%       intensity in current group
%   low_mz_bound
%       low precursor m/z bound
%   high_mz_bound
%       high precursor m/z bound
%   selected_charge
%       current precursor charge
% output:   
%   bhave_non_zeros
%       is there non zero area under XIC
%   idxNonZero
%       the indices of non zero area under XIC
%   area
%       total quantification of each IMP in current group,
%       area under curve of XIC
%   rt_bound
%       the retention time bound, .start and .end
%   idx_selected
%       the indices of finally selected XIC peak

bhave_non_zeros = false;
idxNonZero = [];
auxic = [];
rt_bound = [];
idx_selected = [];
ratio_each_XIC_peak = [];
num_imp = size(current_ratioMatrix,2);

% Preprocess inputs (Sort, Smooth, Denoise)
[sort_rts, sort_inten, sort_ratioMatrix, is_valid] = ...
    CChromatogramUtils.preprocess_ms1_inputs(current_rts, current_inten, current_ratioMatrix, obj.m_minMSMSnum);

if ~is_valid
    bhave_non_zeros = false;
    idxNonZero = [];
    auxic = [];
    rt_bound = [];
    idx_selected = [];
    ratio_each_XIC_peak = [];
    return;
end

% Get Smoothed XIC
[rt_grid, smoothed_intensity, intensity] = ...
    CChromatogramUtils.get_smoothed_xic(obj.m_cMs12DatasetIO, raw_name, low_mz_bound, high_mz_bound, selected_charge);

% Extract the XIC peaks around the identified MSMS precursor
XIC_peaks = CChromatogramUtils.detect_xic_peaks(rt_grid, smoothed_intensity, intensity, sort_rts, obj.m_alpha);

if isempty(XIC_peaks)
    return;
end

% Calculate the ratio on each XIC points using kernel method
esti_ratio = CChromatogramUtils.calculate_kernel_ratio(rt_grid, sort_rts, sort_ratioMatrix, XIC_peaks, true);

% For each XIC peak, calculate the area of each IMP. If an IMP's area is
% less than a threshold (max_area * resFilterThres), remove it from that peak.
% Then, scale the remaining IMPs uniformly within the peak to preserve the
% total peak area (= sum of all pre-filter IMP areas in that peak).
intensityMatrix = esti_ratio.*smoothed_intensity;
for i_Xp = 1:length(XIC_peaks)
    curr_start = XIC_peaks(i_Xp).left_bound;
    curr_end = XIC_peaks(i_Xp).right_bound;

    % Calculate area for each IMP in this peak
    area_filter = zeros(num_imp,1);
    for idx_imp = 1:num_imp
        area_filter(idx_imp) = CChromatogramUtils.calculate_area(...
            rt_grid, intensityMatrix(:,idx_imp), ...
            curr_start, curr_end);
    end
    
    % Filter: keep only IMPs with area >= max_area * threshold
    max_area = max(area_filter);
    keep_mask = area_filter >= max_area * obj.m_resFilterThres;
    esti_ratio(curr_start:curr_end, ~keep_mask) = 0;

    % Normalize rows
    row_sum = sum(esti_ratio(curr_start:curr_end, :), 2);
    row_sum(row_sum == 0) = 1;
    esti_ratio(curr_start:curr_end, :) = esti_ratio(curr_start:curr_end, :) ./ repmat(row_sum, 1, num_imp);
end

% ========================================================================
% Stage 2: Feature Calculation (Vectorized over Peaks)
% ========================================================================

% For each IMP, evaluate all candidate XIC peaks by computing:
%   - imp_max_props, max peak contribution ratio
%   - peak_fwhms: half maximum peak width
%   - ratio_each_XIC_peak: area contribution in each peak

intensityMatrix = esti_ratio.*smoothed_intensity;
auxic = zeros(num_imp,1); 
% rt_bound will be expanded from single_rt_bounds later
idx_selected = zeros(num_imp,1);
ratio_each_XIC_peak = zeros(num_imp,length(XIC_peaks));

num_peaks = length(XIC_peaks);
peak_fwhms = zeros(num_imp, num_peaks);
single_rt_bounds = repmat(struct('start',0,'end',0), 1, num_peaks);
imp_max_props = zeros(num_imp, num_peaks);

for i_Xp = 1:num_peaks
    curr_start = XIC_peaks(i_Xp).left_bound;
    curr_end = XIC_peaks(i_Xp).right_bound;
    
    % 1. XIC-dependent properties (Independent of IMP)
    peak_rts = rt_grid(curr_start:curr_end);
    
    single_rt_bounds(i_Xp).start = rt_grid(curr_start);
    single_rt_bounds(i_Xp).end = rt_grid(curr_end);
    
    % 2. IMP-dependent properties (Vectorized logic)
    ratio_slice = esti_ratio(curr_start:curr_end, :);
    imp_max_props(:, i_Xp) = max(ratio_slice, [], 1)';
    
    % Calculate area for record
    for idx_imp = 1:num_imp
        % Calculate FWHM using IMP-specific intensity
        peak_fwhms(idx_imp, i_Xp) = CChromatogramUtils.get_fwhm(peak_rts, intensityMatrix(curr_start:curr_end, idx_imp));
        
        ratio_each_XIC_peak(idx_imp, i_Xp) = CChromatogramUtils.calculate_area(...
             rt_grid, intensityMatrix(:,idx_imp), curr_start, curr_end);
    end
end

rt_bound = repmat(single_rt_bounds, num_imp, 1);

% ========================================================================
% Stage 3: Peak Selection (Per IMP)
% ========================================================================

for idx_imp = 1:num_imp
    % Get features for this IMP across all peaks
    props = imp_max_props(idx_imp, :);

    % 1. Find all peaks tied for the maximum value (within epsilon tolerance)
    %    Since props cannot exceed max(props), we check props >= max - eps
    candidates = find(props >= max(props) - eps);
    
    % 2. Select the candidate with the largest FWHM
    [~, best_subidx] = max(peak_fwhms(idx_imp, candidates));
    
    idx_selected(idx_imp) = candidates(best_subidx);
end

% ========================================================================
% Stage 4: Global Refinement (Re-distribution based on Selection)
% ========================================================================
% Construct a mask to keep only the valid regions for each IMP.
keep_mask = false(size(esti_ratio));
for idx_imp = 1:num_imp
    % Retrieve the selected peak bounds
    sel_peak_idx = idx_selected(idx_imp);
    
    if sel_peak_idx > 0 && sel_peak_idx <= length(XIC_peaks)
        p_start = XIC_peaks(sel_peak_idx).left_bound;
        p_end   = XIC_peaks(sel_peak_idx).right_bound;
        
        % Mark valid region
        keep_mask(p_start:p_end, idx_imp) = true;
    end
end

% Apply mask: zero out regions not selected by the IMP
esti_ratio(~keep_mask) = 0;

% Re-calculate intensity matrix with the refined ratios
intensityMatrix = esti_ratio .* smoothed_intensity;

% ========================================================================
% Stage 5: Final Area Calculation
% ========================================================================

auxic = zeros(num_imp,1);
for idx_imp = 1:num_imp
    % Retrieve the selected peak bounds
    sel_peak_idx = idx_selected(idx_imp);
    idx_start = XIC_peaks(sel_peak_idx).left_bound;
    idx_end = XIC_peaks(sel_peak_idx).right_bound;
    
    % Calculate final area using the refined intensity matrix
    auxic(idx_imp,1) = CChromatogramUtils.calculate_area(...
        rt_grid, intensityMatrix(:,idx_imp), idx_start, idx_end);
end

% Get the non-zero area under XIC, index and rt_bound
idxNonZero = find(auxic(:,1)~=0);
auxic = auxic(idxNonZero,:);
rt_bound = rt_bound(idxNonZero,:);
idx_selected = idx_selected(idxNonZero,:);
% TODO: use auxic to calculate the ratio_each_XIC_peak instead.
ratio_each_XIC_peak = ratio_each_XIC_peak(idxNonZero,:);
ratio_each_XIC_peak = ratio_each_XIC_peak./sum(ratio_each_XIC_peak,1);
if ~isempty(idxNonZero)
    bhave_non_zeros = true;
end
end