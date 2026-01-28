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
num_iso = size(current_ratioMatrix,2);

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
    area_filter = zeros(num_iso,1);
    for idx_iso = 1:num_iso
        area_filter(idx_iso) = CChromatogramUtils.calculate_area(...
            rt_grid, intensityMatrix(:,idx_iso), ...
            curr_start, curr_end);
    end
    
    % Filter: keep only IMPs with area >= max_area * threshold
    max_area = max(area_filter);
    keep_mask = area_filter >= max_area * obj.m_resFilterThres;
    esti_ratio(curr_start:curr_end, ~keep_mask) = 0;

    % Uniform scaling: preserve total peak area
    sum_keep_area = sum(area_filter(keep_mask));
    if sum_keep_area > 0
        peak_total_area = sum(area_filter);  % Total area before filtering
        scale_factor = peak_total_area / sum_keep_area;
        esti_ratio(curr_start:curr_end, keep_mask) = ...
            esti_ratio(curr_start:curr_end, keep_mask) * scale_factor;
    end
end

% For each IMP, evaluate all candidate XIC peaks by computing:
%   - max_proportions: peak contribution ratio
%   - fwhm: peak width
%   - ratio_each_XIC_peak: area contribution in each peak

intensityMatrix = esti_ratio.*smoothed_intensity;
auxic = zeros(num_iso,1);
rt_bound = repmat(struct('start',0,'end',0), num_iso, length(XIC_peaks));
idx_selected = zeros(num_iso,1);
ratio_each_XIC_peak = zeros(num_iso,length(XIC_peaks));
for idx_iso = 1:num_iso
    max_proportions = zeros(length(XIC_peaks),1);
    fwhm = zeros(length(XIC_peaks),1);
    for i_Xp = 1:length(XIC_peaks)
        % Cache struct fields for performance and readability
        curr_start = XIC_peaks(i_Xp).left_bound;
        curr_end = XIC_peaks(i_Xp).right_bound;

        max_proportions(i_Xp) = max(esti_ratio(curr_start:curr_end,idx_iso));
        % Calculate the fwhm for XIC selection
        peak_rts = rt_grid(curr_start:curr_end);
        peak_intens = intensityMatrix(curr_start:curr_end,idx_iso);
        fwhm(i_Xp) = CChromatogramUtils.get_fwhm(peak_rts, peak_intens);
        % Calculate the ratio of each IMP in each XIC peak
        ratio_each_XIC_peak(idx_iso, i_Xp) = CChromatogramUtils.calculate_area(...
             rt_grid, intensityMatrix(:,idx_iso), ...
             curr_start, curr_end);
        
        % Store rt bounds in the same loop
        rt_bound(idx_iso,i_Xp).start = rt_grid(curr_start);
        rt_bound(idx_iso,i_Xp).end = rt_grid(curr_end);
    end
    % When there are more than one XIC peak counting for an IMP,
    % only record one XIC peak, the choosing critera are:
    % 1. the peak with max proportion, ratio shows the XIC belong more to this IMP;
    % 2. the peak with largest fwhm (full width of half maximum).
    [val_max, idx_selected(idx_iso)] = max(max_proportions);
    if length(find(abs(val_max-max_proportions)<eps)) > 1
        % If the proportion can not select the only one, choose according
        %   to the fwhm
        [~, idx_selected(idx_iso)] = max(fwhm);
    end
    
    idx_start = XIC_peaks(idx_selected(idx_iso)).left_bound;
    idx_end = XIC_peaks(idx_selected(idx_iso)).right_bound;
    
    auxic(idx_iso,1) = CChromatogramUtils.calculate_area(...
        rt_grid, intensityMatrix(:,idx_iso), idx_start, idx_end);
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