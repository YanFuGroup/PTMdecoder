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
% Extract from left to right, each peak is record by 
% [left_bound, right_bound]. Finally record the peaks in cell.
idx_PSM = 1;
XIC_peaks = struct('left_bound',{},'right_bound',{});
i_Xp = 1;
while idx_PSM <= length(sort_rts)
    % index of rt for first identified MS/MS in this peak
    [~, idx_first_rt] = min(abs(rt_grid-sort_rts(idx_PSM)));
    max_peak_inten = smoothed_intensity(idx_first_rt);
    if max_peak_inten == 0
        % An identification with intensity of zero means the a
        %   mis-identification, which can not pass the isotopic filtering.
        %   So these identifications should be skipped.
        idx_PSM = idx_PSM + 1;
        continue;
    end

    % look after the right boundary
    XIC_peaks(i_Xp).right_bound = idx_first_rt;
    min_peak_inten = smoothed_intensity(idx_first_rt);
    min_peak_iter = idx_first_rt;
    for iter_rt = idx_first_rt+1:length(rt_grid)
        % find check whether reach the next rt of identified MS/MS precursor
        if idx_PSM<length(sort_rts) && rt_grid(iter_rt)>sort_rts(idx_PSM+1)
            XIC_peaks(i_Xp).right_bound = iter_rt;
            idx_PSM = idx_PSM + 1;
        end
        % update the local maximum
        if max_peak_inten < smoothed_intensity(iter_rt)
            max_peak_inten = smoothed_intensity(iter_rt);
            % update the min peak to ensure that the local minima is at
            %   right of the left maxima
            min_peak_inten = max_peak_inten;
            min_peak_iter = iter_rt;
        end
        % find the local minimum
        if smoothed_intensity(iter_rt) < min_peak_inten
            min_peak_inten = smoothed_intensity(iter_rt);
            min_peak_iter = iter_rt;
        end
        % find the right bound
        % two criteria: 1. the right is too low; 2. the local minimum is too low
        if smoothed_intensity(iter_rt) < max_peak_inten*obj.m_alpha
            XIC_peaks(i_Xp).right_bound = iter_rt - 1;
            break;
        elseif min_peak_inten < 0.5 * min(smoothed_intensity(iter_rt),max_peak_inten)
            XIC_peaks(i_Xp).right_bound = min_peak_iter;
            break;
        end
    end

    % look after the left boundary
    XIC_peaks(i_Xp).left_bound = idx_first_rt;
    min_peak_inten = smoothed_intensity(idx_first_rt);
    min_peak_iter = idx_first_rt;
    for iter_rt = idx_first_rt-1:-1:1
        % update the local maximum
        if max_peak_inten < smoothed_intensity(iter_rt)
            max_peak_inten = smoothed_intensity(iter_rt);
            % update the min peak to ensure that the local minima is at
            %   right of the left maxima
            min_peak_inten = max_peak_inten;
            min_peak_iter = iter_rt;
        end
        % find the local minimum
        if smoothed_intensity(iter_rt) < min_peak_inten
            min_peak_inten = smoothed_intensity(iter_rt);
            min_peak_iter = iter_rt;
        end
        % find the left bound
        % two criteria: 1. the left is too low; 2. the minimum is too low
        if smoothed_intensity(iter_rt) < max_peak_inten*obj.m_alpha
            XIC_peaks(i_Xp).left_bound = iter_rt + 1;
            break;
        elseif min_peak_inten < 0.5 * min(smoothed_intensity(iter_rt),max_peak_inten)
            XIC_peaks(i_Xp).left_bound = min_peak_iter;
            break;
        end
    end

    % prepare for next peak extraction
    i_Xp = i_Xp + 1;
    idx_PSM = idx_PSM + 1;
end

% remove the XIC peaks only with less than 5 point
for i_Xp = length(XIC_peaks):-1:1
    if sum(intensity(XIC_peaks(i_Xp).left_bound:XIC_peaks(i_Xp).right_bound)~=0) < 5
        XIC_peaks(i_Xp) = [];
    end
end
if isempty(XIC_peaks)
    return;
end

% Calculate the ratio on each XIC points using kernel method, and normalize
% using Nadaraya-Waston kernel averaging method
esti_ratio = zeros(length(rt_grid),num_iso);
for idx_iso = 1:num_iso
    % Gaussian kernel function
    Ker_Gaussian = @(u) (1/sqrt(2*pi))*exp(-0.5*u.^2);
    % Epanechnikov kernel function
%     Ker_Epanechnikov = @(u) (3/4)*(1-u.^2).*(abs(u)<=1);
    % set kernel weights
    for i_Xp = 1:length(XIC_peaks)
        % Collect all of the rts states within current XIC peak
        idxs_rt = sort_rts>=rt_grid(XIC_peaks(i_Xp).left_bound) &...
            sort_rts<=rt_grid(XIC_peaks(i_Xp).right_bound);
        rts_current = sort_rts(idxs_rt);
        % Calculate bandwidth and weights for each XIC peak
        bandwidth = (4/(3*size(rts_current,1)))^0.2*std(rts_current);
        weights = zeros(length(rt_grid), length(rts_current));
        for idx_PSM = 1:length(rts_current)
            if bandwidth == 0
                break;
            end
            weights(XIC_peaks(i_Xp).left_bound:XIC_peaks(i_Xp).right_bound,idx_PSM) = ...
                Ker_Gaussian((rt_grid(XIC_peaks(i_Xp).left_bound:XIC_peaks(i_Xp).right_bound)...
                -rts_current(idx_PSM))/bandwidth);
        end
        % Check if there are nearly no weights in some retention time for all
        %   IMP, or the bandwidth is just zero
        %         if ~all(weights(XIC_peaks(i_Xp).left_bound:XIC_peaks(i_Xp).right_bound,:))
        if bandwidth == 0 ||...
                any(all(weights(XIC_peaks(i_Xp).left_bound:XIC_peaks(i_Xp).right_bound,:)<1e-15/length(sort_rts),2))
            bandwidth = min(rt_grid(XIC_peaks(i_Xp).right_bound)-rt_grid(XIC_peaks(i_Xp).left_bound),1);
            for idx_PSM = 1:length(rts_current)
                weights(XIC_peaks(i_Xp).left_bound:XIC_peaks(i_Xp).right_bound,idx_PSM) = ...
                    Ker_Gaussian((rt_grid(XIC_peaks(i_Xp).left_bound:XIC_peaks(i_Xp).right_bound)...
                    -rts_current(idx_PSM))/bandwidth);
            end
        end
        % Calculate the ratio using normalized weights and ratioMatrix
        esti_ratio(XIC_peaks(i_Xp).left_bound:XIC_peaks(i_Xp).right_bound,idx_iso) = ...
            esti_ratio(XIC_peaks(i_Xp).left_bound:XIC_peaks(i_Xp).right_bound,idx_iso) + ...
            (weights(XIC_peaks(i_Xp).left_bound:XIC_peaks(i_Xp).right_bound,:) * ...
            sort_ratioMatrix(idxs_rt,idx_iso))./ ...
            (sum(weights(XIC_peaks(i_Xp).left_bound:XIC_peaks(i_Xp).right_bound,:),2)+eps);
    end
end
% normalize the ratio in every available retention time
esti_ratio = esti_ratio./(sum(esti_ratio,2)+eps);

% Filter low abundance IMP according to the relative AUXIC
intensityMatrix = esti_ratio.*smoothed_intensity;
for i_Xp = 1:length(XIC_peaks)
    area_filter = zeros(num_iso,1);
    for idx_iso = 1:num_iso
        area_filter(idx_iso) = trapz(rt_grid(XIC_peaks(i_Xp).left_bound:XIC_peaks(i_Xp).right_bound),...
            intensityMatrix(XIC_peaks(i_Xp).left_bound:XIC_peaks(i_Xp).right_bound,idx_iso));
    end
    % If the AUC of the peak of an IMP is less than 10% of the maximum,
    % filter, remove the proportion.
    esti_ratio(XIC_peaks(i_Xp).left_bound:XIC_peaks(i_Xp).right_bound,area_filter<max(area_filter)*obj.m_resFilterThres) = 0;
end
esti_ratio = esti_ratio./(sum(esti_ratio,2)+eps);

% calculate XIC of each IMPs using total XIC and ratio curve
intensityMatrix = esti_ratio.*smoothed_intensity;
auxic = zeros(num_iso,1);
rt_bound = repmat(struct('start',0,'end',0), num_iso, length(XIC_peaks));
idx_selected = zeros(num_iso,1);
ratio_each_XIC_peak = zeros(num_iso,length(XIC_peaks));
for idx_iso = 1:num_iso
    max_proportions = zeros(length(XIC_peaks),1);
    fwhm = zeros(length(XIC_peaks),1);
    for i_Xp = 1:length(XIC_peaks)
        max_proportions(i_Xp) = max(esti_ratio(XIC_peaks(i_Xp).left_bound:XIC_peaks(i_Xp).right_bound,idx_iso));
        % Calculate the fwhm for XIC selection
        peak_rts = rt_grid(XIC_peaks(i_Xp).left_bound:XIC_peaks(i_Xp).right_bound);
        peak_intens = intensityMatrix(XIC_peaks(i_Xp).left_bound:XIC_peaks(i_Xp).right_bound,idx_iso);
        fwhm(i_Xp) = CChromatogramUtils.get_fwhm(peak_rts, peak_intens);
        % Calculate the ratio of each IMP in each XIC peak
        ratio_each_XIC_peak(idx_iso, i_Xp) = ...
            trapz(rt_grid(XIC_peaks(i_Xp).left_bound:XIC_peaks(i_Xp).right_bound),...
            intensityMatrix(XIC_peaks(i_Xp).left_bound:XIC_peaks(i_Xp).right_bound,idx_iso))*60;
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
    final_rt_start = XIC_peaks(idx_selected(idx_iso)).left_bound;
    if final_rt_start ~= 1
        final_rt_start = final_rt_start - 1;
    end
    final_rt_end = XIC_peaks(idx_selected(idx_iso)).right_bound;
    if final_rt_end ~= length(rt_grid)
        final_rt_end = final_rt_end + 1;
    end
    auxic(idx_iso,1) = trapz(rt_grid(final_rt_start:final_rt_end),...
        [0;intensityMatrix(final_rt_start+1:final_rt_end-1,idx_iso);0])*60;
    for idx_p = 1:length(XIC_peaks)
        rt_bound(idx_iso,idx_p).start = rt_grid(XIC_peaks(idx_p).left_bound);
        rt_bound(idx_iso,idx_p).end = rt_grid(XIC_peaks(idx_p).right_bound);
    end
end

% Get the non-zero area under XIC, index and rt_bound
idxNonZero = find(auxic(:,1)~=0);
auxic = auxic(idxNonZero,:);
rt_bound = rt_bound(idxNonZero,:);
idx_selected = idx_selected(idxNonZero,:);
ratio_each_XIC_peak = ratio_each_XIC_peak(idxNonZero,:);
ratio_each_XIC_peak = ratio_each_XIC_peak./sum(ratio_each_XIC_peak,1);
if ~isempty(idxNonZero)
    bhave_non_zeros = true;
end
end