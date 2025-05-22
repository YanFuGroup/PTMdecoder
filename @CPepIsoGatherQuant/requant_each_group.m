function [bhave_non_zeros, idxNonZero, auxic, rt_bound, max_label, ratio_each_XIC_peak]...
    = requant_each_group(obj,raw_name,current_ratioMatrix,current_rts,...
    current_inten, low_mz_bound, high_mz_bound, selected_charge,...
    current_iso_rt_range)
% Re-quantify each group
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
%   current_iso_rt_range
%       retention times of current IMPs
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
%   max_label
%       the max check label of the XIC peaks for each IMP

% System error in saving retention time
eps_rt_print = 1e-6;

bhave_non_zeros = false;
num_iso = size(current_ratioMatrix,2);
rt_error_tol = 1; % RT match tolerance, choose 1 arbitrarily
% A vector showing is needed to skip this IMP.
%   Cannot delete because other filter can also change this vector
is_skip_vec = cellfun(@isempty,current_iso_rt_range);

% Sort MS1 signal (pair of retention time and intensity) by time
sort_rts = [(1:length(current_rts))',current_rts]; % Add an ordered number in front
sort_rts = sortrows(sort_rts,2); % Sort in ascending order by the second column
sort_idx = sort_rts(:,1);
sort_rts = sort_rts(:,2);
sort_inten = current_inten(sort_idx);
sort_ratioMatrix = current_ratioMatrix(sort_idx,:); % Rearrange the matrix according to time order

% Smooth sort_inten, denoise using a relative abundance threshold method
sort_inten = smooth(sort_inten,0.05,'loess'); %0.05,'loess'
maxInten = max(sort_inten);
tmp = sort_inten<0.05*maxInten; % Find results where intensity is less than 0.05 of the maximum value and discard them
sort_inten(tmp) = []; %#ok<NASGU>
sort_rts(tmp) = [];
sort_ratioMatrix(tmp,:) = [];

% find the XIC filtered by m/z of base peak [low_mz_bound, high_mz_bound]
%   and -1,+0,+1,+2,+3
% MS1_index (scan, retention time, peak number, baseline, injection time)
% MS1_peaks (m/z, intensity)
MS1_index = obj.m_cMs12DatasetIO.m_mapNameMS1Index(erase(raw_name,'.mgf'));
MS1_peaks = obj.m_cMs12DatasetIO.m_mapNameMS1Peaks(erase(raw_name,'.mgf'));
rt_grid = MS1_index(:,2);
isotope_num = [-1,0,1,2,3,4];
intensity = zeros(size(MS1_index,1),length(isotope_num)); % retention time -> intensity, XIC
% record the isotopic XIC
for idx_iso = 1:length(isotope_num)
    idxs_target_peaks = find(MS1_peaks(:,1)>low_mz_bound+isotope_num(idx_iso)*CConstant.unitdiff/selected_charge...
        & MS1_peaks(:,1)<high_mz_bound+isotope_num(idx_iso)*CConstant.unitdiff/selected_charge);
    for idx_itp = 1:length(idxs_target_peaks)
        intensity(find(MS1_index(:,3)>idxs_target_peaks(idx_itp),1),idx_iso) = ...
            MS1_peaks(idxs_target_peaks(idx_itp),2);
    end
end
% filter with two criteria:
% 1. the intensity of -1 peak should not greater than monoisotopic
% 2. the intensity of isotopic cluster peak should be enough similar with
%   the IPV matrix.
for idx_inten = 1:size(intensity,1)
    if intensity(idx_inten,1)>intensity(idx_inten,2) % the first criterion
        intensity(idx_inten,:) = 0;
    elseif ~any(intensity(idx_inten,:))
        continue;
    elseif 1-pdist([CConstant.IPV(int64((high_mz_bound+low_mz_bound)/2),:);...
            intensity(idx_inten,2:end)],'cosine') < 0.6
        intensity(idx_inten,:) = 0;
    end
end
intensity = intensity(:,2);
smoothed_intensity = smoothdata(intensity,'movmean',5);
% smoothed_intensity = smoothdata(intensity,'gaussian',1,'SamplePoints',rt_grid);

% Extract the rt bound of XIC peak
final_XIC_peak_for_IMP = repmat(struct('left_bound',0,'right_bound',0), num_iso, 1);
max_label = zeros(num_iso,1);
for idx_iso = 1:num_iso
    % Check if need to consider this IMP
    if is_skip_vec(idx_iso)
        continue;
    end

    % Record all of the check labels
    check_labels = zeros(length(current_iso_rt_range{idx_iso}),1);
    for idx_peak = 1:length(current_iso_rt_range{idx_iso})
        check_labels(idx_peak) = current_iso_rt_range{idx_iso}(idx_peak).check_label;
    end
    % Find the peak with max check label (the first of max peaks)
    [max_label(idx_iso), idx_max] = max(check_labels);
    if max_label(idx_iso) == 0
        % If all check labels are zero, skip this IMP later on
        is_skip_vec(idx_iso)=true;
        continue;
    end
    final_XIC_peak_for_IMP(idx_iso).left_bound = current_iso_rt_range{idx_iso}(idx_max).rt_start;
    final_XIC_peak_for_IMP(idx_iso).right_bound = current_iso_rt_range{idx_iso}(idx_max).rt_end;
end

% Calculate the ratio on each XIC points using kernel method, and normalize
% using Nadaraya-Waston kernel averaging method
esti_ratio = zeros(length(rt_grid),num_iso);
% % bandwidth = (4/(3*size(sort_rts,1)))^0.2*std(sort_rts);
for idx_iso = 1:num_iso
    % Check if need to consider this IMP
    if is_skip_vec(idx_iso)
        continue;
    end

    % Gaussian kernel function
    Ker_Gaussian = @(u) (1/sqrt(2*pi))*exp(-0.5*u.^2);
    % Epanechnikov kernel function
    %     Ker_Epanechnikov = @(u) (3/4)*(1-u.^2).*(abs(u)<=1);
    
    % Retention time bound of current IMP
    cur_iso_rt_left = final_XIC_peak_for_IMP(idx_iso).left_bound;
    cur_iso_rt_right = final_XIC_peak_for_IMP(idx_iso).right_bound;
    [rt_diff, idx_rt_left] = min(abs(rt_grid-cur_iso_rt_left));
    if rt_diff > rt_error_tol
        error(['Cannot find the spectra on the specified retention time: ', ...
            num2str(cur_iso_rt_left)]);
    end
    [rt_diff, idx_rt_right] = min(abs(rt_grid-cur_iso_rt_right));
    if rt_diff > rt_error_tol
        error(['Cannot find the spectra on the specified retention time: ', ...
            num2str(cur_iso_rt_right)]);
    end

    % Collect all of the rts states within current XIC peak
    idxs_ident_rt = sort_rts>=cur_iso_rt_left-eps_rt_print & sort_rts<=cur_iso_rt_right+eps_rt_print;
    rts_current = sort_rts(idxs_ident_rt);

    % Calculate bandwidth and weights for each XIC peak
    bandwidth = (4/(3*size(rts_current,1)))^0.2*std(rts_current);
    weights = zeros(idx_rt_right-idx_rt_left+1, length(rts_current));
    for idx_PSM = 1:length(rts_current)
        if bandwidth == 0
            break;
        end
        % set kernel weights
        weights(:,idx_PSM) = Ker_Gaussian(...
            (rt_grid(idx_rt_left:idx_rt_right)-rts_current(idx_PSM))/bandwidth);
    end
    % Check if there are nearly no weights in some retention time for all
    %   IMP, or the bandwidth is just zero
    if bandwidth == 0 || any(all(weights<1e-15/length(sort_rts),2))
        bandwidth = min(cur_iso_rt_right-cur_iso_rt_left,1);
        for idx_PSM = 1:length(rts_current)
            weights(:,idx_PSM) = Ker_Gaussian(...
                (rt_grid(idx_rt_left:idx_rt_right)-rts_current(idx_PSM))/bandwidth);
        end
    end
    % Calculate the ratio using normalized weights and ratioMatrix
    esti_ratio(idx_rt_left:idx_rt_right,idx_iso) = ...
        (weights * sort_ratioMatrix(idxs_ident_rt,idx_iso))./(sum(weights,2)+eps);
end
% normalize the ratio in every available retention time
esti_ratio = esti_ratio./(sum(esti_ratio,2)+eps);

% Hypothesis: The retention time bounds revised manually cannot be wrong.
%   No need to recognize and delete the low abundance IMP in revised file.
% Filter low abundance IMP according to the relative AUXIC
% intensityMatrix = esti_ratio.*smoothed_intensity;
% for i_Xp = 1:length(XIC_peaks)
%     area_filter = zeros(num_iso,1);
%     for idx_iso = 1:num_iso
%         area_filter(idx_iso) = trapz(rt_grid(XIC_peaks(i_Xp).left_bound:XIC_peaks(i_Xp).right_bound),...
%             intensityMatrix(XIC_peaks(i_Xp).left_bound:XIC_peaks(i_Xp).right_bound,idx_iso));
%     end
%     % If the AUC of the peak of an IMP is less than 10% of the maximum,
%     % filter, remove the proportion.
%     esti_ratio(XIC_peaks(i_Xp).left_bound:XIC_peaks(i_Xp).right_bound,area_filter<max(area_filter)*obj.m_resFilterThres) = 0;
% end
% esti_ratio = esti_ratio./(sum(esti_ratio,2)+eps);

% % calculate XIC of each IMPs using total XIC and ratio curve
% intensityMatrix = esti_ratio.*smoothed_intensity;
% auxic = zeros(num_iso,1);
% rt_bound = repmat(struct('start',0,'end',0), num_iso, length(XIC_peaks));
% idx_selected = zeros(num_iso,1);
% ratio_each_XIC_peak = zeros(num_iso,length(XIC_peaks));
% for idx_iso = 1:num_iso
%     max_proportions = zeros(length(XIC_peaks),1);
%     %     each_area = zeros(length(XIC_peaks),1);
%     fwhm = zeros(length(XIC_peaks),1);
%     for i_Xp = 1:length(XIC_peaks)
%         max_proportions(i_Xp) = max(esti_ratio(XIC_peaks(i_Xp).left_bound:XIC_peaks(i_Xp).right_bound,idx_iso));
%         % calculate the fwhm for XIC selection
%         peak_rt = rt_grid(XIC_peaks(i_Xp).left_bound:XIC_peaks(i_Xp).right_bound);
%         peak_intens = intensityMatrix(XIC_peaks(i_Xp).left_bound:XIC_peaks(i_Xp).right_bound,idx_iso);
%         % the area of each peak
%         % find two ranges of half maxima
%         half_maxima = max(peak_intens)/2;
%         if half_maxima==0
%             continue;
%         end
%         % the two point around the left half maxima, find left half bound
%         left_right = find(peak_intens>half_maxima, 1);
%         left_left = left_right - 1;
%         if left_left > 0
%             left_half_point = peak_rt(left_left)+(half_maxima-peak_intens(left_left))*...
%                 (peak_rt(left_right)-peak_rt(left_left))/(peak_intens(left_right)-peak_intens(left_left));
%         else
%             % when the left half maxima is smaller than the most left point in
%             %   this peak, use the most left point to measure the peak width
%             left_half_point = 1;
%         end
%         % the two point around the right half maxima, find right half bound
%         right_left = find(peak_intens>half_maxima, 1, 'last');
%         right_right = right_left + 1;
%         if right_right < length(peak_rt)
%             right_half_point = peak_rt(right_right)-(half_maxima-peak_intens(right_right))*...
%                 (peak_rt(right_right)-peak_rt(right_left))/(peak_intens(right_left)-peak_intens(right_right));
%         else
%             % when the left half maxima is smaller than the most left point in
%             %   this peak, use the most left point to measure the peak width
%             right_half_point = length(peak_rt);
%         end
%         fwhm(i_Xp) = right_half_point - left_half_point;
%         %         each_area(i_Xp) = trapz(rt_grid(XIC_peaks(i_Xp).left_bound:XIC_peaks(i_Xp).right_bound),...
%         %             intensityMatrix(XIC_peaks(i_Xp).left_bound:XIC_peaks(i_Xp).right_bound,idx_iso))*60;
%         ratio_each_XIC_peak(idx_iso, i_Xp) = ...
%             trapz(rt_grid(XIC_peaks(i_Xp).left_bound:XIC_peaks(i_Xp).right_bound),...
%             intensityMatrix(XIC_peaks(i_Xp).left_bound:XIC_peaks(i_Xp).right_bound,idx_iso))*60;
%     end
%     % When there are more than one XIC peak counting for an IMP,
%     % only record one XIC peak, the choosing critera are:
%     % 1. the peak with max proportion, ratio shows the XIC belong more to this IMP;
%     % 2. the peak with largest fwhm (full width of half maximum).
%     [val_max, idx_selected(idx_iso)] = max(max_proportions);
%     if length(find(abs(val_max-max_proportions)<eps)) > 1
%         % If the proportion can not select the only one, choose according
%         %   to the fwhm
%         [~, idx_selected(idx_iso)] = max(fwhm);
%         %         % If the proportion can not select the only one, choose according
%         %         %   to the AUC
%         %         [~, idx_max] = max(each_area);
%     end
%     final_rt_start = XIC_peaks(idx_selected(idx_iso)).left_bound;
%     if final_rt_start ~= 1
%         final_rt_start = final_rt_start - 1;
%     end
%     final_rt_end = XIC_peaks(idx_selected(idx_iso)).right_bound;
%     if final_rt_end ~= length(rt_grid)
%         final_rt_end = final_rt_end + 1;
%     end
%     auxic(idx_iso,1) = trapz(rt_grid(final_rt_start:final_rt_end),...
%         [0;intensityMatrix(final_rt_start+1:final_rt_end-1,idx_iso);0])*60;
%     for idx_p = 1:length(XIC_peaks)
%         rt_bound(idx_iso,idx_p).start = rt_grid(XIC_peaks(idx_p).left_bound);
%         rt_bound(idx_iso,idx_p).end = rt_grid(XIC_peaks(idx_p).right_bound);
%     end
% end

% Requantification using revised RT
intensityMatrix = esti_ratio.*smoothed_intensity;
auxic = zeros(num_iso,1);
rt_bound = repmat(struct('start',0,'end',0), num_iso, 1);
ratio_each_XIC_peak = zeros(num_iso,1);
for idx_iso = 1:num_iso
    % Check if need to consider this IMP
    if is_skip_vec(idx_iso)
        continue;
    end

    % Get the final rt bound
    rt_bound(idx_iso).start = final_XIC_peak_for_IMP(idx_iso).left_bound;
    rt_bound(idx_iso).end = final_XIC_peak_for_IMP(idx_iso).right_bound;
    [rt_diff, final_rt_start] = min(abs(rt_grid-rt_bound(idx_iso).start));
    if rt_diff > rt_error_tol
        error(['Cannot find the spectra on the specified retention time: ', ...
            num2str(rt_bound(idx_iso).start)]);
    end
    if final_rt_start ~= 1
        final_rt_start = final_rt_start - 1;
    end
    [rt_diff, final_rt_end] = min(abs(rt_grid-rt_bound(idx_iso).end));
    if rt_diff > rt_error_tol
        error(['Cannot find the spectra on the specified retention time: ', ...
            num2str(rt_bound(idx_iso).end)]);
    end
    if final_rt_end ~= length(rt_grid)
        final_rt_end = final_rt_end + 1;
    end
    auxic(idx_iso,1) = trapz(rt_grid(final_rt_start:final_rt_end),...
        [0;intensityMatrix(final_rt_start+1:final_rt_end-1,idx_iso);0])*60;
    total_temp = trapz(rt_grid(final_rt_start:final_rt_end),...
        [0;smoothed_intensity(final_rt_start+1:final_rt_end-1);0])*60;
    ratio_each_XIC_peak(idx_iso,1) = auxic(idx_iso,1) / total_temp;
end

% Get the non-zero area under XIC, index and rt_bound
idxNonZero = find(auxic(:,1)~=0);
auxic = auxic(idxNonZero,:);
rt_bound = rt_bound(idxNonZero,:);
max_label = max_label(idxNonZero,:);
ratio_each_XIC_peak = ratio_each_XIC_peak(idxNonZero,:);
if ~isempty(idxNonZero)
    bhave_non_zeros = true;
end
end

