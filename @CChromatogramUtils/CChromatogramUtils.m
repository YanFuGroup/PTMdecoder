classdef CChromatogramUtils
    % CChromatogramUtils
    % A utility class for chromatogram signal processing, including
    % smoothing, peak detection, and data preprocessing.
    
    methods (Static)
        function [sort_rts, sort_inten, sort_ratioMatrix, is_valid] = preprocess_ms1_inputs(current_rts, current_inten, current_ratioMatrix, minMSMSnum)
            % Preprocess MS1 inputs: Sort by retention time, Smooth intensity, and Denoise
            % Inputs:
            %   current_rts: Array of retention times
            %   current_inten: Array of intensities
            %   current_ratioMatrix: quantification matrix
            %   minMSMSnum: Minimum number of MSMS recurrence required
            
            is_valid = true;
            
            % Sort MS1 signal (pair of retention time and intensity) by time
            sort_rts = [(1:length(current_rts))',current_rts];
            sort_rts = sortrows(sort_rts,2); % Sort in ascending order by retention time column
            sort_idx = sort_rts(:,1);
            sort_rts = sort_rts(:,2);
            sort_inten = current_inten(sort_idx);
            sort_ratioMatrix = current_ratioMatrix(sort_idx,:); % Rearrange the matrix in chronological order
            
            % Sort and denoise using a relative abundance threshold method
            % Using loess smoothing with span 0.05
            if length(sort_inten) > 4 % loess requires enough points
                sort_inten = smooth(sort_inten, 0.05, 'loess');
            end
            
            maxInten = max(sort_inten);
            if isempty(maxInten)
                maxInten = 0;
            end
            
            tmp = sort_inten < 0.05 * maxInten; % Find results where intensity is less than 0.05 of the maximum abundance and discard them
            sort_inten(tmp) = []; 
            sort_rts(tmp) = [];
            sort_ratioMatrix(tmp,:) = [];

            % Check if there are enough rows left
            % Replaces obj.hasMinRows(sort_ratioMatrix, minMSMSnum)
            if size(sort_ratioMatrix, 1) < minMSMSnum
                % If the ratio matrix has less than min rows, skip this group
                is_valid = false;
                sort_inten = [];
                sort_ratioMatrix = [];
            end
            % TODO: The filtering is not only here, but also the XIC integration area is determined later.
            %  Then, within the integration area, it is checked how many fall within this area.
            %  If there are not enough, they will be removed.
        end

        function [rt_grid, smoothed_intensity, intensity] = get_smoothed_xic(ms12DatasetIO, raw_name, low_mz_bound, high_mz_bound, selected_charge)
            % Get Smoothed XIC from MS1 data
            % Inputs:
            %   ms12DatasetIO: Instance of CMS12DatasetIO
            %   raw_name: The raw file name
            %   low_mz_bound: Lower m/z bound
            %   high_mz_bound: Upper m/z bound
            %   selected_charge: Charge state
            
            % MS1_index (scan, retention time, peak number, baseline, injection time)
            % MS1_peaks (m/z, intensity)
            mgf_stem = erase(raw_name,'.mgf');
            
            % Accessing properties from ms12DatasetIO
            % Assuming ms12DatasetIO has m_cMsFileMapper, m_mapNameMS1Index, m_mapNameMS1Peaks
            ms1_stem = ms12DatasetIO.m_cMsFileMapper.get_ms1_stem(mgf_stem);
            MS1_index = ms12DatasetIO.m_mapNameMS1Index(ms1_stem);
            MS1_peaks = ms12DatasetIO.m_mapNameMS1Peaks(ms1_stem);
            
            rt_grid = MS1_index(:,2);
            isotope_num = [-1,0,1,2,3,4];
            intensity = zeros(size(MS1_index,1),length(isotope_num)); % retention time -> intensity, XIC
            
            % record the isotopic XIC
            % Using CConstant.unitdiff
            for idx_iso = 1:length(isotope_num)
                idxs_target_peaks = find(MS1_peaks(:,1) > low_mz_bound + isotope_num(idx_iso) * CConstant.unitdiff / selected_charge...
                    & MS1_peaks(:,1) < high_mz_bound + isotope_num(idx_iso) * CConstant.unitdiff / selected_charge);
                for idx_itp = 1:length(idxs_target_peaks)
                    intensity(find(MS1_index(:,3) > idxs_target_peaks(idx_itp), 1), idx_iso) = ...
                        MS1_peaks(idxs_target_peaks(idx_itp), 2);
                end
            end
            
            % filter with two criteria:
            % 1. the intensity of -1 peak should not greater than monoisotopic
            % 2. the intensity of isotopic cluster peak should be enough similar with
            %   the IPV matrix.
            for idx_inten = 1:size(intensity,1)
                if intensity(idx_inten,1) > intensity(idx_inten,2) % the first criterion
                    intensity(idx_inten,:) = 0;
                elseif ~any(intensity(idx_inten,:))
                    continue;
                else
                    % CConstant.IPV check
                    obs_inten = intensity(idx_inten,2:end);
                    obs_inten = obs_inten / max(obs_inten);
                    cluster_check = [CConstant.IPV(int64((high_mz_bound+low_mz_bound)/2),:); obs_inten];
                    if 1 - pdist(cluster_check, 'cosine') < 0.6
                    intensity(idx_inten,:) = 0;
                    end
                end
            end
            
            intensity = intensity(:,2);
            smoothed_intensity = smoothdata(intensity, 'movmean', 5);
            % smoothed_intensity = smoothdata(intensity,'gaussian',1,'SamplePoints',rt_grid);
        end

        function fwhm = get_fwhm(peak_rts, peak_intens)
            % Get the full width at half maxima (FWHM)
            % Input:
            %   peak_rts
            %       retention times of peaks
            %   peak_intens
            %       intensities of peaks
            % Output:
            %   fwhm
            %       full width at half maxima
            
            % Default, used for zero-intensity peak or null peak
            fwhm = 0;
            
            % Check whether the lengths of peak_rts and peak_intens are equal
            if length(peak_rts) ~= length(peak_intens)
                error('The length of peak_rts (%d) and peak_intens (%d) is not equal',...
                    length(peak_rts), length(peak_intens));
            end
            
            % find two ranges of half maxima
            half_maxima = max(peak_intens)/2;
            if half_maxima==0
                return;
            end
            % the two point around the left half maxima, find left half bound
            left_right = find(peak_intens>half_maxima, 1);
            left_left = left_right - 1;
            if left_left > 0
                left_half_point = peak_rts(left_left)+(half_maxima-peak_intens(left_left))*...
                    (peak_rts(left_right)-peak_rts(left_left))/(peak_intens(left_right)-peak_intens(left_left));
            else
                % when the left half maxima is smaller than the most left point in
                %   this peak, use the most left point to measure the peak width
                left_half_point = peak_rts(1);
            end
            % the two point around the right half maxima, find right half bound
            right_left = find(peak_intens>half_maxima, 1, 'last');
            right_right = right_left + 1;
            if right_right < length(peak_rts)
                right_half_point = peak_rts(right_right)-(half_maxima-peak_intens(right_right))*...
                    (peak_rts(right_right)-peak_rts(right_left))/(peak_intens(right_left)-peak_intens(right_right));
            else
                % when the left half maxima is smaller than the most left point in
                %   this peak, use the most left point to measure the peak width
                right_half_point = peak_rts(end);
            end
            fwhm = right_half_point - left_half_point;
        end
    end
end
