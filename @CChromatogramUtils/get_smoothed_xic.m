function [xic_rt, xic_intensity_smoothed, xic_intensity_raw] = get_smoothed_xic(ms12DatasetIO, raw_name, low_mz_bound, high_mz_bound, selected_charge)
    % Get Smoothed XIC from MS1 data
    % Inputs:
    %   ms12DatasetIO (object)
    %       Instance of CMS12DatasetIO
    %   raw_name (1 x 1 char/string)
    %       The raw file name
    %   low_mz_bound (1 x 1 double) m/z
    %       Lower m/z bound
    %   high_mz_bound (1 x 1 double) m/z
    %       Upper m/z bound
    %   selected_charge (1 x 1 double/int)
    %       Charge state
    % Output:
    %   xic_rt (N x 1 double) minutes
    %       Retention time grid
    %   xic_intensity_smoothed (N x 1 double) intensity
    %       Smoothed XIC intensity
    %   xic_intensity_raw (N x 1 double) intensity
    %       Raw XIC intensity (monoisotopic)
    
    % MS1_index (scan, retention time, peak number, baseline, injection time)
    % MS1_peaks (m/z, intensity)
    mgf_stem = erase(raw_name,'.mgf');
    
    % Accessing properties from ms12DatasetIO
    % Assuming ms12DatasetIO has m_cMsFileMapper, m_mapNameMS1Index, m_mapNameMS1Peaks
    ms1_stem = ms12DatasetIO.m_cMsFileMapper.get_ms1_stem(mgf_stem);
    MS1_index = ms12DatasetIO.m_mapNameMS1Index(ms1_stem);
    MS1_peaks = ms12DatasetIO.m_mapNameMS1Peaks(ms1_stem);
    
    xic_rt = MS1_index(:,2);
    isotope_num = [-1,0,1,2,3,4];
    xic_intensity_raw = zeros(size(MS1_index,1),length(isotope_num)); % retention time -> intensity, XIC
    
    % record the isotopic XIC
    % Using CConstant.unitdiff
    for idx_iso = 1:length(isotope_num)
        idxs_target_peaks = find(MS1_peaks(:,1) > low_mz_bound + isotope_num(idx_iso) * CConstant.unitdiff / selected_charge...
            & MS1_peaks(:,1) < high_mz_bound + isotope_num(idx_iso) * CConstant.unitdiff / selected_charge);
        for idx_itp = 1:length(idxs_target_peaks)
            xic_intensity_raw(find(MS1_index(:,3) > idxs_target_peaks(idx_itp), 1), idx_iso) = ...
                MS1_peaks(idxs_target_peaks(idx_itp), 2);
        end
    end
    
    % filter with two criteria:
    % 1. the intensity of -1 peak should not greater than monoisotopic
    % 2. the intensity of isotopic cluster peak should be enough similar with
    %   the IPV matrix.
    for idx_inten = 1:size(xic_intensity_raw,1)
        if xic_intensity_raw(idx_inten,1) > xic_intensity_raw(idx_inten,2) % the first criterion
            xic_intensity_raw(idx_inten,:) = 0;
        elseif ~any(xic_intensity_raw(idx_inten,:))
            continue;
        else
            % CConstant.IPV check
            obs_inten = xic_intensity_raw(idx_inten,2:end);
            obs_inten = obs_inten / max(obs_inten);
            cluster_check = [CConstant.IPV(int64((high_mz_bound+low_mz_bound)/2),:); obs_inten];
            if 1 - pdist(cluster_check, 'cosine') < 0.6
            xic_intensity_raw(idx_inten,:) = 0;
            end
        end
    end
    
    xic_intensity_raw = xic_intensity_raw(:,2);
    xic_intensity_smoothed = smoothdata(xic_intensity_raw, 'movmean', 5);
    % xic_intensity_smoothed = smoothdata(xic_intensity_raw,'gaussian',1,'SamplePoints',xic_rt);
end
