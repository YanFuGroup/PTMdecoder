% MS1 and MS2 spectrum data class
classdef CMS12DatasetIO<CDatasetIO
    properties(Access=public)
        m_mapNameMS1Index;  % Dictionary from file name to ms1_index
        m_mapNameMS1Peaks;  % Dictionary from file name to ms1_peaks
        m_mapNameMS2Index;  % Dictionary from file name to ms2_index
        m_ms1_tolerance;    % the mass tolerance of MS1
    end

    methods  
        function obj = CMS12DatasetIO(strDatasetFoldname,ms1_tolerance)
            obj.m_strFoldname=strDatasetFoldname;
            obj.m_mapNameMS1Index = containers.Map();
            obj.m_mapNameMS1Peaks = containers.Map();
            obj.m_mapNameMS2Index = containers.Map();
            obj.m_ms1_tolerance = ms1_tolerance;
        end
        % Generate spectrum index MS1_index and peak index MS1_peaks using .ms1 file
        success = load_MS1_file(obj,ms1_fullfile);
        % Generate spectrum index MS2_index and peak index MS2_peaks using .ms2 file
        success = load_MS2_file(obj,ms2_fullfile);
        % Generate mapping between MS1 and MS2 scans using .ms2 file
        success = load_MS1_MS2_mapping(obj, ms2_fullfile);
        % Build a dictionary mapping spectrum names to corresponding index or peaks
        SetMap(obj);
        % Output more accurate mass-to-charge ratio
        acc_mz = get_acc_mz(obj,cen_mz,cur_mz,cur_chg);
    end
end

