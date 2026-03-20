% MS1 and MS2 spectrum data class
classdef CMS12DatasetIO<CDatasetIO
    properties(Access=public)
        m_mapNameMS1Index;  % Dictionary from file name to ms1_index
        m_mapNameMS1Peaks;  % Dictionary from file name to ms1_peaks
        m_mapNameMS1SortedMz;    % Dictionary from file name to globally sorted MS1 m/z
        m_mapNameMS1SortedInt;   % Dictionary from file name to intensity aligned to sorted m/z
        m_mapNameMS1SortedScan;  % Dictionary from file name to scan row aligned to sorted m/z
        m_mapNameMS2Index;  % Dictionary from file name to ms2_index
        m_ms1_tolerance;    % the mass tolerance of MS1
        m_cMsFileMapper;    % Handle to CMsFileMapper for MGF<->MS1 mapping
    end

    methods  
        function obj = CMS12DatasetIO(strDatasetFoldname,ms1_tolerance)
            % Input:
            %   strDatasetFoldname (1 x 1 char/string)
            %       dataset folder path
            %   ms1_tolerance (struct)
            %       MS1 tolerance (fields: value, isppm)
            obj.m_strFoldname=strDatasetFoldname;
            obj.m_mapNameMS1Index = containers.Map();
            obj.m_mapNameMS1Peaks = containers.Map();
            obj.m_mapNameMS1SortedMz = containers.Map();
            obj.m_mapNameMS1SortedInt = containers.Map();
            obj.m_mapNameMS1SortedScan = containers.Map();
            obj.m_mapNameMS2Index = containers.Map();
            obj.m_ms1_tolerance = ms1_tolerance;
            obj.m_cMsFileMapper = CMsFileMapper(strDatasetFoldname);

            % One-step setup: build MS1/MS2 index maps
            obj.SetMap();
        end
        % Generate spectrum index MS1_index and peak index MS1_peaks using .ms1 file
        success = load_MS1_file(obj,ms1_fullfile);
        % Generate spectrum index MS2_index and peak index MS2_peaks using .ms2 file
        success = load_MS2_file(obj,ms2_fullfile);
        % Generate mapping between MS1 and MS2 scans using .ms2 file
        success = load_MS1_MS2_mapping(obj, ms2_fullfile);
        % Output more accurate mass-to-charge ratio
        acc_mz = get_acc_mz(obj,cen_mz,cur_mz,cur_chg);
    end

    methods (Access=protected)
        % Build a dictionary mapping spectrum names to corresponding index or peaks
        SetMap(obj);
    end
end

