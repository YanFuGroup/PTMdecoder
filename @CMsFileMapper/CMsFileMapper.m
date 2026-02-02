classdef CMsFileMapper < handle
    % CMsFileMapper - Maps MGF files to their corresponding MS1 files using suffix rules.
    % This class scans the directory once and provides a lookup service.
    
    properties
        m_specPath;         % Folder path
        m_mgf2ms1_map;      % Map: MGF stem -> MS1 stem
        m_suffixes = {'_HCDFT', '_HCDIT', '_CIDFT', '_CIDIT'};
    end
    
    methods
        function obj = CMsFileMapper(specPath)
            % Input:
            %   specPath (1 x 1 char/string)
            %       spectra folder path
            obj.m_specPath = specPath;
            obj.m_mgf2ms1_map = containers.Map();
            obj.build_mapping();
        end
        
        % Build mgf->ms1 mapping based on file stems and suffix rules
        build_mapping(obj)
        % Returns the stem of the MS1 file corresponding to the given MGF stem
        ms1_stem = get_ms1_stem(obj, mgf_stem)
        % Returns the stem of the MS2 file corresponding to the given MGF stem
        ms2_stem = get_ms2_stem(obj, mgf_stem)
    end
end