classdef CIMPGatherQuant < handle
    % A class for summarizing quantification at the IMP level
    
    properties(Access=public)
        m_buff_length;      % the length of the buffer
        m_prot_names_pos;   % the names of the proteins containing this peptide
        m_cMs12DatasetIO;   % MS1 and MS2 spectrum information index
        m_resFilterThres;   % Threshold for filtering results using relative intensity
        m_ms1_tolerance;    % the tolerance of MS1
        m_alpha;            % the filter threshold factor, thres is max*alpha
        m_outputPath;       % output path of the result file
        m_minMSMSnum;      % Minimum number of MSMS spectra for a peptide to be considered
        m_groupAggregator; % group aggregation helper
    end
    
    properties(Access=private)
        m_mapRawNames;      % map of raw names to index in m_rawIdentData

        % The following property stores per-raw grouped data
        % Each cell element is a CIMPRawIdentStore
        m_rawIdentData;
    end
    
    methods
        function obj = CIMPGatherQuant(prot_names_pos,cMs12DatasetIO,...
            resFilterThres,ms1_tolerance, alpha, outputPath, minMSMSnum)
            % Input:
            %   prot_names_pos (P x 2 cell)
            %       protein name and start position pairs
            %   cMs12DatasetIO (object)
            %       MS1/MS2 dataset IO instance
            %   resFilterThres (1 x 1 double)
            %       threshold for filtering results using relative intensity
            %   ms1_tolerance (struct)
            %       MS1 tolerance (fields: isppm, value)
            %   alpha (1 x 1 double)
            %       peak boundary factor, thres is max*alpha
            %   outputPath (1 x 1 char/string)
            %       output path of the result file
            %   minMSMSnum (1 x 1 double/int)
            %       minimum number of MSMS spectra for a peptide to be considered
            obj.m_buff_length = 50;
            obj.m_prot_names_pos = prot_names_pos;
            obj.m_cMs12DatasetIO = cMs12DatasetIO;
            obj.m_resFilterThres = resFilterThres;
            obj.m_ms1_tolerance = ms1_tolerance;
            obj.m_alpha = alpha;
            obj.m_outputPath = outputPath;
            if nargin < 7
                minMSMSnum = 1; % Default minimum number of MSMS spectra
            end
            obj.m_minMSMSnum = minMSMSnum;
            obj.m_groupAggregator = CIMPGroupAggregator(obj.m_ms1_tolerance);
            obj.m_mapRawNames = containers.Map('KeyType','char','ValueType','any');
            % different in different raw
            obj.m_rawIdentData = {};
        end
        
        % Main entry point for summarizing the quantification of various modified form of a peptide
        block = quantifyIMPs(obj);

        % Re-quantification for gathered peptides using manually-checked rt range
        block = requantifyIMPsWithRTRanges(obj,pep_rtrange_map);

        % Draw the XIC for gathered peptides using manually-checked rt range to dir_save
        % The color_map and legend_map are optional
        drawGather(obj, pep_rtrange_map, dir_save, color_map, legend_map);

        function rawStore = getRawStore(obj, raw_name)
            idx_raw = obj.ensure_raw(raw_name);
            rawStore = obj.get_raw(idx_raw);
        end

        function setRawStore(obj, raw_name, rawStore)
            idx_raw = obj.ensure_raw(raw_name);
            obj.set_raw(idx_raw, rawStore);
        end
    end

    methods (Access=private)
        function idx_raw = ensure_raw(obj, raw_name)
            if ~obj.m_mapRawNames.isKey(raw_name)
                idx_raw = obj.m_mapRawNames.Count + 1;
                obj.m_mapRawNames(raw_name) = idx_raw;
                obj.m_rawIdentData{idx_raw} = obj.init_raw();
            else
                idx_raw = obj.m_mapRawNames(raw_name);
            end
        end

        function raw = get_raw(obj, idx_raw)
            raw = obj.m_rawIdentData{idx_raw};
        end

        function set_raw(obj, idx_raw, raw)
            obj.m_rawIdentData{idx_raw} = raw;
        end

        function raw = init_raw(obj)
            raw = CIMPRawIdentStore(obj.m_buff_length);
        end

        function [raw_names, raw_ident_stores] = getRawEntries(obj)
            raw_names = obj.m_mapRawNames.keys;
            raw_ident_stores = cell(size(raw_names));
            for idx_keys = 1:obj.m_mapRawNames.Count
                idx_r = obj.m_mapRawNames(raw_names{idx_keys});
                raw_ident_stores{idx_keys} = obj.get_raw(idx_r);
            end
        end
    end
end

