classdef CIMPGatherQuant
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
        m_mapRawNames;      % map of raw names to index in m_rawData

        % The following property stores per-raw grouped data
        % Each cell element is a struct with fields:
        %   length, capacity, curRts, curIntens, curMz, curCharge,
        %   ratioMatrix, impNames, mapIMPNames, impMass
        m_rawData;
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
            obj.m_mapRawNames = containers.Map();
            % different in different raw
            obj.m_rawData = {};
        end
        
        % Main entry point for summarizing the quantification of various modified form of a peptide
        runGather(obj);

        % Add one record
        obj = appendOneSpecQuant(obj,raw_name,curRts,curIntens,curMz,cur_ch,cstrIMP,lfMasses,abundance);
        

        % Quantify each group
        [has_nonzero_imp, imp_idx_nonzero, area_imp_final, rt_bound, idx_selected, ratio_each_XIC_peak] = ...
            quant_each_group(obj, raw_name, ratio_raw, rt_raw, ...
            intensity_raw, low_mz_bound, high_mz_bound, selected_charge);

        % Get the m/z bound of ms1 peak
        [low_bound,high_bound, selected_charge, charge_group_idxs] = ...
            get_mz_bound(obj, current_imp_mass, current_charge);

        % Re-quantification for gathered peptides using manually-checked rt range
        rerunGather_quant(obj,pep_rtrange_map);

        % Re-quantify each group
        [has_nonzero_imp, imp_idx_nonzero, area_imp_final, rt_bound, max_label, ratio_each_XIC_peak] = ...
            requant_each_group(obj, raw_name, ratio_raw, rt_raw, ...
            intensity_raw, low_mz_bound, high_mz_bound, selected_charge, ...
            current_imp_rt_range);

        % Draw the XIC for gathered peptides using manually-checked rt range to dir_save
        % The color_map and legend_map are optional
        drawGather(obj, pep_rtrange_map, dir_save, color_map, legend_map);

        % Draw the XIC for each group
        % The color_map and legend_map are optional
        draw_each_group(obj, raw_name, ratio_raw, rt_raw, ...
            intensity_raw, low_mz_bound, high_mz_bound, selected_charge, ...
            current_imp_rt_range, current_imp_name, current_imp_mass, ...
            current_charge, dir_save, color_map, legend_map)

        % Check whether the XIC peaks have at least min_rows PSMs
        has_min_rows = hasMinRows(obj, ratio_matrix, min_rows);
    end

    methods (Access=private)
        function [obj, idx_raw] = ensure_raw(obj, raw_name)
            if ~obj.m_mapRawNames.isKey(raw_name)
                idx_raw = obj.m_mapRawNames.Count + 1;
                obj.m_mapRawNames(raw_name) = idx_raw;
                obj.m_rawData{idx_raw} = obj.init_raw();
            else
                idx_raw = obj.m_mapRawNames(raw_name);
            end
        end

        function raw = get_raw(obj, idx_raw)
            raw = obj.m_rawData{idx_raw};
        end

        function obj = set_raw(obj, idx_raw, raw)
            obj.m_rawData{idx_raw} = raw;
        end

        function raw = init_raw(obj)
            raw = struct();
            raw.length = 0;
            raw.capacity = obj.m_buff_length;
            raw.curRts = zeros(obj.m_buff_length,1);
            raw.curIntens = zeros(obj.m_buff_length,1);
            raw.curMz = zeros(obj.m_buff_length,1);
            raw.curCharge = zeros(obj.m_buff_length,1);
            raw.ratioMatrix = zeros(obj.m_buff_length,0);
            raw.impNames = cell(0,1);
            raw.impMass = [];
            raw.mapIMPNames = containers.Map();
        end
    end
end

