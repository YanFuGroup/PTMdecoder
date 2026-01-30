classdef CPepIsoGatherQuant
    % A class for summarizing quantification at the peptide level
    
    properties(Access=public)
        m_buff_length;      % the length of the buffer
        m_prot_names_pos;   % the names of the proteins containing this peptide
        m_cMs12DatasetIO;   % MS1 and MS2 spectrum information index
        m_resFilterThres;   % Threshold for filtering results using relative intensity
        m_ms1_tolerance;    % the tolerance of MS1
        m_alpha;            % the filter threshold factor, thres is max*alpha
        m_outputPath;       % output path of the result file
        m_minMSMSnum;      % Minimum number of MSMS spectra for a peptide to be considered
        m_resStr;           % result string
        m_mapRawNames;      % the map of raw names to which cell in m_curRts et.al.

        % The following properties are in cell format, each element means
        %   different clusters of different MS files
        m_length;           % the available value of each cluster
        m_capacity;         % the capacity of the each cluster
        m_curRts;           % Current retention time
        m_curIntens;        % Current intensity
        m_curMz;            % observed precursor m/zs of the MS2
        m_curCharge;        % precursor charge of the MS2
        m_ratioMatrix;      % Quantification matrix (can be relative or absolute abundance)
        m_cstrIMPNames; % IMP strings
        m_mapIMPNames;  % Map from IMP to the corresponding column in the quantification matrix
        m_IMPMass;      % the mass of each IMP
    end
    
    methods
        function obj = CPepIsoGatherQuant(prot_names_pos,cMs12DatasetIO,...
                resFilterThres,ms1_tolerance, alpha, outputPath, minMSMSnum)
            % common
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
            obj.m_resStr = [];
            obj.m_mapRawNames = containers.Map();
            % different in different raw
            obj.m_length = {};
            obj.m_capacity = {};
            obj.m_curRts = {};
            obj.m_curIntens = {};
            obj.m_curMz = {};
            obj.m_curCharge = {};
            obj.m_ratioMatrix = {};
            obj.m_cstrIMPNames = {};
            obj.m_mapIMPNames = {};
            obj.m_IMPMass = {};
        end
        
        % Main entry point for summarizing the quantification of various modified form of a peptide
        runGather(obj);

        % Add one record
        obj = appendOneSpecQuant(obj,raw_name,curRts,curIntens,curMz,cur_ch,cstrPepIso,lfMasses,abundance);
        

        % Clear useless rows
        obj = delUselessRaws(obj);

        % Quantify each group
        [bhave_non_zeros, idxNonZero, auxic, rt_bound, idx_selected, ratio_each_XIC_peak] = ...
            quant_each_group(obj, raw_name, current_ratioMatrix, current_rts, ...
            current_inten, low_mz_bound, high_mz_bound, selected_charge);

        % Get the m/z bound of ms1 peak
        [low_bound,high_bound, selected_charge, charge_group_idxs] = ...
            get_mz_bound(obj, current_iso_mass, current_charge);

        % Re-quantification for gathered peptides using manually-checked rt range
        rerunGather_quant(obj,pep_rtrange_map);

        % Re-quantify each group
        [bhave_non_zeros, idxNonZero, auxic, rt_bound, max_label, ratio_each_XIC_peak] = ...
            requant_each_group(obj, raw_name, current_ratioMatrix, current_rts, ...
            current_inten, low_mz_bound, high_mz_bound, selected_charge, ...
            current_iso_rt_range);

        % Draw the XIC for gathered peptides using manually-checked rt range to dir_save
        % The color_map and legend_map are optional
        drawGather(obj, pep_rtrange_map, dir_save, color_map, legend_map);

        % Draw the XIC for each group
        % The color_map and legend_map are optional
        draw_each_group(obj, raw_name, current_ratioMatrix, current_rts, ...
            current_inten, low_mz_bound, high_mz_bound, selected_charge, ...
            current_iso_rt_range, current_iso_name, current_iso_mass, ...
            current_charge, dir_save, color_map, legend_map)

        % Check whether the XIC peaks have at least min_rows PSMs
        is_reserved = hasMinRows(obj, ratio_matrix, min_rows);
    end
end

