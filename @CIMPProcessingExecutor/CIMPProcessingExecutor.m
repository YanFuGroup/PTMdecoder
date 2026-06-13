classdef CIMPProcessingExecutor < handle
    % Execution layer for IMP quant/requant/draw

    properties (Access=private)
        m_ms12DatasetIO
        m_ms1_tolerance
        m_minMSMSnum
        m_alpha
        m_resFilterThres
        m_minXicNonzeroPoints
    end

    methods
        function obj = CIMPProcessingExecutor(ms12DatasetIO, ms1_tolerance, minMSMSnum, alpha, resFilterThres)
            % Construct an executor with shared processing parameters
            if nargin == 1 && isa(ms12DatasetIO, 'CIMPProcessingExecutorConfig')
                cfg = ms12DatasetIO;
                obj.m_ms12DatasetIO = cfg.ms12DatasetIO;
                obj.m_ms1_tolerance = cfg.ms1_tolerance;
                obj.m_minMSMSnum = cfg.minMSMSnum;
                obj.m_alpha = cfg.alpha;
                obj.m_resFilterThres = cfg.resFilterThres;
                obj.m_minXicNonzeroPoints = cfg.minXicNonzeroPoints;
                return;
            end

            obj.m_ms12DatasetIO = ms12DatasetIO;
            obj.m_ms1_tolerance = ms1_tolerance;
            obj.m_minMSMSnum = minMSMSnum;
            obj.m_alpha = alpha;
            obj.m_resFilterThres = resFilterThres;
            obj.m_minXicNonzeroPoints = 5;
        end
    end

    methods (Access=public)
        base_groups = buildBaseGroups(obj, rawIdentManager)

        % Main quantification/re-quantification/drawing methods
        block = quantifyPeptideBlock(obj, prot_names_pos, rawIdentManager, base_groups)
        block = requantifyPeptideBlock(obj, prot_names_pos, rawIdentManager, pep_rtrange_map, base_groups)
        drawImpXicForBlock(obj, rawIdentManager, pep_rtrange_map, dir_save, color_map, legend_map)
    end

    methods (Access=private)
        % Group processing callbacks for quantification, re-quantification, and drawing
        imp_records = onGroupQuant(obj, imp_records, group)
        imp_records = onGroupRequant(obj, imp_records, group)
        state = onGroupDrawXic(obj, state, group)
    end
end
