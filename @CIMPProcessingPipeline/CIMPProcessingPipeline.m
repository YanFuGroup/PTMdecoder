classdef CIMPProcessingPipeline < handle
    % Processing pipeline for IMP quant/requant/draw

    properties (Access=private)
        m_ms12DatasetIO
        m_ms1_tolerance
        m_minMSMSnum
        m_alpha
        m_resFilterThres
    end

    methods
        function obj = CIMPProcessingPipeline(ms12DatasetIO, ms1_tolerance, minMSMSnum, alpha, resFilterThres)
            % Construct a pipeline with shared processing parameters
            obj.m_ms12DatasetIO = ms12DatasetIO;
            obj.m_ms1_tolerance = ms1_tolerance;
            obj.m_minMSMSnum = minMSMSnum;
            obj.m_alpha = alpha;
            obj.m_resFilterThres = resFilterThres;
        end
    end
end
