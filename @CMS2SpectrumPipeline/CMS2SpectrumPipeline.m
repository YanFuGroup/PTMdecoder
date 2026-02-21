classdef CMS2SpectrumPipeline
    % Pipeline for MS2 spectrum-level identification and quantification

    properties
        m_pepSeq;
        m_isProtN;
        m_isProtC;
        m_cMgfDatasetIO;
        m_strDatasetName;
        m_strSpecName;
        m_fixedModNameMass;
        m_variableModNameMass;
        m_model;
        m_method;
        m_lambda;
        m_ms1_tolerance;
        m_ms2_tolerance;
        m_alpha;
        m_resFilterThres;
        m_ionTypes;
        m_enzyme;

        m_case_penalty_intens;
        m_grid_penalty_intens;
        m_case_OLS_intens_weight;

        m_iCharge;
        m_dPrecursorMass;
        m_expPeaks;
    end

    methods
        function obj = CMS2SpectrumPipeline(pepSeq,isProtN,isProtC, ...
            cMgfDatasetIO,strDatasetName,strSpecName,fixedModNameMass, ...
            variableModNameMass,model,method,lambda,ms1_tolerance,ms2_tolerance,...
            alpha,resFilterThres,ionTypes,enzyme,case_penalty_intens,grid_penalty_intens,case_OLS_intens_weight)

            obj.m_pepSeq = pepSeq;
            obj.m_isProtN = isProtN;
            obj.m_isProtC = isProtC;
            obj.m_cMgfDatasetIO = cMgfDatasetIO;
            obj.m_strDatasetName = strDatasetName;
            obj.m_strSpecName = strSpecName;
            obj.m_fixedModNameMass = fixedModNameMass;
            obj.m_variableModNameMass = variableModNameMass;
            obj.m_model = model;
            obj.m_method = method;
            obj.m_lambda = lambda;
            obj.m_ms1_tolerance = ms1_tolerance;
            obj.m_ms2_tolerance = ms2_tolerance;
            obj.m_alpha = alpha;
            obj.m_resFilterThres = resFilterThres;
            obj.m_ionTypes = ionTypes;
            obj.m_enzyme = enzyme;

            if nargin < 18 || isempty(case_penalty_intens)
                obj.m_case_penalty_intens = 'intens_sum';
            else
                obj.m_case_penalty_intens = case_penalty_intens;
            end
            if nargin < 19 || isempty(grid_penalty_intens)
                obj.m_grid_penalty_intens = 'intens_sum';
            else
                obj.m_grid_penalty_intens = grid_penalty_intens;
            end
            if nargin < 20 || isempty(case_OLS_intens_weight)
                obj.m_case_OLS_intens_weight = 'none';
            else
                obj.m_case_OLS_intens_weight = case_OLS_intens_weight;
            end
        end

        [bSuccess,cstrIMP,abundance,ionTypePosCharge,ionIntens,frageff,warning_msg,is_X_not_full_column_rank] = processSpectrum(obj);
    end
end
