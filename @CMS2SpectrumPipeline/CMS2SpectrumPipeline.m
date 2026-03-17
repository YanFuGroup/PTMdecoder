classdef CMS2SpectrumPipeline
    % Pipeline for MS2 spectrum-level identification and quantification

    properties
        m_fixedModNameMass;
        m_variableModNameMass;
        m_model;
        m_method;
        m_lambda;
        m_ms1_tolerance;
        m_ms2_tolerance;
        m_alpha;
        m_resFilterThres;
        m_max_mod_per_peptide;
        m_ionTypes;
        m_enzyme;

        m_case_penalty_intens;
        m_grid_penalty_intens;
        m_case_OLS_intens_weight;
    end

    methods
        function obj = CMS2SpectrumPipeline(fixedModNameMass, variableModNameMass, cfg)
            if nargin ~= 3
                error('CMS2SpectrumPipeline:InvalidConstructorInput', ...
                    'Constructor requires fixedModNameMass, variableModNameMass, and cfg.');
            end
            if ~isa(cfg, 'CMS2SpectrumPipelineConfig')
                error('CMS2SpectrumPipeline:InvalidConfig', ...
                    'cfg must be a CMS2SpectrumPipelineConfig instance.');
            end

            obj.m_fixedModNameMass = fixedModNameMass;
            obj.m_variableModNameMass = variableModNameMass;
            obj.m_model = cfg.model;
            obj.m_method = cfg.method;
            obj.m_lambda = cfg.lambda;
            obj.m_ms1_tolerance = cfg.ms1_tolerance;
            obj.m_ms2_tolerance = cfg.ms2_tolerance;
            obj.m_alpha = cfg.alpha;
            obj.m_resFilterThres = cfg.resFilterThres;
            obj.m_max_mod_per_peptide = cfg.max_mod_per_peptide;
            obj.m_ionTypes = cfg.ionTypes;
            obj.m_enzyme = cfg.enzyme;
            obj.m_case_penalty_intens = cfg.case_penalty_intens;
            obj.m_grid_penalty_intens = cfg.grid_penalty_intens;
            obj.m_case_OLS_intens_weight = cfg.case_OLS_intens_weight;
        end

        [bSuccess,cstrIMP,abundance,ionTypePosCharge,ionIntens,frageff,is_X_not_full_column_rank,solver_diag,noise_model_fit_inputs,stability_cache] = runBaselineSpectrumStage(obj, peptideCtx, spectrumCtx);
    end
end
