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
            variableModNameMass,varargin)

            if nargin == 9 && isa(varargin{1}, 'CMS2SpectrumPipelineConfig')
                cfg = varargin{1};

                model = cfg.model;
                method = cfg.method;
                lambda = cfg.lambda;
                ms1_tolerance = cfg.ms1_tolerance;
                ms2_tolerance = cfg.ms2_tolerance;
                alpha = cfg.alpha;
                resFilterThres = cfg.resFilterThres;
                ionTypes = cfg.ionTypes;
                enzyme = cfg.enzyme;
                case_penalty_intens = cfg.case_penalty_intens;
                grid_penalty_intens = cfg.grid_penalty_intens;
                case_OLS_intens_weight = cfg.case_OLS_intens_weight;
            else
                if numel(varargin) < 9
                    error('CMS2SpectrumPipeline requires either cfg-mode or full-argument mode.');
                end

                model = varargin{1};
                method = varargin{2};
                lambda = varargin{3};
                ms1_tolerance = varargin{4};
                ms2_tolerance = varargin{5};
                alpha = varargin{6};
                resFilterThres = varargin{7};
                ionTypes = varargin{8};
                enzyme = varargin{9};

                if numel(varargin) >= 10 && ~isempty(varargin{10})
                    case_penalty_intens = varargin{10};
                else
                    case_penalty_intens = 'intens_sum';
                end
                if numel(varargin) >= 11 && ~isempty(varargin{11})
                    grid_penalty_intens = varargin{11};
                else
                    grid_penalty_intens = 'intens_sum';
                end
                if numel(varargin) >= 12 && ~isempty(varargin{12})
                    case_OLS_intens_weight = varargin{12};
                else
                    case_OLS_intens_weight = 'none';
                end
            end

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
            obj.m_case_penalty_intens = case_penalty_intens;
            obj.m_grid_penalty_intens = grid_penalty_intens;
            obj.m_case_OLS_intens_weight = case_OLS_intens_weight;
        end

        [bSuccess,cstrIMP,abundance,ionTypePosCharge,ionIntens,frageff,is_X_not_full_column_rank] = processSpectrumWithContext(obj, peptideCtx, spectrumCtx);
    end
end
