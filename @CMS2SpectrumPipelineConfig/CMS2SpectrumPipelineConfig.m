classdef CMS2SpectrumPipelineConfig
    % Config for CMS2SpectrumPipeline
    % Recommended external usage:
    %   cfg = CMS2SpectrumPipelineConfig.fromTaskParam(...)
    % Constructor supports struct and positional argument styles.

    properties
        model
        method
        lambda
        ms1_tolerance
        ms2_tolerance
        alpha
        resFilterThres
        ionTypes
        enzyme
        case_penalty_intens
        grid_penalty_intens
        case_OLS_intens_weight
    end

    methods
        function obj = CMS2SpectrumPipelineConfig(varargin)
            % Input style 1:
            %   CMS2SpectrumPipelineConfig(struct_with_fields)
            % Input style 2:
            %   CMS2SpectrumPipelineConfig(model, method, lambda, ms1_tolerance, ms2_tolerance, ...
            %       alpha, resFilterThres, ionTypes, enzyme, case_penalty_intens, ...
            %       grid_penalty_intens, case_OLS_intens_weight)

            if nargin == 1 && isstruct(varargin{1})
                cfg = varargin{1};
            else
                if nargin < 9
                    error('CMS2SpectrumPipelineConfig requires at least 9 arguments or 1 struct argument.');
                end
                cfg = struct(...
                    'model', varargin{1}, ...
                    'method', varargin{2}, ...
                    'lambda', varargin{3}, ...
                    'ms1_tolerance', varargin{4}, ...
                    'ms2_tolerance', varargin{5}, ...
                    'alpha', varargin{6}, ...
                    'resFilterThres', varargin{7}, ...
                    'ionTypes', varargin{8}, ...
                    'enzyme', varargin{9});

                if nargin >= 10
                    cfg.case_penalty_intens = varargin{10};
                end
                if nargin >= 11
                    cfg.grid_penalty_intens = varargin{11};
                end
                if nargin >= 12
                    cfg.case_OLS_intens_weight = varargin{12};
                end
            end

            cfg = CMS2SpectrumPipelineConfig.finalize(cfg);

            obj.model = cfg.model;
            obj.method = cfg.method;
            obj.lambda = cfg.lambda;
            obj.ms1_tolerance = cfg.ms1_tolerance;
            obj.ms2_tolerance = cfg.ms2_tolerance;
            obj.alpha = cfg.alpha;
            obj.resFilterThres = cfg.resFilterThres;
            obj.ionTypes = cfg.ionTypes;
            obj.enzyme = cfg.enzyme;
            obj.case_penalty_intens = cfg.case_penalty_intens;
            obj.grid_penalty_intens = cfg.grid_penalty_intens;
            obj.case_OLS_intens_weight = cfg.case_OLS_intens_weight;
        end
    end

    methods (Static)
        function obj = fromTaskParam(taskParam, overrides)
            % Build config from CTaskParam with optional overrides.
            % Input:
            %   taskParam (CTaskParam)
            %   overrides (struct, optional)

            if nargin < 1 || isempty(taskParam)
                error('CMS2SpectrumPipelineConfig:MissingInput', 'taskParam is required.');
            end
            if nargin < 2 || isempty(overrides)
                overrides = struct();
            end

            % Keep explicit finalize() validation even though values come from CTaskParam.
            % This provides a single defensive gate for runtime safety now; for future C++
            % migration we can redesign this validation boundary with stronger typed configs.

            cfg = struct();
            cfg.model = taskParam.m_model;
            cfg.method = taskParam.m_method;
            cfg.lambda = taskParam.m_lambda;
            cfg.ms1_tolerance = taskParam.m_ms1_tolerance;
            cfg.ms2_tolerance = taskParam.m_ms2_tolerance;
            cfg.alpha = taskParam.m_alpha;
            cfg.resFilterThres = taskParam.m_result_filter_threshold;

            % The following are pipeline defaults (not part of current CTaskParam contract)
            cfg.ionTypes = [1,2];
            cfg.enzyme = struct('name', taskParam.m_enzyme_name, 'limits', taskParam.m_enzyme_limits);
            cfg.case_penalty_intens = 'intens_sum';
            cfg.grid_penalty_intens = 'intens_sum';
            cfg.case_OLS_intens_weight = 'none';

            cfg = CMS2SpectrumPipelineConfig.applyOverrides(cfg, overrides);
            obj = CMS2SpectrumPipelineConfig(cfg);
        end
    end

    methods (Static, Access = private)
        function cfg = finalize(cfg)
            if ~isfield(cfg, 'model') || isempty(cfg.model)
                cfg.model = [];
            end
            if ~isfield(cfg, 'method') || isempty(cfg.method)
                cfg.method = [];
            end
            if ~isfield(cfg, 'lambda') || isempty(cfg.lambda)
                cfg.lambda = 0;
            end
            if ~isfield(cfg, 'ms1_tolerance') || isempty(cfg.ms1_tolerance)
                cfg.ms1_tolerance = [];
            end
            if ~isfield(cfg, 'ms2_tolerance') || isempty(cfg.ms2_tolerance)
                cfg.ms2_tolerance = [];
            end
            if ~isfield(cfg, 'alpha') || isempty(cfg.alpha)
                cfg.alpha = 0;
            end
            if ~isfield(cfg, 'resFilterThres') || isempty(cfg.resFilterThres)
                cfg.resFilterThres = 0;
            end
            if ~isfield(cfg, 'ionTypes') || isempty(cfg.ionTypes)
                cfg.ionTypes = [1,2];
            end
            if ~isfield(cfg, 'enzyme') || isempty(cfg.enzyme)
                cfg.enzyme = [];
            end
            if ~isfield(cfg, 'case_penalty_intens') || isempty(cfg.case_penalty_intens)
                cfg.case_penalty_intens = 'intens_sum';
            end
            if ~isfield(cfg, 'grid_penalty_intens') || isempty(cfg.grid_penalty_intens)
                cfg.grid_penalty_intens = 'intens_sum';
            end
            if ~isfield(cfg, 'case_OLS_intens_weight') || isempty(cfg.case_OLS_intens_weight)
                cfg.case_OLS_intens_weight = 'none';
            end

            if isempty(cfg.model)
                error('CMS2SpectrumPipelineConfig:InvalidModel', 'model must be provided.');
            end
            if isempty(cfg.method)
                error('CMS2SpectrumPipelineConfig:InvalidMethod', 'method must be provided.');
            end
            if ~isscalar(cfg.lambda) || ~isnumeric(cfg.lambda) || cfg.lambda < 0
                error('CMS2SpectrumPipelineConfig:InvalidLambda', 'lambda must be a numeric scalar >= 0.');
            end
            if isempty(cfg.ms1_tolerance) || ~isstruct(cfg.ms1_tolerance) || ...
                    ~isfield(cfg.ms1_tolerance, 'value') || ~isfield(cfg.ms1_tolerance, 'isppm')
                error('CMS2SpectrumPipelineConfig:InvalidMs1Tolerance', ...
                    'ms1_tolerance must be a struct with fields: value, isppm.');
            end
            if isempty(cfg.ms2_tolerance) || ~isnumeric(cfg.ms2_tolerance) 
                error('CMS2SpectrumPipelineConfig:InvalidMs2Tolerance', ...
                    'ms2_tolerance must be a numeric scalar >= 0.');
            end
            if ~isscalar(cfg.alpha) || ~isnumeric(cfg.alpha) || cfg.alpha < 0
                error('CMS2SpectrumPipelineConfig:InvalidAlpha', 'alpha must be a numeric scalar >= 0.');
            end
            if ~isscalar(cfg.resFilterThres) || ~isnumeric(cfg.resFilterThres) || cfg.resFilterThres < 0
                error('CMS2SpectrumPipelineConfig:InvalidResFilterThres', ...
                    'resFilterThres must be a numeric scalar >= 0.');
            end
            if isempty(cfg.ionTypes) || ~isnumeric(cfg.ionTypes)
                error('CMS2SpectrumPipelineConfig:InvalidIonTypes', 'ionTypes must be a numeric array.');
            end
            if isempty(cfg.enzyme) || ~isstruct(cfg.enzyme) || ...
                    ~isfield(cfg.enzyme, 'name') || ~isfield(cfg.enzyme, 'limits')
                error('CMS2SpectrumPipelineConfig:InvalidEnzyme', ...
                    'enzyme must be a struct with fields: name, limits.');
            end
        end

        function cfg = applyOverrides(cfg, overrides)
            override_fields = fieldnames(overrides);
            for idx = 1:numel(override_fields)
                cfg.(override_fields{idx}) = overrides.(override_fields{idx});
            end
        end
    end
end
