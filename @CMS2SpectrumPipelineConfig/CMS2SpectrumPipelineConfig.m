classdef CMS2SpectrumPipelineConfig
    % Config for CMS2SpectrumPipeline
    % Constructor accepts struct input style.

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
        function obj = CMS2SpectrumPipelineConfig(cfg)
            % Input:
            %   cfg (struct)
            %       config struct with CMS2 spectrum pipeline fields
            if nargin < 1 || ~isstruct(cfg)
                error('CMS2SpectrumPipelineConfig:InvalidConstructorArgs', ...
                    'Expected a config struct.');
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

    end
end
