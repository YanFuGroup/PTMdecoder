classdef CIMPProcessingPipelineConfig
    % Config for CIMPProcessingPipeline
    % Constructor is kept mainly for backward compatibility.

    properties
        ms12DatasetIO
        ms1_tolerance
        minMSMSnum
        alpha
        resFilterThres
    end

    methods
        function obj = CIMPProcessingPipelineConfig(varargin)
            % Backward-compatible constructor.
            % Input style 1:
            %   CIMPProcessingPipelineConfig(ms12DatasetIO, ms1_tolerance, minMSMSnum, alpha, resFilterThres)
            % Input style 2:
            %   CIMPProcessingPipelineConfig(struct_with_fields)
            if nargin == 1 && isstruct(varargin{1})
                cfg = varargin{1};
            else
                if nargin < 5
                    error('CIMPProcessingPipelineConfig requires 5 arguments or 1 struct argument.');
                end
                cfg = struct(...
                    'ms12DatasetIO', varargin{1}, ...
                    'ms1_tolerance', varargin{2}, ...
                    'minMSMSnum', varargin{3}, ...
                    'alpha', varargin{4}, ...
                    'resFilterThres', varargin{5});
            end

            cfg = CIMPProcessingPipelineConfig.finalize(cfg);

            obj.ms12DatasetIO = cfg.ms12DatasetIO;
            obj.ms1_tolerance = cfg.ms1_tolerance;
            obj.minMSMSnum = cfg.minMSMSnum;
            obj.alpha = cfg.alpha;
            obj.resFilterThres = cfg.resFilterThres;
        end
    end

    methods (Static, Access = private)
        function cfg = finalize(cfg)
            if ~isfield(cfg, 'ms12DatasetIO') || isempty(cfg.ms12DatasetIO)
                cfg.ms12DatasetIO = [];
            end
            if ~isfield(cfg, 'ms1_tolerance') || isempty(cfg.ms1_tolerance)
                cfg.ms1_tolerance = [];
            end
            if ~isfield(cfg, 'minMSMSnum') || isempty(cfg.minMSMSnum)
                cfg.minMSMSnum = 1;
            end
            if ~isfield(cfg, 'alpha') || isempty(cfg.alpha)
                cfg.alpha = 0;
            end
            if ~isfield(cfg, 'resFilterThres') || isempty(cfg.resFilterThres)
                cfg.resFilterThres = 0;
            end

            if isempty(cfg.ms12DatasetIO)
                error('CIMPProcessingPipelineConfig:InvalidMs12DatasetIO', ...
                    'ms12DatasetIO must be provided.');
            end

            if isempty(cfg.ms1_tolerance) || ~isstruct(cfg.ms1_tolerance) || ...
                    ~isfield(cfg.ms1_tolerance, 'value') || ~isfield(cfg.ms1_tolerance, 'isppm')
                error('CIMPProcessingPipelineConfig:InvalidMs1Tolerance', ...
                    'ms1_tolerance must be a struct with fields: value, isppm.');
            end

            if ~isscalar(cfg.minMSMSnum) || ~isnumeric(cfg.minMSMSnum) || cfg.minMSMSnum < 1
                error('CIMPProcessingPipelineConfig:InvalidMinMSMSnum', ...
                    'minMSMSnum must be a numeric scalar >= 1.');
            end

            if ~isscalar(cfg.alpha) || ~isnumeric(cfg.alpha) || cfg.alpha < 0
                error('CIMPProcessingPipelineConfig:InvalidAlpha', ...
                    'alpha must be a numeric scalar >= 0.');
            end

            if ~isscalar(cfg.resFilterThres) || ~isnumeric(cfg.resFilterThres) || cfg.resFilterThres < 0
                error('CIMPProcessingPipelineConfig:InvalidResFilterThres', ...
                    'resFilterThres must be a numeric scalar >= 0.');
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
