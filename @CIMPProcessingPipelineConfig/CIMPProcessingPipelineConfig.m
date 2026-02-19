classdef CIMPProcessingPipelineConfig
    % Config for CIMPProcessingPipeline
    % Recommended external usage:
    %   cfg = CIMPProcessingPipelineConfig.fromTaskParam(...)
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

    methods (Static)
        function obj = fromTaskParam(taskParam, ms12DatasetIO, overrides)
            % Preferred factory method for business-layer calls.
            % Build config from CTaskParam with optional overrides
            % Input:
            %   taskParam (CTaskParam)
            %       task parameter object
            %   ms12DatasetIO (CMS12DatasetIO)
            %       initialized MS1/MS2 dataset IO
            %   overrides (struct, optional)
            %       override fields: minMSMSnum, alpha, resFilterThres, ms1_tolerance, ms12DatasetIO

            if nargin < 2
                error('CIMPProcessingPipelineConfig:MissingInput', ...
                    'taskParam and ms12DatasetIO are required.');
            end
            if nargin < 3 || isempty(overrides)
                overrides = struct();
            end

            cfg = struct();
            cfg.ms12DatasetIO = ms12DatasetIO;
            cfg.ms1_tolerance = taskParam.m_ms1_tolerance;
            cfg.minMSMSnum = taskParam.m_min_MSMS_num;
            cfg.alpha = taskParam.m_alpha;
            cfg.resFilterThres = taskParam.m_result_filter_threshold;

            cfg = CIMPProcessingPipelineConfig.applyOverrides(cfg, overrides);
            obj = CIMPProcessingPipelineConfig(cfg);
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
