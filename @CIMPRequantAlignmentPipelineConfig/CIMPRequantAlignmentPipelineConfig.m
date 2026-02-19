classdef CIMPRequantAlignmentPipelineConfig
    % Config for CIMPRequantAlignmentPipeline
    % Recommended external usage:
    %   cfg = CIMPRequantAlignmentPipelineConfig.fromTaskParam(...)
    % Constructor is kept mainly for backward compatibility.

    properties
        ms12DatasetIO
        ms1_tolerance
        minMSMSnum
        alpha
        resFilterThres
        aligner
        align_strategy
        align_options
    end

    methods
        function obj = CIMPRequantAlignmentPipelineConfig(varargin)
            % Backward-compatible constructor.
            % Input style 1:
            %   CIMPRequantAlignmentPipelineConfig(ms12DatasetIO, ms1_tolerance, minMSMSnum, alpha, resFilterThres, aligner, align_strategy, align_options)
            % Input style 2:
            %   CIMPRequantAlignmentPipelineConfig(struct_with_fields)
            if nargin == 1 && isstruct(varargin{1})
                cfg = varargin{1};
            else
                if nargin < 7
                    error('CIMPRequantAlignmentPipelineConfig requires at least 7 arguments or 1 struct argument.');
                end
                if nargin < 8
                    align_options = struct();
                else
                    align_options = varargin{8};
                end
                cfg = struct(...
                    'ms12DatasetIO', varargin{1}, ...
                    'ms1_tolerance', varargin{2}, ...
                    'minMSMSnum', varargin{3}, ...
                    'alpha', varargin{4}, ...
                    'resFilterThres', varargin{5}, ...
                    'aligner', varargin{6}, ...
                    'align_strategy', varargin{7}, ...
                    'align_options', align_options);
            end

            cfg = CIMPRequantAlignmentPipelineConfig.finalize(cfg);

            obj.ms12DatasetIO = cfg.ms12DatasetIO;
            obj.ms1_tolerance = cfg.ms1_tolerance;
            obj.minMSMSnum = cfg.minMSMSnum;
            obj.alpha = cfg.alpha;
            obj.resFilterThres = cfg.resFilterThres;
            obj.aligner = cfg.aligner;
            obj.align_strategy = cfg.align_strategy;
            obj.align_options = cfg.align_options;
        end
    end

    methods (Static)
        function obj = fromTaskParam(taskParam, ms12DatasetIO, aligner, align_strategy, align_options, overrides)
            % Preferred factory method for business-layer calls.
            % Build config from CTaskParam with optional overrides
            % Input:
            %   taskParam (CTaskParam)
            %       task parameter object
            %   ms12DatasetIO (CMS12DatasetIO)
            %       initialized MS1/MS2 dataset IO
            %   aligner (CXICAligner)
            %       aligner instance
            %   align_strategy (CRunAlignStrategy)
            %       alignment strategy
            %   align_options (struct, optional)
            %       alignment options
            %   overrides (struct, optional)
            %       override fields of this config

            if nargin < 4
                error('CIMPRequantAlignmentPipelineConfig:MissingInput', ...
                    'taskParam, ms12DatasetIO, aligner, and align_strategy are required.');
            end
            if nargin < 5 || isempty(align_options)
                align_options = struct();
            end
            if nargin < 6 || isempty(overrides)
                overrides = struct();
            end

            cfg = struct();
            cfg.ms12DatasetIO = ms12DatasetIO;
            cfg.ms1_tolerance = taskParam.m_ms1_tolerance;
            cfg.minMSMSnum = taskParam.m_min_MSMS_num;
            cfg.alpha = taskParam.m_alpha;
            cfg.resFilterThres = taskParam.m_result_filter_threshold;
            cfg.aligner = aligner;
            cfg.align_strategy = align_strategy;
            cfg.align_options = align_options;

            cfg = CIMPRequantAlignmentPipelineConfig.applyOverrides(cfg, overrides);
            obj = CIMPRequantAlignmentPipelineConfig(cfg);
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
            if ~isfield(cfg, 'aligner') || isempty(cfg.aligner)
                cfg.aligner = [];
            end
            if ~isfield(cfg, 'align_strategy') || isempty(cfg.align_strategy)
                cfg.align_strategy = [];
            end
            if ~isfield(cfg, 'align_options') || isempty(cfg.align_options)
                cfg.align_options = struct();
            end

            if isempty(cfg.ms12DatasetIO)
                error('CIMPRequantAlignmentPipelineConfig:InvalidMs12DatasetIO', ...
                    'ms12DatasetIO must be provided.');
            end

            if isempty(cfg.ms1_tolerance) || ~isstruct(cfg.ms1_tolerance) || ...
                    ~isfield(cfg.ms1_tolerance, 'value') || ~isfield(cfg.ms1_tolerance, 'isppm')
                error('CIMPRequantAlignmentPipelineConfig:InvalidMs1Tolerance', ...
                    'ms1_tolerance must be a struct with fields: value, isppm.');
            end

            if ~isscalar(cfg.minMSMSnum) || ~isnumeric(cfg.minMSMSnum) || cfg.minMSMSnum < 1
                error('CIMPRequantAlignmentPipelineConfig:InvalidMinMSMSnum', ...
                    'minMSMSnum must be a numeric scalar >= 1.');
            end

            if ~isscalar(cfg.alpha) || ~isnumeric(cfg.alpha) || cfg.alpha < 0
                error('CIMPRequantAlignmentPipelineConfig:InvalidAlpha', ...
                    'alpha must be a numeric scalar >= 0.');
            end

            if ~isscalar(cfg.resFilterThres) || ~isnumeric(cfg.resFilterThres) || cfg.resFilterThres < 0
                error('CIMPRequantAlignmentPipelineConfig:InvalidResFilterThres', ...
                    'resFilterThres must be a numeric scalar >= 0.');
            end

            if isempty(cfg.aligner)
                error('CIMPRequantAlignmentPipelineConfig:InvalidAligner', ...
                    'aligner must be provided.');
            end

            if isempty(cfg.align_strategy)
                error('CIMPRequantAlignmentPipelineConfig:InvalidAlignStrategy', ...
                    'align_strategy must be provided.');
            end

            if ~isstruct(cfg.align_options)
                error('CIMPRequantAlignmentPipelineConfig:InvalidAlignOptions', ...
                    'align_options must be a struct.');
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
