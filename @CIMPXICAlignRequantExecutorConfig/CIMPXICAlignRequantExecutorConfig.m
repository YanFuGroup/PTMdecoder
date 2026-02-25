classdef CIMPXICAlignRequantExecutorConfig
    % Config for CIMPXICAlignRequantExecutor

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
        function obj = CIMPXICAlignRequantExecutorConfig(varargin)
            % Constructor.
            % Input style 1:
            %   CIMPXICAlignRequantExecutorConfig(ms12DatasetIO, ms1_tolerance, minMSMSnum, alpha, resFilterThres, aligner, align_strategy, align_options)
            % Input style 2:
            %   CIMPXICAlignRequantExecutorConfig(struct_with_fields)
            if nargin == 1 && isstruct(varargin{1})
                cfg = varargin{1};
            else
                if nargin < 7
                    error('CIMPXICAlignRequantExecutorConfig requires at least 7 arguments or 1 struct argument.');
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

            cfg = CIMPXICAlignRequantExecutorConfig.finalize(cfg);

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
                error('CIMPXICAlignRequantExecutorConfig:InvalidMs12DatasetIO', ...
                    'ms12DatasetIO must be provided.');
            end

            if isempty(cfg.ms1_tolerance) || ~isstruct(cfg.ms1_tolerance) || ...
                    ~isfield(cfg.ms1_tolerance, 'value') || ~isfield(cfg.ms1_tolerance, 'isppm')
                error('CIMPXICAlignRequantExecutorConfig:InvalidMs1Tolerance', ...
                    'ms1_tolerance must be a struct with fields: value, isppm.');
            end

            if ~isscalar(cfg.minMSMSnum) || ~isnumeric(cfg.minMSMSnum) || cfg.minMSMSnum < 1
                error('CIMPXICAlignRequantExecutorConfig:InvalidMinMSMSnum', ...
                    'minMSMSnum must be a numeric scalar >= 1.');
            end

            if ~isscalar(cfg.alpha) || ~isnumeric(cfg.alpha) || cfg.alpha < 0
                error('CIMPXICAlignRequantExecutorConfig:InvalidAlpha', ...
                    'alpha must be a numeric scalar >= 0.');
            end

            if ~isscalar(cfg.resFilterThres) || ~isnumeric(cfg.resFilterThres) || cfg.resFilterThres < 0
                error('CIMPXICAlignRequantExecutorConfig:InvalidResFilterThres', ...
                    'resFilterThres must be a numeric scalar >= 0.');
            end

            if isempty(cfg.aligner)
                error('CIMPXICAlignRequantExecutorConfig:InvalidAligner', ...
                    'aligner must be provided.');
            end

            if isempty(cfg.align_strategy)
                error('CIMPXICAlignRequantExecutorConfig:InvalidAlignStrategy', ...
                    'align_strategy must be provided.');
            end

            if ~isstruct(cfg.align_options)
                error('CIMPXICAlignRequantExecutorConfig:InvalidAlignOptions', ...
                    'align_options must be a struct.');
            end
        end
    end
end
