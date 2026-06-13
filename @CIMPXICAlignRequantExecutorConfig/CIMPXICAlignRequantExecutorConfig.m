classdef CIMPXICAlignRequantExecutorConfig
    % Config for CIMPXICAlignRequantExecutor

    properties
        ms12DatasetIO
        ms1_tolerance
        minMSMSnum
        resFilterThres
        aligner
        align_strategy
        align_options
    end

    methods
        function obj = CIMPXICAlignRequantExecutorConfig(cfg)
            % Constructor.
            % Input:
            %   cfg (struct)
            %       config struct for XIC align-requant executor
            if nargin < 1 || ~isstruct(cfg)
                CLogger.error(['[CIMPXICAlignRequantExecutorConfig:InvalidConstructorArgs] ', ...
                    'Expected a config struct.']);
            end

            cfg = CIMPXICAlignRequantExecutorConfig.finalize(cfg);

            obj.ms12DatasetIO = cfg.ms12DatasetIO;
            obj.ms1_tolerance = cfg.ms1_tolerance;
            obj.minMSMSnum = cfg.minMSMSnum;
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
                CLogger.error(['[CIMPXICAlignRequantExecutorConfig:InvalidMs12DatasetIO] ', ...
                    'ms12DatasetIO must be provided.']);
            end

            if isempty(cfg.ms1_tolerance) || ~isstruct(cfg.ms1_tolerance) || ...
                    ~isfield(cfg.ms1_tolerance, 'value') || ~isfield(cfg.ms1_tolerance, 'isppm')
                CLogger.error(['[CIMPXICAlignRequantExecutorConfig:InvalidMs1Tolerance] ', ...
                    'ms1_tolerance must be a struct with fields: value, isppm.']);
            end

            if ~isscalar(cfg.minMSMSnum) || ~isnumeric(cfg.minMSMSnum) || cfg.minMSMSnum < 1
                CLogger.error(['[CIMPXICAlignRequantExecutorConfig:InvalidMinMSMSnum] ', ...
                    'minMSMSnum must be a numeric scalar >= 1.']);
            end

            if ~isscalar(cfg.resFilterThres) || ~isnumeric(cfg.resFilterThres) || cfg.resFilterThres < 0
                CLogger.error(['[CIMPXICAlignRequantExecutorConfig:InvalidResFilterThres] ', ...
                    'resFilterThres must be a numeric scalar >= 0.']);
            end

            if isempty(cfg.aligner)
                CLogger.error(['[CIMPXICAlignRequantExecutorConfig:InvalidAligner] ', ...
                    'aligner must be provided.']);
            end

            if isempty(cfg.align_strategy)
                CLogger.error(['[CIMPXICAlignRequantExecutorConfig:InvalidAlignStrategy] ', ...
                    'align_strategy must be provided.']);
            end

            if ~isstruct(cfg.align_options)
                CLogger.error(['[CIMPXICAlignRequantExecutorConfig:InvalidAlignOptions] ', ...
                    'align_options must be a struct.']);
            end
        end
    end
end
