classdef CIMPProcessingExecutorConfig
    % Config for CIMPProcessingExecutor

    properties
        ms12DatasetIO
        ms1_tolerance
        minMSMSnum
        alpha
        resFilterThres
        minXicNonzeroPoints
    end

    methods
        function obj = CIMPProcessingExecutorConfig(cfg)
            % Constructor.
            % Input:
            %   cfg (struct)
            %       config struct for IMP processing executor
            if nargin < 1 || ~isstruct(cfg)
                error('CIMPProcessingExecutorConfig:InvalidConstructorArgs', ...
                    'Expected a config struct.');
            end

            cfg = CIMPProcessingExecutorConfig.finalize(cfg);

            obj.ms12DatasetIO = cfg.ms12DatasetIO;
            obj.ms1_tolerance = cfg.ms1_tolerance;
            obj.minMSMSnum = cfg.minMSMSnum;
            obj.alpha = cfg.alpha;
            obj.resFilterThres = cfg.resFilterThres;
            obj.minXicNonzeroPoints = cfg.minXicNonzeroPoints;
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
            if ~isfield(cfg, 'minXicNonzeroPoints') || isempty(cfg.minXicNonzeroPoints)
                cfg.minXicNonzeroPoints = 5;
            end

            if isempty(cfg.ms12DatasetIO)
                error('CIMPProcessingExecutorConfig:InvalidMs12DatasetIO', ...
                    'ms12DatasetIO must be provided.');
            end

            if isempty(cfg.ms1_tolerance) || ~isstruct(cfg.ms1_tolerance) || ...
                    ~isfield(cfg.ms1_tolerance, 'value') || ~isfield(cfg.ms1_tolerance, 'isppm')
                error('CIMPProcessingExecutorConfig:InvalidMs1Tolerance', ...
                    'ms1_tolerance must be a struct with fields: value, isppm.');
            end

            if ~isscalar(cfg.minMSMSnum) || ~isnumeric(cfg.minMSMSnum) || cfg.minMSMSnum < 1
                error('CIMPProcessingExecutorConfig:InvalidMinMSMSnum', ...
                    'minMSMSnum must be a numeric scalar >= 1.');
            end

            if ~isscalar(cfg.alpha) || ~isnumeric(cfg.alpha) || cfg.alpha < 0
                error('CIMPProcessingExecutorConfig:InvalidAlpha', ...
                    'alpha must be a numeric scalar >= 0.');
            end

            if ~isscalar(cfg.resFilterThres) || ~isnumeric(cfg.resFilterThres) || cfg.resFilterThres < 0
                error('CIMPProcessingExecutorConfig:InvalidResFilterThres', ...
                    'resFilterThres must be a numeric scalar >= 0.');
            end
            if ~isscalar(cfg.minXicNonzeroPoints) || ~isnumeric(cfg.minXicNonzeroPoints) || ...
                    ~isfinite(cfg.minXicNonzeroPoints) || cfg.minXicNonzeroPoints < 1 || ...
                    floor(cfg.minXicNonzeroPoints) ~= cfg.minXicNonzeroPoints
                error('CIMPProcessingExecutorConfig:InvalidMinXicNonzeroPoints', ...
                    'minXicNonzeroPoints must be a positive integer.');
            end
        end
    end
end
