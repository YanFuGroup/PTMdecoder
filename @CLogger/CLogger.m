classdef CLogger
    % Lightweight logger facade for workflow/service diagnostics.

    methods (Static)
        function cfg = getConfig()
            cfg = CLogger.getCore().getConfig();
        end

        function configure(cfg)
            core = CLogger.getCore();
            core.configure(cfg);
        end

        function flush()
            CLogger.getCore().flush();
        end

        function debug(msg, varargin)
            CLogger.getCore().log('DEBUG', msg, varargin{:});
        end

        function info(msg, varargin)
            CLogger.getCore().log('INFO', msg, varargin{:});
        end

        function warn(msg, varargin)
            CLogger.getCore().log('WARN', msg, varargin{:});
        end

        function error(msg, varargin)
            % Log and then raise an exception. Do not swallow errors.
            % Input patterns:
            %   1) error(msgFmt, args...)
            %      - logs formatted message and raises CLogger:LoggedError
            %   2) error(me, contextFmt, args...)
            %      - logs context + original exception, then rethrows original me
            if isa(msg, 'MException')
                me = msg;
                if ~isempty(varargin)
                    CLogger.getCore().log('ERROR', varargin{1}, varargin{2:end});
                end
                CLogger.getCore().log('ERROR', 'Original exception: [%s] %s', me.identifier, me.message);
                throwAsCaller(me);
            else
                CLogger.getCore().log('ERROR', msg, varargin{:});
                raised_message = sprintf(msg, varargin{:});
                builtin('error', 'CLogger:LoggedError', '%s', raised_message);
            end
        end

        function progress(label, current_step, total_step)
            % Render in-place progress on console without writing percentage logs.
            CLogger.getCore().progress(label, current_step, total_step);
        end

        function resetForTests()
            CLogger.getCore().resetForTests();
        end
    end

    methods (Static, Access = private)
        function core = getCore()
            persistent logger_core;
            if isempty(logger_core) || ~isvalid(logger_core)
                logger_core = CLoggerCore();
            end
            core = logger_core;
        end
    end
end
