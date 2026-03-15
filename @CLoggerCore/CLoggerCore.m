classdef CLoggerCore < handle
    % Stateful logger core with buffered file writes.

    properties (Access = private)
        enabled
        file_level
        file_level_rank
        console_level
        console_level_rank
        to_console
        file_path
        buffer
        buffer_size
        last_progress_text_length
    end

    methods
        function obj = CLoggerCore()
            obj.resetDefaults();
        end

        function cfg = getConfig(obj)
            % Retrieve current runtime settings.
            cfg = struct();
            cfg.enabled = obj.enabled;
            cfg.file_level = obj.file_level;
            cfg.console_level = obj.console_level;
            cfg.to_console = obj.to_console;
            cfg.file_path = obj.file_path;
            cfg.buffer_size = obj.buffer_size;
        end

        function configure(obj, cfg)
            % Configure logger runtime settings.
            % Input:
            %   cfg (struct)
            %       optional logger config with fields:
            %       enabled, file_level, console_level, to_console, file_path, buffer_size
            obj.flush();
            obj.resetDefaults();

            if nargin >= 2 && ~isempty(cfg)
                if isfield(cfg, 'enabled')
                    obj.enabled = obj.ensureLogicalScalar(cfg.enabled, 'enabled');
                end
                if isfield(cfg, 'file_level') && ~isempty(cfg.file_level)
                    obj.file_level = upper(strtrim(char(cfg.file_level)));
                end
                if isfield(cfg, 'file_path') && ~isempty(cfg.file_path)
                    obj.file_path = char(cfg.file_path);
                end
                if isfield(cfg, 'to_console')
                    obj.to_console = obj.ensureLogicalScalar(cfg.to_console, 'to_console');
                end
                if isfield(cfg, 'console_level') && ~isempty(cfg.console_level)
                    obj.console_level = upper(strtrim(char(cfg.console_level)));
                end
                if isfield(cfg, 'buffer_size') && ~isempty(cfg.buffer_size)
                    obj.buffer_size = obj.toPositiveInt(cfg.buffer_size, obj.buffer_size);
                end
            end

            if isempty(obj.file_path)
                obj.file_path = obj.defaultLogFilePath();
            end

            obj.file_level_rank = obj.levelRank(obj.file_level);
            obj.console_level_rank = obj.levelRank(obj.console_level);
            obj.ensureParentDir(obj.file_path);
            obj.log('INFO', 'Logger configured. file_level=%s, console_level=%s, file=%s', ...
                obj.file_level, obj.console_level, obj.file_path);
        end

        function log(obj, level, msg, varargin)
            if ~obj.enabled
                return;
            end
            level_rank = obj.levelRank(level);

            message = sprintf(msg, varargin{:});
            ts = datestr(now, 'yyyy-mm-dd HH:MM:SS.FFF');
            line = sprintf('[%s] [%s] %s', ts, upper(level), message);

            if level_rank >= obj.file_level_rank
                obj.buffer{end + 1} = line;
                if numel(obj.buffer) >= obj.buffer_size || strcmpi(level, 'ERROR')
                    obj.flush();
                end
            end
            if obj.to_console && level_rank >= obj.console_level_rank
                obj.endProgressLineIfNeeded();
                fprintf('[%s] %s\n', upper(level), message);
            end
        end

        function progress(obj, label, current_step, total_step)
            % Render progress in-place on console without writing to log file.
            if ~obj.enabled || ~obj.to_console
                return;
            end
            if nargin < 4 || isempty(total_step) || total_step <= 0
                total_step = 1;
            end
            if nargin < 3 || isempty(current_step)
                current_step = 0;
            end
            if nargin < 2 || isempty(label)
                label = 'progress';
            end

            pct = floor(double(current_step) / double(total_step) * 100);
            pct = max(0, min(100, pct));
            progress_text = sprintf('%s: %d/%d (%d%%)', char(label), current_step, total_step, pct);

            if obj.last_progress_text_length > 0
                fprintf(repmat('\b', 1, obj.last_progress_text_length));
            end

            fprintf('%s', progress_text);

            pad_len = max(0, obj.last_progress_text_length - length(progress_text));
            if pad_len > 0
                fprintf(repmat(' ', 1, pad_len));
                fprintf(repmat('\b', 1, pad_len));
            end

            obj.last_progress_text_length = length(progress_text);

            if current_step >= total_step
                fprintf('\n');
                obj.last_progress_text_length = 0;
            end
        end

        function flush(obj)
            % Force flush buffered file logs to disk.
            if isempty(obj.buffer)
                return;
            end
            obj.appendLines(obj.file_path, obj.buffer);
            obj.buffer = {};
        end

        function resetForTests(obj)
            obj.flush();
            obj.resetDefaults();
        end
    end

    methods (Access = private)
        function appendLines(~, path, lines)
            if isempty(lines)
                return;
            end
            fid = fopen(path, 'a');
            if fid <= 0
                return;
            end
            for idx_line = 1:numel(lines)
                fprintf(fid, '%s\n', lines{idx_line});
            end
            fclose(fid);
        end

        function ensureParentDir(~, path)
            [folder, ~, ~] = fileparts(path);
            if isempty(folder)
                return;
            end
            if ~exist(folder, 'dir')
                mkdir(folder);
            end
        end

        function rank = levelRank(~, level)
            switch upper(level)
                case 'DEBUG'
                    rank = 0;
                case 'INFO'
                    rank = 1;
                case 'WARN'
                    rank = 2;
                case 'ERROR'
                    rank = 3;
                otherwise
                    rank = 1;
            end
        end

        function value = toPositiveInt(~, raw, default_value)
            if nargin < 3
                default_value = 50;
            end

            if isnumeric(raw)
                parsed = double(raw);
            elseif isstring(raw) && isscalar(raw)
                parsed = str2double(char(raw));
            elseif ischar(raw)
                parsed = str2double(raw);
            else
                value = default_value;
                return;
            end

            if isnan(parsed) || ~isfinite(parsed) || parsed < 1
                value = default_value;
                return;
            end

            value = floor(parsed);
        end

        function value = ensureLogicalScalar(~, raw, field_name)
            if ~islogical(raw) || ~isscalar(raw)
                error('CLoggerCore:InvalidLogicalOption', ...
                    'Logger option ''%s'' must be a logical scalar.', field_name);
            end
            value = logical(raw);
        end

        function resetDefaults(obj)
            obj.enabled = true;
            obj.file_level = 'DEBUG';
            obj.file_level_rank = 0;
            obj.console_level = 'INFO';
            obj.console_level_rank = 1;
            obj.to_console = true;
            obj.file_path = obj.defaultLogFilePath();
            obj.buffer = {};
            obj.buffer_size = 50;
            obj.last_progress_text_length = 0;
        end

        function endProgressLineIfNeeded(obj)
            if obj.last_progress_text_length > 0
                fprintf('\n');
                obj.last_progress_text_length = 0;
            end
        end

        function path = defaultLogFilePath(~)
            % Build a timestamped default log file path in the working directory.
            timestamp = datestr(now, 'yyyymmdd_HHMMSS');
            path = fullfile(pwd, sprintf('ptmdecoder_%s.log', timestamp));
        end
    end
end
