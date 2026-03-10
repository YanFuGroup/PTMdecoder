classdef CParamMapUtils
    % Utilities for containers.Map-based parameter retrieval.

    methods (Static)
        function value = getRequired(param_map, key_name, explain, err_id_prefix)
            % Get a required parameter value from map, error if missing.
            % Input:
            %   param_map (containers.Map)
            %       parameter key-value map
            %   key_name (1 x 1 char/string)
            %       parameter key name
            %   explain (1 x 1 char/string)
            %       description for error message
            %   err_id_prefix (1 x 1 char/string, optional)
            %       error id prefix (default: CParamMapUtils)
            % Output:
            %   value (any)
            %       required parameter value
            if nargin < 4 || isempty(err_id_prefix)
                err_id_prefix = 'CParamMapUtils';
            end
            if ~isa(param_map, 'containers.Map')
                error([err_id_prefix, ':InvalidParamMap'], ...
                    'Expected param_map to be a containers.Map.');
            end
            if ~param_map.isKey(key_name)
                error([err_id_prefix, ':MissingRequiredParam'], ...
                    'Required param ''%s'' is missing (%s).', key_name, explain);
            end
            value = param_map(key_name);
        end

        function value = getOptional(param_map, key_name, default_value, err_id_prefix)
            % Get an optional parameter value from map, use default if missing.
            % Input:
            %   param_map (containers.Map)
            %       parameter key-value map
            %   key_name (1 x 1 char/string)
            %       parameter key name
            %   default_value (any, optional)
            %       default value if key not present;
            %       if omitted, returns [] when key missing/empty
            %   err_id_prefix (1 x 1 char/string, optional)
            %       error id prefix (default: CParamMapUtils)
            % Output:
            %   value (any)
            %       parameter value, default value, or [] if no default
            if nargin < 3
                default_value = [];
            end
            if nargin < 4 || isempty(err_id_prefix)
                err_id_prefix = 'CParamMapUtils';
            end
            if ~isa(param_map, 'containers.Map')
                error([err_id_prefix, ':InvalidParamMap'], ...
                    'Expected param_map to be a containers.Map.');
            end
            if param_map.isKey(key_name) && ~isempty(param_map(key_name))
                value = param_map(key_name);
            else
                value = default_value;
            end
        end

        function value = getRequiredNumber(param_map, key_name, explain, err_id_prefix)
            % Get a required numeric parameter value from map.
            % Input:
            %   param_map (containers.Map)
            %       parameter key-value map
            %   key_name (1 x 1 char/string)
            %       parameter key name
            %   explain (1 x 1 char/string)
            %       description for error message
            %   err_id_prefix (1 x 1 char/string, optional)
            %       error id prefix (default: CParamMapUtils)
            % Output:
            %   value (double)
            %       required numeric parameter value
            if nargin < 4 || isempty(err_id_prefix)
                err_id_prefix = 'CParamMapUtils';
            end
            raw_value = CParamMapUtils.getRequired(param_map, key_name, explain, err_id_prefix);
            value = CParamMapUtils.toNumber(raw_value, key_name, err_id_prefix);
        end

        function value = getOptionalNumber(param_map, key_name, default_value, err_id_prefix)
            % Get an optional numeric parameter value from map.
            % Input:
            %   param_map (containers.Map)
            %       parameter key-value map
            %   key_name (1 x 1 char/string)
            %       parameter key name
            %   default_value (numeric/char/string, optional)
            %       default numeric value if key not present;
            %       if omitted, returns [] when key missing/empty
            %   err_id_prefix (1 x 1 char/string, optional)
            %       error id prefix (default: CParamMapUtils)
            % Output:
            %   value (double or [])
            %       optional numeric parameter value, or [] if no default
            if nargin < 4 || isempty(err_id_prefix)
                err_id_prefix = 'CParamMapUtils';
            end
            if nargin < 3
                default_value = [];
            end

            raw_value = CParamMapUtils.getOptional(param_map, key_name, default_value, err_id_prefix);
            if isempty(raw_value)
                value = [];
                return;
            end

            value = CParamMapUtils.toNumber(raw_value, key_name, err_id_prefix);
        end

        function value = getOptionalLogical(param_map, key_name, default_value, err_id_prefix)
            % Get an optional logical parameter value from map.
            % Input:
            %   param_map (containers.Map)
            %       parameter key-value map
            %   key_name (1 x 1 char/string)
            %       parameter key name
            %   default_value (logical/numeric/char/string, optional)
            %       default logical value if key not present;
            %       if omitted, returns false when key missing/empty
            %   err_id_prefix (1 x 1 char/string, optional)
            %       error id prefix (default: CParamMapUtils)
            % Output:
            %   value (logical)
            %       optional logical parameter value
            if nargin < 4 || isempty(err_id_prefix)
                err_id_prefix = 'CParamMapUtils';
            end
            if nargin < 3
                default_value = false;
            end

            raw_value = CParamMapUtils.getOptional(param_map, key_name, default_value, err_id_prefix);
            value = CParamMapUtils.toLogical(raw_value, key_name, err_id_prefix);
        end

        function values = parseQuotedList(raw_value, delimiter)
            % Parse quoted values from a delimiter-separated string.
            % Input:
            %   raw_value (1 x 1 char/string)
            %       source string containing quoted values
            %   delimiter (1 x 1 char/string, optional)
            %       segment delimiter (default: ';')
            % Output:
            %   values (cell)
            %       parsed values as 1 x N cell array
            if nargin < 2 || isempty(delimiter)
                delimiter = ';';
            end

            if isempty(raw_value)
                values = {};
                return;
            end

            if isstring(raw_value) && isscalar(raw_value)
                raw_value = char(raw_value);
            end
            if isstring(delimiter) && isscalar(delimiter)
                delimiter = char(delimiter);
            end

            if ~ischar(raw_value)
                error('CParamMapUtils:InvalidQuotedListType', ...
                    'Expected raw_value to be a char/string scalar.');
            end

            values = {};
            segments = strsplit(raw_value, delimiter);
            token_groups = cellfun(@(s) regexp(s, '"([^"]*)"', 'tokens'), ...
                segments, 'UniformOutput', false);
            for idx = 1:length(token_groups)
                if isempty(token_groups{idx})
                    continue;
                end
                values = [values, token_groups{idx}{1}]; %#ok<AGROW>
            end
        end

    end

    methods (Static, Access = private)
        function value = toNumber(raw_value, key_name, err_id_prefix)
            % Convert raw value to finite scalar double.
            if isnumeric(raw_value) || islogical(raw_value)
                if isempty(raw_value) || ~isscalar(raw_value)
                    error([err_id_prefix, ':InvalidNumericParam'], ...
                        'Param ''%s'' must be a numeric scalar.', key_name);
                end
                value = double(raw_value);
            elseif ischar(raw_value) || (isstring(raw_value) && isscalar(raw_value))
                value = str2double(raw_value);
            else
                error([err_id_prefix, ':InvalidNumericParamType'], ...
                    'Param ''%s'' must be numeric or numeric string.', key_name);
            end

            if isnan(value) || ~isfinite(value)
                error([err_id_prefix, ':InvalidNumericParamValue'], ...
                    'Param ''%s'' must be a finite numeric value.', key_name);
            end
        end

        function value = toLogical(raw_value, key_name, err_id_prefix)
            % Convert raw value to logical scalar.
            if islogical(raw_value)
                if isempty(raw_value) || ~isscalar(raw_value)
                    error([err_id_prefix, ':InvalidLogicalParam'], ...
                        'Param ''%s'' must be a logical scalar.', key_name);
                end
                value = logical(raw_value);
                return;
            end

            if isnumeric(raw_value)
                if isempty(raw_value) || ~isscalar(raw_value)
                    error([err_id_prefix, ':InvalidLogicalParam'], ...
                        'Param ''%s'' must be a numeric/logical scalar.', key_name);
                end
                value = raw_value ~= 0;
                return;
            end

            if isstring(raw_value) && isscalar(raw_value)
                raw_value = char(raw_value);
            end
            if ischar(raw_value)
                normalized = lower(strtrim(raw_value));
                if any(strcmp(normalized, {'1', 'true', 'yes', 'on'}))
                    value = true;
                    return;
                end
                if any(strcmp(normalized, {'0', 'false', 'no', 'off'}))
                    value = false;
                    return;
                end
            end

            error([err_id_prefix, ':InvalidLogicalParamValue'], ...
                'Param ''%s'' must be logical or a logical-like string (1/0/true/false/yes/no/on/off).', key_name);
        end
    end
end
