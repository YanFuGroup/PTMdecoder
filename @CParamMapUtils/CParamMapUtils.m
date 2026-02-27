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

        function value = getOptional(param_map, key_name, default_value)
            % Get an optional parameter value from map, use default if missing.
            % Input:
            %   param_map (containers.Map)
            %       parameter key-value map
            %   key_name (1 x 1 char/string)
            %       parameter key name
            %   default_value (any)
            %       default value if key not present
            % Output:
            %   value (any)
            %       parameter value or default value
            if ~isa(param_map, 'containers.Map')
                error('CParamMapUtils:InvalidParamMap', ...
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
            %   default_value (numeric/char/string)
            %       default numeric value if key not present
            %   err_id_prefix (1 x 1 char/string, optional)
            %       error id prefix (default: CParamMapUtils)
            % Output:
            %   value (double)
            %       optional numeric parameter value
            if nargin < 4 || isempty(err_id_prefix)
                err_id_prefix = 'CParamMapUtils';
            end
            raw_value = CParamMapUtils.getOptional(param_map, key_name, default_value);
            value = CParamMapUtils.toNumber(raw_value, key_name, err_id_prefix);
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
    end
end
