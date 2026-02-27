classdef CStructOptionUtils
    % Utilities for struct options

    methods (Static)
        function value = get(options, field_name, default_value)
            % Get an option value with a default fallback.
            % Input:
            %   options (struct)
            %       Options struct
            %   field_name (1 x 1 char/string)
            %       Field name to query
            %   default_value (any)
            %       Value to return if missing/empty
            % Output:
            %   value (any)
            %       Option value

            if nargin < 1 || isempty(options)
                options = struct();
            end
            if isfield(options, field_name)
                value = options.(field_name);
                if isempty(value)
                    value = default_value;
                end
            else
                value = default_value;
            end
        end
    end
end
