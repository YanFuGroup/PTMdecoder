classdef CWorkflowParamParser
    % Utility parser for workflow parameter files.

    methods (Static)
        task_param_map = parseFileToMap(param_file_path);
    end

    methods (Static, Access = private)
        function str_line = removeComments(str_line)
            % Remove comments from a line string.
            %   str_line (1 x 1 char/string)
            %       input line string
            idx = strfind(str_line, '#');
            if ~isempty(idx)
                str_line(idx(1):end) = [];
            end
        end
    end
end
