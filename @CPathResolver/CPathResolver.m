classdef CPathResolver
    % Static path resolver utilities.
    % This class is intentionally stateless.

    methods (Static)
        % Resolve file path with override support.
        path = resolveFilePath(file_dir, default_name, override_path)

        % Ensure output directory exists.
        ensureDir(file_dir)
    end
end