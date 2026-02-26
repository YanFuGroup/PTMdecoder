function path = resolveFilePath(file_dir, default_name, override_path)
% Resolve file path with override support.
% Input:
%   file_dir (1 x 1 char/string)
%       output directory
%   default_name (1 x 1 char/string)
%       default file name
%   override_path (1 x 1 char/string)
%       override file path
% Output:
%   path (1 x 1 char)
%       resolved output path

if nargin < 3
    override_path = '';
end

if ~isempty(override_path)
    path = char(string(override_path));
    return;
end

if isempty(file_dir)
    path = char(string(default_name));
else
    path = fullfile(char(string(file_dir)), char(string(default_name)));
end
end
