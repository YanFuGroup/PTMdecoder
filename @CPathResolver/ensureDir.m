function ensureDir(file_dir)
% Ensure output directory exists.
% Input:
%   file_dir (1 x 1 char/string)
%       output directory path

if isempty(file_dir)
    return;
end
file_dir = char(string(file_dir));
if ~isfolder(file_dir)
    mkdir(file_dir);
end
end
