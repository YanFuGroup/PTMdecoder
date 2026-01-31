function check_files_existence(~, varargin)
% Check if the files exist.
% If yes, do nothing.
% If not, throw an error.
% Input:
%   varargin (1 x N cell)
%       list of file paths

for i = 1:length(varargin)
    if ~exist(varargin{i}, 'file')
        error('File %s does not exist.', varargin{i});
    end
end
end