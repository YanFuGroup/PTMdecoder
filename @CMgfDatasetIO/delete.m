function delete(obj)
% Destructor for MGF dataset IO
% Input:
%   obj (CMgfDatasetIO)
%       dataset IO instance
if isempty(obj) || isempty(obj.m_mapFid)
    return;
end

obj.CloseAllFile();
end
