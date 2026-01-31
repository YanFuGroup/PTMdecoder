function Init(obj,strDatasetFoldname)
% Initialize three data structures
% Input:
%   obj (CMgfDatasetIO)
%       dataset IO instance
%   strDatasetFoldname (1 x 1 char/string)
%       dataset folder path
obj.m_mapDatasetIdx=containers.Map();% Create an instance of containers.Map
% obj.m_mapSpecIdx=containers.Map();
obj.m_mapFid=containers.Map();% Create an instance of containers.Map
obj.m_strFoldname=strDatasetFoldname;
end
