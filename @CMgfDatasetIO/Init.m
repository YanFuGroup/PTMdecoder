function Init(obj,strDatasetFoldname)
% Initialize three data structures
obj.m_mapDatasetIdx=containers.Map();% Create an instance of containers.Map
% obj.m_mapSpecIdx=containers.Map();
obj.m_mapFid=containers.Map();% Create an instance of containers.Map
obj.m_strFoldname=strDatasetFoldname;
end
