function SetFidmap(obj)
% Build up fid for each data set
% Input:
%   obj (CMgfDatasetIO)
%       dataset IO instance
sDataset=dir(fullfile(obj.m_strFoldname,'*.mgf'));

for i=1:length(sDataset)
    obj.m_mapFid(sDataset(i).name)=obj.OpenFile(fullfile(obj.m_strFoldname,sDataset(i).name));
end

end