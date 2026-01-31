function CloseAllFile(obj)
% Close all file in the mgf dictionary
% Input:
%   obj (CMgfDatasetIO)
%       dataset IO instance
Fids=values(obj.m_mapFid);
for i=1:length(Fids)
    fclose(Fids{i});
end

end